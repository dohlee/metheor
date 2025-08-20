use rust_htslib::bam::Read;
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::fs;
use std::io::Write;
use std::str;

use crate::{bamutil, progressbar, readutil};

#[derive(Eq)]
struct PDRResult {
    pos: readutil::CpGPosition,
    n_concordant: u32,
    n_discordant: u32,
}

impl PDRResult {
    fn new(pos: readutil::CpGPosition) -> Self {
        Self {
            pos,
            n_concordant: 0,
            n_discordant: 0,
        }
    }

    fn inc_concordant(&mut self) {
        self.n_concordant += 1;
    }

    fn inc_discordant(&mut self) {
        self.n_discordant += 1;
    }

    fn get_n_concordant(&self) -> u32 {
        self.n_concordant
    }

    fn get_n_discordant(&self) -> u32 {
        self.n_discordant
    }

    fn get_coverage(&self) -> u32 {
        self.n_concordant + self.n_discordant
    }

    fn compute_pdr(&self) -> f32 {
        (self.n_discordant as f32) / (self.n_concordant as f32 + self.n_discordant as f32)
    }
}

impl fmt::Display for PDRResult {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let pdr =
            (self.n_discordant as f32) / (self.n_concordant as f32 + self.n_discordant as f32);
        write!(
            f,
            "{}\t{}\t{}\t{}",
            self.pos, pdr, self.n_concordant, self.n_discordant
        )
    }
}

impl PartialEq for PDRResult {
    fn eq(&self, other: &Self) -> bool {
        self.pos == other.pos
    }
}

impl PartialOrd for PDRResult {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PDRResult {
    fn cmp(&self, other: &Self) -> Ordering {
        self.pos.cmp(&other.pos)
    }
}

pub fn compute(
    input: &str,
    output: &str,
    min_depth: u32,
    min_cpgs: usize,
    min_qual: u8,
    cpg_set: &Option<String>,
) {
    let reader = bamutil::get_reader(input);
    let header = bamutil::get_header(&reader);

    let result = compute_helper(input, min_depth, min_cpgs, min_qual, cpg_set);

    let mut out = fs::OpenOptions::new()
        .create(true)
        .read(true)
        .write(true)
        .truncate(true)
        .open(output)
        .unwrap();
    for (cpg, (pdr, n_concordant, n_discordant)) in result.iter() {
        let chrom = bamutil::tid2chrom(cpg.tid, &header);

        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}",
            chrom,
            cpg.pos,
            cpg.pos + 2,
            pdr,
            n_concordant,
            n_discordant
        )
        .expect("Error writing to output file.");
    }
}

pub fn compute_helper(
    input: &str,
    min_depth: u32,
    min_cpgs: usize,
    min_qual: u8,
    cpg_set: &Option<String>,
) -> BTreeMap<readutil::CpGPosition, (f32, u32, u32)> {
    let mut reader = bamutil::get_reader(input);
    let header = bamutil::get_header(&reader);

    let target_cpgs = &readutil::get_target_cpgs(cpg_set, &header);

    let mut cpg2reads: HashMap<readutil::CpGPosition, PDRResult> = HashMap::new();

    let mut readcount = 0;
    let mut valid_readcount = 0;

    let mut result: BTreeMap<readutil::CpGPosition, (f32, u32, u32)> = BTreeMap::new();
    let bar = progressbar::ProgressBar::new();

    for r in reader.records().map(|r| r.unwrap()) {
        let mut br = readutil::BismarkRead::new(&r);

        if let Some(target_cpgs) = target_cpgs {
            br.filter_isin(target_cpgs); // cpg_set is specified
        }

        readcount += 1;
        if br.get_num_cpgs() < min_cpgs {
            continue;
        }
        if r.mapq() < min_qual {
            continue;
        } // Read filtering: Minimum quality should be >= min_qual.

        let mut cpg_positions = br.get_cpg_positions();
        if cpg_positions.is_empty() {
            continue;
        } // Read filtering: Ignore reads without CpGs.

        if let Some(first_cpg_position) = br.get_first_cpg_position() {
            cpg2reads.retain(|&cpg, reads| {
                // if cpg < first_cpg_position {
                if cpg.is_before(&first_cpg_position, 150) {
                    if reads.get_coverage() >= min_depth {
                        result.insert(
                            cpg,
                            (
                                reads.compute_pdr(),
                                reads.get_n_concordant(),
                                reads.get_n_discordant(),
                            ),
                        );
                    }
                    false
                } else {
                    true
                }
            }); // Finalize and compute metric for the CpGs before the first CpG in this read.
        }

        for cpg_position in cpg_positions.iter_mut() {
            let r = cpg2reads
                .entry(*cpg_position)
                .or_insert(PDRResult::new(*cpg_position));

            let concordance_state = br.get_concordance_state();

            match concordance_state {
                readutil::ReadConcordanceState::Concordant => r.inc_concordant(),
                readutil::ReadConcordanceState::Discordant => r.inc_discordant(),
            }
        }

        valid_readcount += 1;
        if readcount % 10000 == 0 {
            bar.update(readcount, valid_readcount)
        };
    }

    for (&cpg, reads) in cpg2reads.iter() {
        if reads.get_coverage() >= min_depth {
            result.insert(
                cpg,
                (
                    reads.compute_pdr(),
                    reads.get_n_concordant(),
                    reads.get_n_discordant(),
                ),
            );
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test1() {
        let input = "tests/test1.bam";
        let min_depth = 0;
        let min_cpgs = 0;
        let min_qual = 10;
        let cpg_set = None;

        let target_pdrs = [14.0 / 16.0; 4];
        let target_n_concordant = [2; 4];
        let target_n_discordant = [14; 4];

        let result = compute_helper(input, min_depth, min_cpgs, min_qual, &cpg_set);

        assert_eq!(result.len(), 4);
        for (i, (_, (pdr, n_concordant, n_discordant))) in result.iter().enumerate() {
            assert_eq!(*pdr, target_pdrs[i]);
            assert_eq!(*n_concordant, target_n_concordant[i]);
            assert_eq!(*n_discordant, target_n_discordant[i]);
        }
    }

    #[test]
    fn test2() {
        let input = "tests/test2.bam";
        let min_depth = 0;
        let min_cpgs = 0;
        let min_qual = 10;
        let cpg_set = None;

        let target_pdrs = [0.0; 4];
        let target_n_concordant = [16; 4];
        let target_n_discordant = [0; 4];

        let result = compute_helper(input, min_depth, min_cpgs, min_qual, &cpg_set);

        assert_eq!(result.len(), 4);
        for (i, (_, (pdr, n_concordant, n_discordant))) in result.iter().enumerate() {
            assert_eq!(*pdr, target_pdrs[i]);
            assert_eq!(*n_concordant, target_n_concordant[i]);
            assert_eq!(*n_discordant, target_n_discordant[i]);
        }
    }

    #[test]
    fn test3() {
        let input = "tests/test3.bam";

        let min_depth = 0;
        let min_cpgs = 0;
        let min_qual = 10;
        let cpg_set = None;

        let target_pdrs = [0.0; 4];
        let target_n_concordant = [2; 4];
        let target_n_discordant = [0; 4];

        let result = compute_helper(input, min_depth, min_cpgs, min_qual, &cpg_set);

        assert_eq!(result.len(), 4);
        for (i, (_, (pdr, n_concordant, n_discordant))) in result.iter().enumerate() {
            assert_eq!(*pdr, target_pdrs[i]);
            assert_eq!(*n_concordant, target_n_concordant[i]);
            assert_eq!(*n_discordant, target_n_discordant[i]);
        }
    }

    #[test]
    fn test4() {
        // Two read groups, each has 16 reads, aligned to two distinct position.
        // Reads have four CpGs, all the 16 methylation patterns are present. (Max entropy)
        let input = "tests/test4.bam";

        let min_depth = 0;
        let min_cpgs = 0;
        let min_qual = 10;
        let cpg_set = None;

        let target_pdrs = [14.0 / 16.0; 8]; // Two concordant patterns (0000, 1111)
        let target_n_concordant = [2; 8];
        let target_n_discordant = [14; 8];

        let result = compute_helper(input, min_depth, min_cpgs, min_qual, &cpg_set);

        assert_eq!(result.len(), 8);
        for (i, (_, (pdr, n_concordant, n_discordant))) in result.iter().enumerate() {
            assert_eq!(*pdr, target_pdrs[i]);
            assert_eq!(*n_concordant, target_n_concordant[i]);
            assert_eq!(*n_discordant, target_n_discordant[i]);
        }
    }

    #[test]
    fn test5() {
        // No reads pass quality cutoff.
        let input = "tests/test5.bam";

        let min_depth = 0;
        let min_cpgs = 0;
        let min_qual = 10;
        let cpg_set = None;

        let result = compute_helper(input, min_depth, min_cpgs, min_qual, &cpg_set);

        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test6_mincpg_low() {
        // Similar to test4, but reads have only one CpGs.
        // There are 2 unique positions of CpGs.
        // Here we set min_cpgs to 1, so every read passes the filter.
        let input = "tests/test6.bam";

        let min_depth = 0;
        let min_cpgs = 1;
        let min_qual = 10;
        let cpg_set = None;

        let target_pdrs = [0.0; 2];
        let target_n_concordant = [16; 2];
        let target_n_discordant = [0; 2];

        let result = compute_helper(input, min_depth, min_cpgs, min_qual, &cpg_set);

        assert_eq!(result.len(), 2);
        for (i, (_, (pdr, n_concordant, n_discordant))) in result.iter().enumerate() {
            assert_eq!(*pdr, target_pdrs[i]);
            assert_eq!(*n_concordant, target_n_concordant[i]);
            assert_eq!(*n_discordant, target_n_discordant[i]);
        }
    }

    #[test]
    fn test6_mincpg_high() {
        // Similar to test4, but reads have only one CpGs.
        // There are 2 unique positions of CpGs.
        // Here we set min_cpgs to 2, so no reads pass the filter.
        let input = "tests/test6.bam";

        let min_depth = 0;
        let min_cpgs = 2;
        let min_qual = 10;
        let cpg_set = None;

        let result = compute_helper(input, min_depth, min_cpgs, min_qual, &cpg_set);
        assert_eq!(result.len(), 0); // No CpGs participate in the PDR calculation.
    }
}
