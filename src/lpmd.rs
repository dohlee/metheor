use rust_htslib::{bam, bam::Read};
use std::{collections::HashMap, fs};
use std::{io::Write, str, vec::Vec};

use crate::{bamutil, progressbar, readutil};

pub struct LPMDResult {
    header: bam::HeaderView,
    n_read: i32,
    n_valid_read: i32,
    n_concordant: i32,
    n_discordant: i32,
    pair2n_concordant: HashMap<(readutil::CpGPosition, readutil::CpGPosition), i32>,
    pair2n_discordant: HashMap<(readutil::CpGPosition, readutil::CpGPosition), i32>,
}

impl LPMDResult {
    fn new(header: bam::HeaderView) -> Self {
        let pair2n_concordant: HashMap<(readutil::CpGPosition, readutil::CpGPosition), i32> =
            HashMap::new();
        let pair2n_discordant: HashMap<(readutil::CpGPosition, readutil::CpGPosition), i32> =
            HashMap::new();

        Self {
            header,
            n_read: 0,
            n_valid_read: 0,
            n_concordant: 0,
            n_discordant: 0,
            pair2n_concordant,
            pair2n_discordant,
        }
    }

    fn inc_n_read(&mut self, i: i32) {
        self.n_read += i;
    }

    fn inc_n_valid_read(&mut self, i: i32) {
        self.n_valid_read += i;
    }

    fn inc_n_concordant(&mut self, i: i32) {
        self.n_concordant += i;
    }

    fn inc_n_discordant(&mut self, i: i32) {
        self.n_discordant += i;
    }

    fn compute_lpmd(&self) -> f32 {
        let lpmd: f32 =
            (self.n_discordant as f32) / ((self.n_concordant + self.n_discordant) as f32);
        lpmd
    }

    fn progress_string(&self) -> String {
        let lpmd = self.compute_lpmd();

        format!(
            "Processed {} reads, found {} valid reads. LPMD={:.4} ({}/{})",
            self.n_read,
            self.n_valid_read,
            lpmd,
            self.n_discordant,
            self.n_concordant + self.n_discordant
        )
    }

    fn add_pair_concordance(
        &mut self,
        pos1: &readutil::CpGPosition,
        pos2: &readutil::CpGPosition,
        concordance: &readutil::ReadConcordanceState,
    ) {
        let n_concordant = self.pair2n_concordant.entry((*pos1, *pos2)).or_insert(0);
        let n_discordant = self.pair2n_discordant.entry((*pos1, *pos2)).or_insert(0);

        match concordance {
            readutil::ReadConcordanceState::Concordant => {
                *n_concordant += 1;
            }
            readutil::ReadConcordanceState::Discordant => {
                *n_discordant += 1;
            }
        }
    }

    fn print_pair_statistics(&self, output: &str) {
        let mut pairs: Vec<&(readutil::CpGPosition, readutil::CpGPosition)> = self
            .pair2n_concordant
            .keys()
            .collect::<Vec<&(readutil::CpGPosition, readutil::CpGPosition)>>();
        pairs.sort();

        let mut out = fs::OpenOptions::new()
            .create(true)
            .read(true)
            .write(true)
            .truncate(true)
            .open(output)
            .unwrap();

        writeln!(out, "chrom\tcpg1\tcpg2\tlpmd\tn_concordant\tn_discordant")
            .expect("Error writing to output file.");

        for (cpg1, cpg2) in pairs {
            let k = (*cpg1, *cpg2);
            let n_concordant = self.pair2n_concordant[&k];
            let n_discordant = self.pair2n_discordant[&k];
            let lpmd = (n_discordant as f32) / (n_concordant as f32 + n_discordant as f32);

            let chrom = bamutil::tid2chrom(cpg1.tid, &self.header);

            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}",
                chrom, cpg1.pos, cpg2.pos, lpmd, n_concordant, n_discordant
            )
            .expect("Error writing to output file.");
        }
    }
}

pub fn compute(
    input: &str,
    output: &str,
    min_distance: i32,
    max_distance: i32,
    min_qual: u8,
    cpg_set: &Option<String>,
    pairs: &Option<String>,
) {
    let result = compute_helper(input, min_distance, max_distance, min_qual, cpg_set);
    let lpmd = result.compute_lpmd();

    let mut out = fs::OpenOptions::new()
        .create(true)
        .read(true)
        .write(true)
        .truncate(true)
        .open(output)
        .unwrap();

    writeln!(out, "name\tlpmd").expect("Error writing to output file.");

    writeln!(out, "{}\t{}", input, lpmd).expect("Error writing to output file.");

    if let Some(f) = pairs {
        result.print_pair_statistics(f);
    }
}

fn compute_helper(
    input: &str,
    min_distance: i32,
    max_distance: i32,
    min_qual: u8,
    cpg_set: &Option<String>,
) -> LPMDResult {
    eprintln!(
        "Computing subset-LPMD with parameters input={}, min_distance={}, max_distance={}",
        input, min_distance, max_distance
    );
    let mut reader = bamutil::get_reader(input);
    let header = bamutil::get_header(&reader);

    eprint!("Processing target CpG set... ");
    let target_cpgs = &readutil::get_target_cpgs(cpg_set, &header);

    let mut res = LPMDResult::new(header);
    let bar = progressbar::ProgressBar::new();

    // Iterate over reads and compute LPMD.
    for r in reader.records().map(|r| r.unwrap()) {
        res.inc_n_read(1);
        if r.mapq() < min_qual {
            continue;
        }

        let mut br = readutil::BismarkRead::new(&r);
        if let Some(target_cpgs) = target_cpgs {
            br.filter_isin(target_cpgs);
        }

        let (c, d, pair2concordance) =
            br.compute_pairwise_cpg_concordance_discordance(min_distance, max_distance);

        res.inc_n_valid_read(1);
        res.inc_n_concordant(c);
        res.inc_n_discordant(d);
        for (cpg1, cpg2, concordance) in &pair2concordance {
            res.add_pair_concordance(cpg1, cpg2, concordance);
        }

        if res.n_read % 10000 == 0 {
            bar.update_lpmd(res.progress_string());
        }
    }

    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test1() {
        let input = "tests/test1.bam";

        let min_distance = 2;
        let max_distance = 16;
        let min_qual = 10;
        let cpg_set = None;

        let result = compute_helper(input, min_distance, max_distance, min_qual, &cpg_set);

        assert_eq!(result.compute_lpmd(), 0.5);
    }
    #[test]
    fn test2() {
        let input = "tests/test2.bam";
        let min_distance = 2;
        let max_distance = 16;
        let min_qual = 10;
        let cpg_set = None;

        let result = compute_helper(input, min_distance, max_distance, min_qual, &cpg_set);

        assert_eq!(result.compute_lpmd(), 0.0);
    }
    #[test]
    fn test3() {
        let input = "tests/test3.bam";
        let min_distance = 2;
        let max_distance = 16;
        let min_qual = 10;
        let cpg_set = None;

        let result = compute_helper(input, min_distance, max_distance, min_qual, &cpg_set);

        assert_eq!(result.compute_lpmd(), 0.0);
    }
    #[test]
    fn test4() {
        let input = "tests/test4.bam";
        let min_distance = 2;
        let max_distance = 16;
        let min_qual = 10;
        let cpg_set = None;

        let result = compute_helper(input, min_distance, max_distance, min_qual, &cpg_set);

        assert_eq!(result.compute_lpmd(), 0.5);
    }
    #[test]
    fn test5() {
        // No reads pass quality cutoff.
        let input = "tests/test5.bam";
        let min_distance = 2;
        let max_distance = 16;
        let min_qual = 10;
        let cpg_set = None;

        let result = compute_helper(input, min_distance, max_distance, min_qual, &cpg_set);

        assert!(result.compute_lpmd().is_nan());
    }
}
