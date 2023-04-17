use rust_htslib::bam::Read;
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::io::Write;
use std::str;
use std::vec::Vec;

use crate::{bamutil, progressbar, readutil};

#[derive(Eq)]
struct AssociatedReads {
    pos: readutil::CpGPosition,
    stretch_info: HashMap<i32, i32>,
    num_cpgs: Vec<i32>,
    max_num_cpgs: usize,
}

impl AssociatedReads {
    fn new(pos: readutil::CpGPosition) -> Self {
        let stretch_info: HashMap<i32, i32> = HashMap::new();
        let num_cpgs: Vec<i32> = Vec::new();
        let max_num_cpgs = 0;
        Self {
            pos,
            stretch_info,
            num_cpgs,
            max_num_cpgs,
        }
    }

    fn get_coverage(&self) -> u32 {
        self.num_cpgs.len() as u32
    }

    fn add_stretch_info(&mut self, stretch_info: HashMap<i32, i32>) {
        for (l, count) in stretch_info.iter() {
            let curr_count = self.stretch_info.entry(*l).or_insert(0);
            *curr_count += count;
        }
    }

    fn compute_mhl(&self) -> f32 {
        let mut mhl = 0.0;
        let mut l_sum = 0.0;
        for l in 1..self.max_num_cpgs + 1 {
            l_sum += l as f32;
        }

        for (&l, count) in self.stretch_info.iter() {
            let dom = *count as f32;

            let mut denom = 0.0;
            for num_cpg in self.num_cpgs.iter() {
                if num_cpg >= &l {
                    denom += (num_cpg - l + 1) as f32;
                }
            }

            assert!(
                denom > 0.0,
                "denom <= 0!, max_num_cpgs={}, num_cpgs={:?}, l={}",
                self.max_num_cpgs,
                self.num_cpgs,
                l
            );

            mhl += (l as f32 * dom) / denom;
        }

        mhl /= l_sum;
        mhl
    }

    fn add_num_cpgs(&mut self, num_cpgs: usize) {
        self.num_cpgs.push(num_cpgs as i32);
        if num_cpgs >= self.max_num_cpgs {
            self.max_num_cpgs = num_cpgs;
        }
    }
}

impl PartialEq for AssociatedReads {
    fn eq(&self, other: &Self) -> bool {
        self.pos == other.pos
    }
}

impl PartialOrd for AssociatedReads {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(&other))
    }
}

impl Ord for AssociatedReads {
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
    let reader = bamutil::get_reader(&input);
    let header = bamutil::get_header(&reader);

    let result = compute_helper(input, min_depth, min_cpgs, min_qual, cpg_set);

    let mut out = fs::OpenOptions::new()
        .create(true)
        .read(true)
        .write(true)
        .truncate(true)
        .open(output)
        .unwrap();

    for (cpg, mhl) in result.iter() {
        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            bamutil::tid2chrom(cpg.tid, &header),
            cpg.pos,
            cpg.pos + 2,
            mhl
        )
        .ok()
        .expect("Error writing to output file.");
    }
}

pub fn compute_helper(
    input: &str,
    min_depth: u32,
    min_cpgs: usize,
    min_qual: u8,
    cpg_set: &Option<String>,
) -> BTreeMap<readutil::CpGPosition, f32> {
    let mut reader = bamutil::get_reader(&input);
    let header = bamutil::get_header(&reader);

    let target_cpgs = &readutil::get_target_cpgs(cpg_set, &header);

    let mut cpg2reads: HashMap<readutil::CpGPosition, AssociatedReads> = HashMap::new();
    let mut result: BTreeMap<readutil::CpGPosition, f32> = BTreeMap::new();

    let mut readcount = 0;
    let mut valid_readcount = 0;

    let bar = progressbar::ProgressBar::new();

    for r in reader.records().map(|r| r.unwrap()) {
        let mut br = readutil::BismarkRead::new(&r);

        match target_cpgs {
            Some(target_cpgs) => br.filter_isin(target_cpgs),
            None => {}
        }

        match br.get_first_cpg_position() {
            Some(first_cpg_position) => {
                cpg2reads.retain(|&cpg, reads| {
                    let retain = {
                        if cpg < first_cpg_position {
                            if reads.get_coverage() >= min_depth {
                                result.insert(cpg, reads.compute_mhl());
                            }
                            false
                        } else {
                            true
                        }
                    };
                    retain
                }); // Finalize and compute metric for the CpGs before the first CpG in this read.
            }
            None => {}
        }

        readcount += 1;
        if r.mapq() < min_qual {
            continue;
        } // Read filtering: Minimum quality should be >= min_qual.

        let mut cpg_positions = br.get_cpg_positions();
        if br.get_num_cpgs() < min_cpgs {
            continue;
        } // Read filtering: Ignore reads with few CpGs.

        for cpg_position in cpg_positions.iter_mut() {
            let r = cpg2reads
                .entry(*cpg_position)
                .or_insert(AssociatedReads::new(*cpg_position));

            r.add_num_cpgs(br.get_num_cpgs());
            r.add_stretch_info(br.get_stretch_info());
        }

        valid_readcount += 1;
        if readcount % 10000 == 0 {
            bar.update(readcount, valid_readcount)
        };
    }

    // Flush remaining CpGs.
    for (&cpg, reads) in cpg2reads.iter_mut() {
        if reads.get_coverage() >= min_depth {
            result.insert(cpg, reads.compute_mhl());
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::super::bamutil;
    use super::*;

    fn startup(input: &str) -> HashMap<readutil::CpGPosition, AssociatedReads> {
        let mut reader = bamutil::get_reader(&input);
        // let header = bamutil::get_header(&reader);

        let min_qual = 10;

        let mut cpg2reads: HashMap<readutil::CpGPosition, AssociatedReads> = HashMap::new();

        for r in reader.records().map(|r| r.unwrap()) {
            let br = readutil::BismarkRead::new(&r);
            if r.mapq() < min_qual {
                continue;
            } // Read filtering: Minimum quality should be >= min_qual.

            let mut cpg_positions = br.get_cpg_positions();

            for cpg_position in cpg_positions.iter_mut() {
                let r = cpg2reads
                    .entry(*cpg_position)
                    .or_insert(AssociatedReads::new(*cpg_position));

                r.add_num_cpgs(br.get_num_cpgs());
                r.add_stretch_info(br.get_stretch_info());
            }
        }

        cpg2reads
    }
    #[test]
    fn test_stretch_info() {
        let input = "tests/test1.bam";
        let cpg2reads = startup(input);

        for (_, reads) in cpg2reads.iter() {
            assert_eq!(reads.compute_mhl(), 0.1625);
        }
    }
    #[test]
    fn test1() {
        let input = "tests/test1.bam";
        let cpg2reads = startup(input);

        for (_, reads) in cpg2reads.iter() {
            assert_eq!(reads.compute_mhl(), 0.1625);
        }
    }
    #[test]
    fn test2() {
        let input = "tests/test2.bam";
        let cpg2reads = startup(input);

        for (_, reads) in cpg2reads.iter() {
            assert_eq!(reads.compute_mhl(), 0.5);
        }
    }
    #[test]
    fn test3() {
        let input = "tests/test3.bam";
        let cpg2reads = startup(input);

        for (_, reads) in cpg2reads.iter() {
            assert_eq!(reads.compute_mhl(), 0.5);
        }
    }
    #[test]
    fn test4() {
        let input = "tests/test4.bam";
        let cpg2reads = startup(input);

        for (_, reads) in cpg2reads.iter() {
            assert_eq!(reads.compute_mhl(), 0.1625);
        }
    }
    #[test]
    fn test5() {
        // No reads pass quality cutoff.
        let input = "tests/test5.bam";
        let cpg2reads = startup(input);

        assert_eq!(cpg2reads.len(), 0);
    }
}
