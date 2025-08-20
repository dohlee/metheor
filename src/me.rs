use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::str;

use crate::{bamutil, progressbar, readutil};

pub struct QuartetStat {
    pos1: readutil::CpGPosition,
    pos2: readutil::CpGPosition,
    pos3: readutil::CpGPosition,
    pos4: readutil::CpGPosition,
    quartet_pattern_counts: [u32; 16],
}

impl QuartetStat {
    fn new(q: readutil::Quartet) -> Self {
        let pos1 = q.pos1;
        let pos2 = q.pos2;
        let pos3 = q.pos3;
        let pos4 = q.pos4;

        let quartet_pattern_counts = [0; 16];
        Self {
            pos1,
            pos2,
            pos3,
            pos4,
            quartet_pattern_counts,
        }
    }

    fn get_read_depth(&self) -> u32 {
        self.quartet_pattern_counts.iter().sum()
    }

    fn add_quartet_pattern(&mut self, p: readutil::QuartetPattern) {
        self.quartet_pattern_counts[p] += 1;
    }

    fn compute_me(&self) -> f32 {
        let mut me: f32 = 0.0;

        let total: u32 = self.quartet_pattern_counts.iter().sum();
        for count in self.quartet_pattern_counts.iter() {
            let p: f32 = (*count as f32) / (total as f32);
            if *count > 0 {
                me += p * p.log2();
            }
        }
        me *= -0.25;

        me
    }

    fn to_bedgraph_field(&self, header: &bam::HeaderView) -> String {
        let chrom = bamutil::tid2chrom(self.pos1.tid, header);
        let me = self.compute_me();

        format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            chrom, self.pos1.pos, self.pos2.pos, self.pos3.pos, self.pos4.pos, me
        )
    }
}

pub fn compute(input: &str, output: &str, min_depth: u32, min_qual: u8, cpg_set: &Option<String>) {
    let reader = bamutil::get_reader(input);
    let header = bamutil::get_header(&reader);

    let result = compute_helper(input, min_qual, cpg_set);

    let mut out = fs::OpenOptions::new()
        .create(true)
        .read(true)
        .write(true)
        .truncate(true)
        .open(output)
        .unwrap();
    for stat in result.values() {
        if stat.get_read_depth() < min_depth {
            continue;
        }
        writeln!(out, "{}", stat.to_bedgraph_field(&header))
            .expect("Error writing to output file.");
    }
}

pub fn compute_helper(
    input: &str,
    min_qual: u8,
    cpg_set: &Option<String>,
) -> HashMap<readutil::Quartet, QuartetStat> {
    let mut reader = bamutil::get_reader(input);
    let header = bamutil::get_header(&reader);

    let target_cpgs = &readutil::get_target_cpgs(cpg_set, &header);
    let mut quartet2stat: HashMap<readutil::Quartet, QuartetStat> = HashMap::new();

    let mut readcount = 0;
    let mut valid_readcount = 0;

    let bar = progressbar::ProgressBar::new();

    for r in reader.records().map(|r| r.unwrap()) {
        let mut br = readutil::BismarkRead::new(&r);

        if let Some(target_cpgs) = target_cpgs {
            br.filter_isin(target_cpgs);
        }

        readcount += 1;

        if r.mapq() < min_qual {
            continue;
        }
        valid_readcount += 1;

        let (quartets, patterns) = br.get_cpg_quartets_and_patterns();
        for (q, p) in quartets.iter().zip(patterns.iter()) {
            let stat = quartet2stat.entry(*q).or_insert(QuartetStat::new(*q));

            stat.add_quartet_pattern(*p);
        }

        if readcount % 10000 == 0 {
            bar.update(readcount, valid_readcount)
        };
    }
    quartet2stat
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test1() {
        let input = "tests/test1.bam";
        let min_qual = 10;
        let cpg_set = None;

        let quartet2stat = compute_helper(input, min_qual, &cpg_set);

        assert_eq!(quartet2stat.len(), 1);

        for (_, reads) in quartet2stat.iter() {
            assert_eq!(reads.get_read_depth(), 16);
            assert_eq!(reads.compute_me(), 1.0);
        }
    }

    #[test]
    fn test2() {
        let input = "tests/test2.bam";
        let min_qual = 10;
        let cpg_set = None;

        let quartet2stat = compute_helper(input, min_qual, &cpg_set);

        assert_eq!(quartet2stat.len(), 1);

        for (_, reads) in quartet2stat.iter() {
            assert_eq!(reads.compute_me(), 0.25);
        }
    }
    #[test]
    fn test3() {
        let input = "tests/test3.bam";
        let min_qual = 10;
        let cpg_set = None;

        let quartet2stat = compute_helper(input, min_qual, &cpg_set);

        assert_eq!(quartet2stat.len(), 1);

        for (_, reads) in quartet2stat.iter() {
            assert_eq!(reads.compute_me(), 0.25);
        }
    }
    #[test]
    fn test4() {
        let input = "tests/test4.bam";
        let min_qual = 10;
        let cpg_set = None;

        let quartet2stat = compute_helper(input, min_qual, &cpg_set);

        assert_eq!(quartet2stat.len(), 2);

        for (_, reads) in quartet2stat.iter() {
            assert_eq!(reads.compute_me(), 1.0);
        }
    }
    #[test]
    fn test5() {
        // No reads pass quality cutoff.
        let input = "tests/test5.bam";

        let min_qual = 10;
        let cpg_set = None;

        let quartet2stat = compute_helper(input, min_qual, &cpg_set);

        assert_eq!(quartet2stat.len(), 0);
    }
}
