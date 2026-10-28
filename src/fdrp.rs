use rust_htslib::bam::Read;
use std::collections::BTreeMap;
use std::fs;
use std::io::Write;

use crate::{bamutil, progressbar, readutil};

use crate::pairwise::AssociatedReads;

pub fn compute(
    input: &str,
    output: &str,
    min_qual: u8,
    min_depth: usize,
    max_depth: usize,
    min_overlap: i32,
    cpg_set: &Option<String>,
) {
    let result = compute_helper(input, min_qual, min_depth, max_depth, min_overlap, cpg_set);

    let reader = bamutil::get_reader(input);
    let header = bamutil::get_header(&reader);

    let mut out = fs::OpenOptions::new()
        .create(true)
        .read(true)
        .write(true)
        .truncate(true)
        .open(output)
        .unwrap();
    for (cpg, fdrp) in result.iter() {
        let chrom = bamutil::tid2chrom(cpg.tid, &header);
        writeln!(out, "{}\t{}\t{}\t{}", chrom, cpg.pos, cpg.pos + 2, fdrp)
            .expect("Error writing to output file.");
    }
}

fn compute_helper(
    input: &str,
    min_qual: u8,
    min_depth: usize,
    max_depth: usize,
    min_overlap: i32,
    cpg_set: &Option<String>,
) -> BTreeMap<readutil::CpGPosition, f32> {
    let mut reader = bamutil::get_reader(input);
    let header = bamutil::get_header(&reader);

    let mut readcount = 0;
    let mut valid_readcount = 0;

    let target_cpgs = &readutil::get_target_cpgs(cpg_set, &header);

    let bar = progressbar::ProgressBar::new();

    let mut cpg2reads: BTreeMap<readutil::CpGPosition, AssociatedReads> = BTreeMap::new();
    let mut result: BTreeMap<readutil::CpGPosition, f32> = BTreeMap::new();

    for r in reader.records().map(|r| r.unwrap()) {
        let mut br = readutil::BismarkRead::new(&r);

        if let Some(target_cpgs) = target_cpgs {
            br.filter_isin(target_cpgs);
        }

        readcount += 1;
        if r.mapq() < min_qual {
            continue;
        }
        if br.get_num_cpgs() == 0 {
            continue;
        }

        // Reverse-strand CpGs can begin one base before the next alignment.
        let read_start = readutil::CpGPosition::new(r.tid(), r.pos() as i32);
        cpg2reads.retain(|&cpg, reads| {
            if cpg.is_before(&read_start, 1) {
                if reads.get_num_reads() >= min_depth {
                    result.insert(cpg, reads.compute_fdrp(min_overlap));
                }
                false
            } else {
                true
            }
        });

        for cpg_position in br.get_cpg_positions().iter() {
            let r = cpg2reads
                .entry(*cpg_position)
                .or_insert(AssociatedReads::new(max_depth));

            r.add_read(&br);
        }
        valid_readcount += 1;
        if readcount % 10000 == 0 {
            bar.update(readcount, valid_readcount)
        };
    }

    // Flush remaining CpGs.
    for (cpg, reads) in cpg2reads.iter_mut() {
        if reads.get_num_reads() >= min_depth {
            result.insert(*cpg, reads.compute_fdrp(min_overlap));
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
        let min_qual = 0;
        let min_depth = 2;
        let max_depth = 40;
        let min_overlap = 4;
        let cpg_set = None;

        let cpg_positions = [0, 2, 4, 6];

        let result = compute_helper(input, min_qual, min_depth, max_depth, min_overlap, &cpg_set);
        for (i, (cpg, fdrp)) in result.iter().enumerate() {
            assert_eq!(cpg.pos, cpg_positions[i]);
            assert_eq!(*fdrp, 1.0);
        }
    }
    #[test]
    fn test2() {
        let input = "tests/test2.bam";
        let min_qual = 0;
        let min_depth = 2;
        let max_depth = 40;
        let min_overlap = 4;
        let cpg_set = None;

        let cpg_positions = [0, 2, 4, 6];

        let result = compute_helper(input, min_qual, min_depth, max_depth, min_overlap, &cpg_set);
        for (i, (cpg, fdrp)) in result.iter().enumerate() {
            assert_eq!(cpg.pos, cpg_positions[i]);
            assert!((*fdrp - (1.0 - 56.0 / 120.0)).abs() < 1e-4); // Approximately same.
        }
    }
    #[test]
    fn test3() {
        let input = "tests/test3.bam";
        let min_qual = 1;
        let min_depth = 2;
        let max_depth = 40;
        let min_overlap = 4;
        let cpg_set = None;

        let cpg_positions = [0, 2, 4, 6];

        let result = compute_helper(input, min_qual, min_depth, max_depth, min_overlap, &cpg_set);
        for (i, (cpg, fdrp)) in result.iter().enumerate() {
            assert_eq!(cpg.pos, cpg_positions[i]);
            assert_eq!(*fdrp, 1.0);
        }
    }
    #[test]
    fn test4() {
        let input = "tests/test4.bam";
        let min_qual = 1;
        let min_depth = 2;
        let max_depth = 40;
        let min_overlap = 4;
        let cpg_set = None;

        let cpg_positions = [0, 2, 4, 6, 13, 15, 17, 19];

        let result = compute_helper(input, min_qual, min_depth, max_depth, min_overlap, &cpg_set);
        for (i, (cpg, fdrp)) in result.iter().enumerate() {
            assert_eq!(cpg.pos, cpg_positions[i]);
            assert_eq!(*fdrp, 1.0);
        }
    }
    #[test]
    fn test5() {
        // No reads pass quality cutoff.
        let input = "tests/test5.bam";
        let min_qual = 1;
        let min_depth = 2;
        let max_depth = 40;
        let min_overlap = 4;
        let cpg_set = None;

        let result = compute_helper(input, min_qual, min_depth, max_depth, min_overlap, &cpg_set);
        assert_eq!(result.len(), 0);
    }
}
