use itertools::Itertools;
use rand::Rng;
use rust_htslib::bam::Read;
use std::collections::BTreeMap;
use std::fs;
use std::io::Write;

use crate::{bamutil, progressbar, readutil};

const MAX_READ_LEN: i32 = 201;

struct AssociatedReads {
    // Use compact representation of reads.
    // Position "MAX_READ_LEN" represents this CpG, and positions of other CpGs are
    // determined according to the fixed position "MAX_READ_LEN".
    // Each position in the array is filled with three-bit representation of reads.
    // 000 (0 in decimal) : read does not span this potiion.
    // 001 (1 in decimal) : read covers this position, but the base at this position is not C of CpG.
    // 011 (3 in decimal) : read covers this position, but CpG at this position is not methylated.
    // 111 (7 in decimal) : read covers this position, and CpG at this position is methylated.
    pos: readutil::CpGPosition,
    reads: Vec<[u8; (MAX_READ_LEN * 2 + 1) as usize]>,
    num_total_read: i32,
    num_sampled_read: i32,
    max_depth: usize,
}

impl AssociatedReads {
    fn new(pos: readutil::CpGPosition, max_depth: usize) -> Self {
        let reads: Vec<[u8; (MAX_READ_LEN * 2 + 1) as usize]> = Vec::new();
        let num_total_read = 0;
        let num_sampled_read = 0;

        Self {
            pos,
            reads,
            num_total_read,
            num_sampled_read,
            max_depth,
        }
    }

    fn get_relative_position(&self, other_pos: readutil::CpGPosition) -> usize {
        (MAX_READ_LEN + (other_pos.pos - self.pos.pos)) as usize
    }

    fn get_num_reads(&self) -> usize {
        self.num_sampled_read as usize
    }

    fn add_read(&mut self, br: &readutil::BismarkRead) {
        let mut new_read: [u8; (MAX_READ_LEN * 2 + 1) as usize] =
            [0; (MAX_READ_LEN * 2 + 1) as usize];

        let start_relative_pos = MAX_READ_LEN + (br.get_start_pos() - self.pos.pos);
        let end_relative_pos = MAX_READ_LEN + (br.get_end_pos() - self.pos.pos);

        if start_relative_pos < 0 { return }
        if end_relative_pos >= MAX_READ_LEN * 2 + 1 { return }

        for relative_pos in start_relative_pos..end_relative_pos + 1 {
            new_read[relative_pos as usize] |= 1;
        }

        for cpg in br.get_cpgs().iter() {
            let relative_pos = self.get_relative_position(cpg.abspos);

            new_read[relative_pos] |= 2;

            if cpg.methylated {
                new_read[relative_pos] |= 4;
            }
        }

        // Reservoir sampling.
        // Fill if current reads are fewer than specified maximum depth.
        if self.num_total_read < self.max_depth as i32 {
            self.num_total_read += 1;
            self.num_sampled_read += 1;
            self.reads.push(new_read);
        }
        // Sample jth element and replace with current read with probability 1/num_total_read.
        else {
            self.num_total_read += 1;

            let j = rand::thread_rng().gen_range(1..self.num_total_read + 1);
            if j <= self.max_depth as i32 {
                self.reads[(j - 1) as usize] = new_read;
            }
        }
    }

    fn get_num_overlap_bases(&self, i: usize, j: usize) -> i32 {
        let r1 = self.reads[i];
        let r2 = self.reads[j];

        let mut num_overlap_bases = 0;
        for p in 0..MAX_READ_LEN * 2 + 1 {
            num_overlap_bases += ((r1[p as usize] & r2[p as usize]) & 1) as i32;
        }

        num_overlap_bases
    }

    fn get_num_overlap_cpgs(&self, i: usize, j: usize) -> i32 {
        let r1 = self.reads[i];
        let r2 = self.reads[j];

        let mut num_overlap_cpgs = 0;
        for p in 0..MAX_READ_LEN * 2 + 1 {
            num_overlap_cpgs += (((r1[p as usize] >> 1) & (r2[p as usize] >> 1)) & 1) as i32;
        }

        num_overlap_cpgs
    }

    fn hamming_distance(&self, i: usize, j: usize) -> f32 {
        let r1 = self.reads[i];
        let r2 = self.reads[j];

        let mut dist = 0.0;
        for p in 0..MAX_READ_LEN * 2 + 1 {
            if (r1[p as usize] & r2[p as usize]) & 3 == 3 {
                if ((r1[p as usize] ^ r2[p as usize]) & 4) >> 2 == 1 {
                    dist += 1.0;
                }
            }
        }

        dist
    }

    fn compute_qfdrp(&mut self, min_overlap: i32) -> f32 {
        let num_reads = self.get_num_reads();

        let mut qfdrp = 0.0;
        for comb in (0..num_reads).combinations(2) {
            let i = comb[0];
            let j = comb[1];

            // Read pair filtering.
            let num_overlap_bases = self.get_num_overlap_bases(i, j);
            let num_overlap_cpgs = self.get_num_overlap_cpgs(i, j);
            if num_overlap_bases < min_overlap {
                continue;
            }

            qfdrp += self.hamming_distance(i, j) / num_overlap_cpgs as f32;
        }

        qfdrp /= (num_reads * (num_reads - 1)) as f32 / 2.0;
        qfdrp
    }
}

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

    let reader = bamutil::get_reader(&input);
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
            .ok()
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
    let mut reader = bamutil::get_reader(&input);
    let header = bamutil::get_header(&reader);

    let mut readcount = 0;
    let mut valid_readcount = 0;

    let bar = progressbar::ProgressBar::new();

    let target_cpgs = &readutil::get_target_cpgs(cpg_set, &header);

    let mut cpg2reads: BTreeMap<readutil::CpGPosition, AssociatedReads> = BTreeMap::new();
    let mut result: BTreeMap<readutil::CpGPosition, f32> = BTreeMap::new();

    for r in reader.records().map(|r| r.unwrap()) {
        let mut br = readutil::BismarkRead::new(&r);
        match target_cpgs {
            Some(target_cpgs) => br.filter_isin(target_cpgs),
            None => {}
        }

        readcount += 1;
        if r.mapq() < min_qual {
            continue;
        }
        if br.get_num_cpgs() == 0 {
            continue;
        }

        match br.get_first_cpg_position() {
            Some(first_cpg_position) => {
                cpg2reads.retain(|&cpg, reads| {
                    let retain = {
                        if cpg < first_cpg_position {
                            if reads.get_num_reads() >= min_depth {
                                result.insert(cpg, reads.compute_qfdrp(min_overlap));
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

        for cpg_position in br.get_cpg_positions().iter() {
            let r = cpg2reads
                .entry(*cpg_position)
                .or_insert(AssociatedReads::new(*cpg_position, max_depth));

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
            result.insert(*cpg, reads.compute_qfdrp(min_overlap));
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::super::bamutil;
    use super::*;

    fn assert_approximately_equal(a: f32, b: f32) {
        assert_eq!((a - b) < 1e-5, true);
    }

    #[test]
    fn test_hamming_distance() {
        let input = "tests/test1.bam";
        let max_depth = 40;

        let mut reader = bamutil::get_reader(&input);

        let mut cpg2reads: BTreeMap<readutil::CpGPosition, AssociatedReads> = BTreeMap::new();

        for r in reader.records().map(|r| r.unwrap()) {
            let br = readutil::BismarkRead::new(&r);

            for cpg_position in br.get_cpg_positions().iter() {
                let r = cpg2reads
                    .entry(*cpg_position)
                    .or_insert(AssociatedReads::new(*cpg_position, max_depth));

                r.add_read(&br);
            }
        }

        for (_, reads) in cpg2reads.iter() {
            assert_eq!(reads.hamming_distance(0, 1), 1.0);
            assert_eq!(reads.hamming_distance(0, 2), 1.0);
            assert_eq!(reads.hamming_distance(0, 3), 2.0);
            assert_eq!(reads.hamming_distance(0, 4), 1.0);
            assert_eq!(reads.hamming_distance(0, 5), 2.0);
            assert_eq!(reads.hamming_distance(0, 6), 2.0);
            assert_eq!(reads.hamming_distance(0, 7), 3.0);
            assert_eq!(reads.hamming_distance(0, 8), 1.0);
            assert_eq!(reads.hamming_distance(0, 9), 2.0);
            assert_eq!(reads.hamming_distance(0, 10), 2.0);
            assert_eq!(reads.hamming_distance(0, 11), 3.0);
            assert_eq!(reads.hamming_distance(0, 12), 2.0);
            assert_eq!(reads.hamming_distance(0, 13), 3.0);
            assert_eq!(reads.hamming_distance(0, 14), 3.0);
            assert_eq!(reads.hamming_distance(0, 15), 4.0);
            break;
        }
    }
    #[test]
    fn test_get_num_reads() {
        let input = "tests/test1.bam";
        let max_depth = 40;

        let mut reader = bamutil::get_reader(&input);

        let mut cpg2reads: BTreeMap<readutil::CpGPosition, AssociatedReads> = BTreeMap::new();

        for r in reader.records().map(|r| r.unwrap()) {
            let br = readutil::BismarkRead::new(&r);

            for cpg_position in br.get_cpg_positions().iter() {
                let r = cpg2reads
                    .entry(*cpg_position)
                    .or_insert(AssociatedReads::new(*cpg_position, max_depth));

                r.add_read(&br);
            }
        }

        for (_, reads) in cpg2reads.iter() {
            assert_eq!(reads.get_num_reads(), 16);
        }
    }
    #[test]
    fn test_num_overlap_cpgs() {
        let input = "tests/test1.bam";
        let max_depth = 40;

        let mut reader = bamutil::get_reader(&input);

        let mut cpg2reads: BTreeMap<readutil::CpGPosition, AssociatedReads> = BTreeMap::new();

        for r in reader.records().map(|r| r.unwrap()) {
            let br = readutil::BismarkRead::new(&r);

            for cpg_position in br.get_cpg_positions().iter() {
                let r = cpg2reads
                    .entry(*cpg_position)
                    .or_insert(AssociatedReads::new(*cpg_position, max_depth));

                r.add_read(&br);
            }
        }

        for (_, reads) in cpg2reads.iter() {
            assert_eq!(reads.get_num_overlap_cpgs(0, 1), 4);
            break;
        }
    }
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
        for (i, (cpg, qfdrp)) in result.iter().enumerate() {
            assert_eq!(cpg.pos, cpg_positions[i]);
            assert_approximately_equal(*qfdrp, 8.0 / 15.0);
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
        for (i, (cpg, qfdrp)) in result.iter().enumerate() {
            assert_eq!(cpg.pos, cpg_positions[i]);
            assert_eq!(*qfdrp, 8.0 / 15.0);
            assert_approximately_equal(*qfdrp, 8.0 / 15.0);
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
        for (i, (cpg, qfdrp)) in result.iter().enumerate() {
            assert_eq!(cpg.pos, cpg_positions[i]);
            assert_eq!(*qfdrp, 1.0);
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
        for (i, (cpg, qfdrp)) in result.iter().enumerate() {
            assert_eq!(cpg.pos, cpg_positions[i]);
            assert_eq!(*qfdrp, 8.0 / 15.0);
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
