use rand::Rng;
use std::cmp::Ordering;

use crate::readutil::BismarkRead;

pub(crate) struct AssociatedReads {
    reads: Vec<BismarkRead>,
    num_total_reads: usize,
    max_depth: usize,
}

impl AssociatedReads {
    pub fn new(max_depth: usize) -> Self {
        Self {
            reads: Vec::new(),
            num_total_reads: 0,
            max_depth,
        }
    }

    pub fn get_num_reads(&self) -> usize {
        self.reads.len()
    }

    pub fn add_read(&mut self, read: &BismarkRead) {
        self.num_total_reads += 1;
        if self.reads.len() < self.max_depth {
            self.reads.push(read.clone());
        } else {
            let index = rand::thread_rng().gen_range(0..self.num_total_reads);
            if index < self.max_depth {
                self.reads[index] = read.clone();
            }
        }
    }

    fn pair_counts(&self, i: usize, j: usize) -> (usize, usize) {
        let left = self.reads[i].get_cpgs();
        let right = self.reads[j].get_cpgs();
        let (mut a, mut b, mut common, mut discordant) = (0, 0, 0, 0);
        while a < left.len() && b < right.len() {
            match left[a].abspos.cmp(&right[b].abspos) {
                Ordering::Less => a += 1,
                Ordering::Greater => b += 1,
                Ordering::Equal => {
                    common += 1;
                    discordant += usize::from(left[a].methylated != right[b].methylated);
                    a += 1;
                    b += 1;
                }
            }
        }
        (common, discordant)
    }

    #[cfg(test)]
    pub fn get_num_overlap_cpgs(&self, i: usize, j: usize) -> usize {
        self.pair_counts(i, j).0
    }

    #[cfg(test)]
    pub fn hamming_distance(&self, i: usize, j: usize) -> f32 {
        self.pair_counts(i, j).1 as f32
    }

    fn compute(&self, min_overlap: i32, quantitative: bool) -> f32 {
        let mut score = 0.0;
        for i in 0..self.reads.len() {
            for j in i + 1..self.reads.len() {
                let left = &self.reads[i];
                let right = &self.reads[j];
                let start = i64::from(left.get_start_pos().max(right.get_start_pos()));
                let end = i64::from(left.get_end_pos().min(right.get_end_pos()));
                if end - start + 1 < i64::from(min_overlap) {
                    continue;
                }
                let (common, discordant) = self.pair_counts(i, j);
                if common > 0 {
                    score += if quantitative {
                        discordant as f32 / common as f32
                    } else {
                        f32::from(discordant > 0)
                    };
                }
            }
        }
        let n = self.reads.len();
        if n < 2 {
            return f32::NAN;
        }
        score / (n * (n - 1)) as f32 * 2.0
    }

    pub fn compute_fdrp(&self, min_overlap: i32) -> f32 {
        self.compute(min_overlap, false)
    }
    pub fn compute_qfdrp(&self, min_overlap: i32) -> f32 {
        self.compute(min_overlap, true)
    }
}
