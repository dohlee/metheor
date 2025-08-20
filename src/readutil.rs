use rust_htslib::{
    bam,
    bam::ext::BamRecordExtensions,
    bam::record::{Aux, Record},
};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::fs;

use crate::bamutil;

pub type QuartetPattern = usize;

pub struct BismarkRead {
    start_pos: i32,
    end_pos: i32,
    // Read is defined as an array of CpG methylation states
    // and their relative/absolute positions.
    cpgs: Vec<CpG>,
}

impl BismarkRead {
    pub fn new(r: &Record) -> Self {
        let mut start_pos = -1;
        let mut end_pos = -1;

        for abspos in r.reference_positions_full().flatten() {
            if start_pos == -1 {
                start_pos = abspos as i32;
            }
            end_pos = abspos as i32;
        }

        match r.aux(b"XM") {
            Ok(value) => {
                if let Aux::String(xm) = value {
                    // if value is a type of Aux::String, run:
                    let cpgs = get_cpgs(r, xm);
                    Self {
                        start_pos,
                        end_pos,
                        cpgs,
                    }
                } else {
                    panic!("Error reading XM tag in BAM record. Make sure the reads are aligned using Bismark!");
                }
            }
            Err(_) => {
                panic!("Error reading XM tag in BAM record. Make sure the reads are aligned using Bismark!");
            }
        }
    }

    pub fn get_first_cpg_position(&self) -> Option<CpGPosition> {
        match self.get_num_cpgs() {
            0 => None,
            _ => Some(self.cpgs[0].abspos),
        }
    }

    pub fn get_start_pos(&self) -> i32 {
        self.start_pos
    }

    pub fn get_end_pos(&self) -> i32 {
        self.end_pos
    }

    pub fn get_num_cpgs(&self) -> usize {
        self.cpgs.len()
    }

    pub fn get_cpgs(&self) -> &Vec<CpG> {
        &self.cpgs
    }

    pub fn get_cpg_positions(&self) -> Vec<CpGPosition> {
        let mut s: Vec<CpGPosition> = Vec::new();
        for cpg in &self.cpgs {
            s.push(cpg.abspos);
        }

        s
    }

    pub fn filter_isin(&mut self, target_cpgs: &HashSet<CpGPosition>) {
        let mut new_cpgs: Vec<CpG> = Vec::new();

        for cpg in self.cpgs.iter().filter(|x| target_cpgs.contains(&x.abspos)) {
            new_cpgs.push(*cpg);
        }

        self.cpgs = new_cpgs;
    }

    pub fn get_cpg_quartets_and_patterns(&self) -> (Vec<Quartet>, Vec<QuartetPattern>) {
        let mut quartets: Vec<Quartet> = Vec::new();
        let mut patterns: Vec<QuartetPattern> = Vec::new();

        if self.get_num_cpgs() < 4 {
            return (quartets, patterns);
        }

        for i in 0..self.get_num_cpgs() - 3 {
            let q = Quartet {
                pos1: self.cpgs[i].abspos,
                pos2: self.cpgs[i + 1].abspos,
                pos3: self.cpgs[i + 2].abspos,
                pos4: self.cpgs[i + 3].abspos,
            };
            let mut p = 0;

            if self.cpgs[i].methylated {
                p += 8;
            }
            if self.cpgs[i + 1].methylated {
                p += 4;
            }
            if self.cpgs[i + 2].methylated {
                p += 2;
            }
            if self.cpgs[i + 3].methylated {
                p += 1;
            }

            quartets.push(q);
            patterns.push(p);
        }

        (quartets, patterns)
    }

    pub fn get_concordance_state(&self) -> ReadConcordanceState {
        let init_methylated = self.cpgs[0].methylated;
        let mut res = ReadConcordanceState::Concordant;

        for cpg in &self.cpgs {
            if cpg.methylated != init_methylated {
                res = ReadConcordanceState::Discordant;
            }
        }

        res
    }

    pub fn get_stretch_info(&self) -> HashMap<i32, i32> {
        let mut stretch_info: HashMap<i32, i32> = HashMap::new();
        let mut curr_stretch_length = 0;

        for cpg in &self.cpgs {
            if cpg.methylated {
                curr_stretch_length += 1;
                for l in 1..curr_stretch_length + 1 {
                    let v = stretch_info.entry(l).or_insert(0);
                    *v += 1;
                }
            } else {
                curr_stretch_length = 0;
            }
        }

        stretch_info
    }

    pub fn compute_pairwise_cpg_concordance_discordance(
        &self,
        min_distance: i32,
        max_distance: i32,
    ) -> (
        i32,
        i32,
        Vec<(CpGPosition, CpGPosition, ReadConcordanceState)>,
    ) {
        let mut anchors: Vec<CpG> = Vec::new();
        let mut pair2concordance: Vec<(CpGPosition, CpGPosition, ReadConcordanceState)> =
            Vec::new();
        let mut min_anchor_pos = -1;
        let mut n_concordant = 0;
        let mut n_discordant = 0;

        for cpg in &self.cpgs {
            if min_anchor_pos != -1 {
                while (cpg.relpos - min_anchor_pos > max_distance) && !anchors.is_empty() {
                    anchors.remove(0);

                    if !anchors.is_empty() {
                        min_anchor_pos = anchors[0].relpos;
                    } else {
                        min_anchor_pos = -1;
                    }
                }
            }

            for anchor in anchors.iter() {
                if cpg.relpos - anchor.relpos < min_distance {
                    continue;
                }

                if anchor.methylated == cpg.methylated {
                    n_concordant += 1;
                    pair2concordance.push((
                        anchor.abspos,
                        cpg.abspos,
                        ReadConcordanceState::Concordant,
                    ));
                } else {
                    n_discordant += 1;
                    pair2concordance.push((
                        anchor.abspos,
                        cpg.abspos,
                        ReadConcordanceState::Discordant,
                    ));
                }
            }

            if min_anchor_pos == -1 {
                min_anchor_pos = cpg.relpos;
            }
            anchors.push(*cpg);
        }

        (n_concordant, n_discordant, pair2concordance)
    }
}

#[derive(Eq, PartialEq, Hash, Copy)]
pub struct Quartet {
    pub pos1: CpGPosition,
    pub pos2: CpGPosition,
    pub pos3: CpGPosition,
    pub pos4: CpGPosition,
}

impl Clone for Quartet {
    fn clone(&self) -> Self {
        *self
    }
}

pub enum ReadConcordanceState {
    Concordant,
    Discordant,
}

#[derive(Copy)]
pub struct CpG {
    pub relpos: i32,
    pub abspos: CpGPosition,
    pub methylated: bool,
}

impl CpG {
    fn new(relpos: i32, abspos: CpGPosition, c: char) -> Self {
        Self {
            relpos,
            abspos,
            methylated: c == 'Z',
        }
    }
}

impl Clone for CpG {
    fn clone(&self) -> Self {
        *self
    }
}

impl fmt::Display for CpG {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}, {}, {}, {}",
            self.relpos, self.abspos.tid, self.abspos.pos, self.methylated
        )
    }
}

#[derive(Eq, PartialEq, Hash, Copy)]
pub struct CpGPosition {
    pub tid: i32,
    pub pos: i32,
}

impl CpGPosition {
    pub fn new(tid: i32, pos: i32) -> Self {
        Self { tid, pos }
    }

    pub fn is_before(&self, other: &Self, distance: i32) -> bool {
        match self.tid.cmp(&other.tid) {
            Ordering::Greater => false,
            Ordering::Less => true,
            Ordering::Equal => self.pos + distance < other.pos,
        }
    }
}

impl fmt::Display for CpGPosition {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\t{}\t{}", self.tid, self.pos, self.pos + 2)
    }
}

impl PartialOrd for CpGPosition {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for CpGPosition {
    fn cmp(&self, other: &Self) -> Ordering {
        self.tid.cmp(&other.tid).then(self.pos.cmp(&other.pos))
    }
}

impl Clone for CpGPosition {
    fn clone(&self) -> Self {
        *self
    }
}

fn get_cpgs(r: &Record, xm: &str) -> Vec<CpG> {
    let mut cpgs: Vec<CpG> = Vec::new();

    for (relpos, (abspos, c)) in r.reference_positions_full().zip(xm.chars()).enumerate() {
        if (c != 'z') && (c != 'Z') {
            continue;
        }

        if let Some(abspos) = abspos {
            if (r.flags() == 0) || (r.flags() == 99) || (r.flags() == 147) {
                // Forward
                let cpgpos = CpGPosition::new(r.tid(), abspos as i32);
                cpgs.push(CpG::new(relpos as i32, cpgpos, c));
            } else {
                // Reverse
                let cpgpos = CpGPosition::new(r.tid(), (abspos - 1) as i32);
                cpgs.push(CpG::new(relpos as i32, cpgpos, c));
            }
        }
    }

    cpgs
}

pub fn get_target_cpgs(
    cpg_set: &Option<String>,
    header: &bam::HeaderView,
) -> Option<HashSet<CpGPosition>> {
    match cpg_set {
        Some(cpg_set) => {
            eprint!("Processing target CpG set... ");
            let mut target_cpgs: HashSet<CpGPosition> = HashSet::new();

            let contents = fs::read_to_string(cpg_set).expect("Could not read target CpG file.");

            for line in contents.lines() {
                let tokens: Vec<&str> = line.split("\t").collect();

                let chrom = tokens[0];
                let pos = tokens[1].parse::<i32>().unwrap();

                target_cpgs.insert(CpGPosition {
                    tid: bamutil::chrom2tid(chrom.as_bytes(), header) as i32,
                    pos,
                });
            }

            Some(target_cpgs)
        }
        None => None,
    }
}

#[cfg(test)]
mod tests {
    use super::super::bamutil;
    use super::*;
    use rust_htslib::bam::Read;

    #[test]
    fn test_bismarkread_constructor() {
        let input = "tests/test1.bam";
        let mut reader = bamutil::get_reader(input);
        for r in reader.records() {
            let r = r.unwrap();
            let _br = BismarkRead::new(&r);
        }
    }

    #[test]
    fn test_bismarkread_get_concordance_state() {
        let input = "tests/test1.bam";
        let mut reader = bamutil::get_reader(input);
        for r in reader.records() {
            let r = r.unwrap();
            let br = BismarkRead::new(&r);

            br.get_concordance_state();
        }
    }

    #[test]
    fn pairwise_concordance_discordance_no_pair() {
        let min_distance = 2;
        let max_distance = 16;

        let input = "tests/test1.bam";
        let mut reader = bamutil::get_reader(input);
        for r in reader.records() {
            let r = r.unwrap();

            let br = BismarkRead::new(&r);

            br.compute_pairwise_cpg_concordance_discordance(min_distance, max_distance);
        }
    }

    #[test]
    fn test_test1_pdr() {
        let input = "tests/test1.bam";
        let mut reader = bamutil::get_reader(input);
        let mut n_read = 0;
        let mut n_discordant_read = 0;
        for r in reader.records() {
            let r = r.unwrap();

            let br = BismarkRead::new(&r);

            n_read += 1;
            match br.get_concordance_state() {
                ReadConcordanceState::Concordant => {}
                ReadConcordanceState::Discordant => n_discordant_read += 1,
            }
        }

        assert_eq!(n_read, 16);
        assert_eq!(n_discordant_read, 14);
    }

    #[test]
    fn test_cpgposition_eq() {
        let pos1 = CpGPosition { tid: 0, pos: 1 };
        let pos2 = CpGPosition { tid: 0, pos: 1 };
        let pos3 = CpGPosition { tid: 0, pos: 2 };
        let pos4 = CpGPosition { tid: 1, pos: 1 };

        assert!(pos1 == pos2);
        assert!(pos1 != pos3);
        assert!(pos1 != pos4);
    }

    #[test]
    fn test_cpgposition_ordering() {
        let pos1 = CpGPosition { tid: 0, pos: 1 };
        let pos2 = CpGPosition { tid: 0, pos: 1 };
        let pos3 = CpGPosition { tid: 0, pos: 2 };
        let pos4 = CpGPosition { tid: 1, pos: 1 };

        assert!(pos1 <= pos2);
        assert!(pos1 >= pos2);

        assert!(pos1 < pos3);
        assert!(pos1 <= pos3);

        assert!(pos1 < pos4);
        assert!(pos3 < pos4);
    }
}
