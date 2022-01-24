use rust_htslib::{bam::Read, bam::ext::BamRecordExtensions, bam::record::{Aux, Record}};

use std::collections::{HashSet};
use std::fmt;

pub struct BismarkRead {
    // Read is defined as an array of CpG methylation states
    // and their relative/absolute positions.
    cpgs: Vec<CpG>
}

impl BismarkRead {
    pub fn new(r: &Record) -> Self {
        
        match r.aux(b"XM") {
            Ok(value) => {
                if let Aux::String(xm) = value {  // if value is a type of Aux::String, run:
                    let cpgs = get_cpgs(r, xm);

                    Self { cpgs: cpgs }
                } else {
                    panic!("Error reading XM tag in BAM record. Make sure the reads are aligned using Bismark!");
                }
            }
            Err(_) => {
                panic!("Error reading XM tag in BAM record. Make sure the reads are aligned using Bismark!");
            }
        }
    }

    pub fn get_num_cpgs(&self) -> u32 {
        self.cpgs.len() as u32
    }

    pub fn get_cpg_positions(&self) -> Vec<CpGPosition> {
        let mut s: Vec<CpGPosition> = Vec::new();
        for cpg in &self.cpgs {
            s.push(cpg.abspos);
        }

        s
    }

    pub fn get_first_cpg_position(&self) -> &CpGPosition {
        &self.cpgs[0].abspos
    }

    pub fn get_concordance_state(&self) -> ReadConcordanceState {
        let init_methylated = self.cpgs[0].methylated;
        let mut res = ReadConcordanceState::Concordant;

        for cpg in &self.cpgs {
            if cpg.methylated != init_methylated { res = ReadConcordanceState::Discordant; }
        }

        res
    }
}

pub enum ReadConcordanceState {
    Concordant,
    Discordant, 
}

#[derive(Eq, Hash, Copy)]
pub struct CpGPosition {
    tid: i32,
    pos: i32,
}

impl CpGPosition {
    pub fn new(tid: i32, pos: i32) -> Self {
        Self { tid: tid, pos: pos }
    }
}

impl fmt::Display for CpGPosition {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\t{}\t{}", self.tid, self.pos, self.pos+2)
    }
}

impl PartialEq for CpGPosition {
    fn eq(&self, other: &Self) -> bool {
        (self.tid == other.tid) && (self.pos == other.pos)
    }
}

impl Clone for CpGPosition {
    fn clone(&self) -> Self {
        *self
    }
}

struct CpG {
    relpos: i32,
    abspos: CpGPosition,
    methylated: bool
}

impl CpG {
    fn new(relpos: i32, abspos: CpGPosition, c: char) -> Self {
        Self { relpos:relpos, abspos: abspos, methylated: c == 'Z' }
    }
}

impl ToString for CpG {
    fn to_string(&self) -> String {
        format!("{}, {}, {}, {}", self.relpos, self.abspos.tid, self.abspos.pos, self.methylated)
    }
}

fn get_cpgs(r: &Record, xm: &str) -> Vec<CpG> {
    let mut cpgs:Vec<CpG> = Vec::new();

    for (relpos, (abspos, c)) in r.reference_positions_full().zip(xm.chars()).enumerate() {
    
        if (c != 'z') && (c != 'Z') { continue; }

        match abspos {
            Some(abspos) => {
                if (r.flags() == 99) || (r.flags() == 147) { // Forward
                    let cpgpos = CpGPosition::new(r.tid() as i32, abspos as i32);
                    cpgs.push(CpG::new(relpos as i32, cpgpos, c));
                } else {
                    let cpgpos = CpGPosition::new(r.tid() as i32, (abspos - 1) as i32);
                    cpgs.push(CpG::new(relpos as i32, cpgpos, c));
                }
            },
            None => {}
        }
    }

    return cpgs
}

pub fn count_z(meth_str: &str) -> i32 {
    (meth_str.matches("z").count() + meth_str.matches("Z").count()) as i32
}

pub fn compute_pairwise_concordance_discordance_from_read(cpgs: Vec<(i32, char)>, min_distance: i32, max_distance: i32) -> (i32, i32) {
    let mut anchors: Vec<CpG> = Vec::new();
    let mut min_anchor_pos = -1;
    let mut n_concordant = 0;
    let mut n_discordant = 0;

    // for (i, c) in meth_str.chars().enumerate().filter(|(_i, c)| (*c == 'z') || (*c == 'Z')) {
    for (i, c) in cpgs.iter(){
        let dummy_cpgpos = CpGPosition::new(-1, -1);
        let cpg_state = CpG::new(*i as i32, dummy_cpgpos, *c);

        if min_anchor_pos != -1 {
            while (cpg_state.relpos - min_anchor_pos > max_distance) && (anchors.len() > 0) {
                anchors.remove(0);

                if anchors.len() > 0 {
                    min_anchor_pos = anchors[0].relpos;
                } else {
                    min_anchor_pos = -1;
                }
            }
        }

        for anchor in anchors.iter() {
            if cpg_state.relpos - anchor.relpos < min_distance {
                continue;
            }
            
            if anchor.methylated == cpg_state.methylated {
                n_concordant += 1;
            } else {
                n_discordant += 1;
            }
        }

        if min_anchor_pos == -1 { min_anchor_pos = cpg_state.relpos; }
        anchors.push(cpg_state);
    }

    // println!("{}, {}", n_concordant, n_discordant);
    (n_concordant, n_discordant)
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::bamutil;

    #[test]
    fn pairwise_concordance_discordance_no_pair() {
        let min_distance = 2;
        let max_distance = 16;
        let tmp_cpgs: Vec<(i32, char)> = Vec::from([(2, 'z'), (100, 'z')]);
        let (n_con, n_dis) = compute_pairwise_concordance_discordance_from_read(tmp_cpgs, min_distance, max_distance);

        assert_eq!(n_con, 0);
        assert_eq!(n_dis, 0);
    }

    #[test]
    fn pairwise_concordance_discordance_complex() {
        let min_distance = 2;
        let max_distance = 16;
        let tmp_cpgs: Vec<(i32, char)> = Vec::from([(2, 'z'), (8, 'Z'), (20, 'z'), (32, 'z'), (34, 'z'), (100, 'Z'), (116, 'z')]);
        let (n_con, n_dis) = compute_pairwise_concordance_discordance_from_read(tmp_cpgs, min_distance, max_distance);

        assert_eq!(n_con, 3);
        assert_eq!(n_dis, 3);
    }

    #[test]
    fn test_bismarkread_constructor() {
        let input = "tests/test1.bam";
        let mut reader = bamutil::get_reader(&input);
        for r in reader.records() {
            let mut r = r.unwrap();

            let br = BismarkRead::new(&r);
        }
    }

    #[test]
    fn test_bismarkread_get_concordance_state() {
        let input = "tests/test1.bam";
        let mut reader = bamutil::get_reader(&input);
        for r in reader.records() {
            let mut r = r.unwrap();

            let br = BismarkRead::new(&r);

            br.get_concordance_state();
        }
    }

    #[test]
    fn test_test1_pdr() {
        let input = "tests/test1.bam";
        let mut reader = bamutil::get_reader(&input);
        let mut n_read = 0;
        let mut n_discordant_read = 0;
        for r in reader.records() {
            let mut r = r.unwrap();

            let br = BismarkRead::new(&r);

            n_read += 1;
            match br.get_concordance_state() {
                ReadConcordanceState::Concordant => {},
                ReadConcordanceState::Discordant => { n_discordant_read += 1 }
            }
        }

        assert_eq!(n_read, 16);
        assert_eq!(n_discordant_read, 14);
    }
}
