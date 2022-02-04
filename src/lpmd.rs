use rust_htslib::{bam, bam::Read, bam::ext::BamRecordExtensions, bam::record::{Record}};
use bio_types::genome::AbstractInterval;
use std::fs;
use std::io::Write;
use std::vec::Vec;
use std::str;
use std::collections::{HashSet, HashMap};
use indicatif::{HumanDuration};

use crate::{readutil, bamutil, progressbar};

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
        let pair2n_concordant: HashMap<(readutil::CpGPosition, readutil::CpGPosition), i32> = HashMap::new();
        let pair2n_discordant: HashMap<(readutil::CpGPosition, readutil::CpGPosition), i32> = HashMap::new();

        Self {
            header: header,
            n_read: 0,
            n_valid_read: 0,
            n_concordant: 0,
            n_discordant: 0,
            pair2n_concordant: pair2n_concordant,
            pair2n_discordant: pair2n_discordant,
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
        let lpmd: f32 = (self.n_discordant as f32) / ((self.n_concordant + self.n_discordant) as f32);
        lpmd
    }

    fn progress_string(&self) -> String {
        let lpmd = self.compute_lpmd();

        format!("Processed {} reads, found {} valid reads. LPMD={:.4} ({}/{})",
            self.n_read, self.n_valid_read, lpmd, self.n_discordant, self.n_concordant + self.n_discordant)
    }

    fn done_string(&self, t: HumanDuration) -> String {
        let lpmd = self.compute_lpmd();

        format!("Processed {} reads, found {} valid reads. LPMD={:.4} ({}/{}). Done in {}", 
            self.n_read, self.n_valid_read, lpmd, self.n_discordant, self.n_concordant + self.n_discordant, t)
    }

    fn add_pair_concordance(&mut self, pos1: &readutil::CpGPosition, pos2: &readutil::CpGPosition, concordance: &readutil::ReadConcordanceState) {
        
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

    fn print_pair_statistics(&self) {
        let mut pairs: Vec<&(readutil::CpGPosition, readutil::CpGPosition)> = self.pair2n_concordant.keys().collect::<Vec<&(readutil::CpGPosition, readutil::CpGPosition)>>();
        pairs.sort();

        for (cpg1, cpg2) in pairs {

            let k = (*cpg1, *cpg2);
            let n_concordant = self.pair2n_concordant[&k];
            let n_discordant = self.pair2n_discordant[&k];
            let lpmd = (n_discordant as f32) / (n_concordant as f32 + n_discordant as f32);

            let chrom = bamutil::tid2chrom(cpg1.tid, &self.header);
            println!("{}\t{}\t{}\t{}\t{}\t{}", chrom, cpg1.pos, cpg2.pos, lpmd, n_concordant, n_discordant);
        }
    }
}

pub fn compute(input: &str, output: &str, min_distance: i32, max_distance: i32, cpg_set: &str, min_qual: u8) -> LPMDResult {

    if cpg_set.is_empty() {
        run_all(input, output, min_distance, max_distance, min_qual)
    } else {
        run_subset(input, output, min_distance, max_distance, cpg_set, min_qual)
    }
}

fn run_all(input: &str, output: &str, min_distance: i32, max_distance: i32, min_qual: u8) -> LPMDResult {
    eprintln!("Computing LPMD with parameters input={}, output={}, min_distance={}, max_distance={}", input, output, min_distance, max_distance);

    let res = compute_all(input, min_distance, max_distance, min_qual);
    res.print_pair_statistics();

    let lpmd: f32 = res.compute_lpmd();
    let mut out = fs::OpenOptions::new().create(true).read(true).write(true).open(output).unwrap();

    write!(out, "name\tlpmd\n")
        .ok()
        .expect("Error writing to output file.");

    write!(out, "{}\t{}\n", input, lpmd)
        .ok()
        .expect("Error writing to output file.");

    res
}

fn run_subset(input: &str, output: &str, min_distance: i32, max_distance: i32, cpg_set: &str, min_qual: u8) -> LPMDResult {
    eprintln!("Computing subset-LPMD with parameters input={}, output={}, min_distance={}, max_distance={}", input, output, min_distance, max_distance);

    eprint!("Processing target CpG set... ");
    let mut target_cpgs: HashSet<(&str, i32)> = HashSet::new();

    let contents = fs::read_to_string(cpg_set)
                    .expect("Could not read target CpG file.");

    for line in contents.lines() {
        let tokens: Vec<&str> = line.split("\t").collect();

        let chrom = tokens[0];
        let pos = tokens[1].parse::<i32>().unwrap();
        
        target_cpgs.insert((chrom, pos));
    }

    eprintln!("Analyzing {} CpGs in total.", target_cpgs.len());

    let mut reader = bamutil::get_reader(&input);
    let header = bamutil::get_header(&reader);
    let mut res = LPMDResult::new(header);

    let _flag_counter: HashMap<u16, i32> = HashMap::new();

    let bar = progressbar::ProgressBar::new();

    // Iterate over reads and compute LPMD.
    for r in reader.records().map(|r| r.unwrap()) {
        res.inc_n_read(1);
        if r.mapq() < min_qual { continue; }

        let br = readutil::BismarkRead::new(&r);
        let (c, d, pair2concordance) = br.compute_pairwise_cpg_concordance_discordance(min_distance, max_distance);

        res.inc_n_valid_read(1);
        res.inc_n_concordant(c);
        res.inc_n_discordant(d);
        for (cpg1, cpg2, concordance) in &pair2concordance {
            res.add_pair_concordance(cpg1, cpg2, concordance);
        }
    
        if res.n_read % 10000 == 0 { bar.update_lpmd(res.progress_string()); }
    }

    let mut out = fs::OpenOptions::new().create(true).read(true).write(true).open(output).unwrap();
    let lpmd = res.compute_lpmd();

    write!(out, "name\tlpmd\n")
        .ok()
        .expect("Error writing to output file.");

    write!(out, "{}\t{}\n", input, lpmd)
        .ok()
        .expect("Error writing to output file.");

    res
}


fn compute_all(input: &str, min_distance: i32, max_distance: i32, min_qual: u8) -> LPMDResult {
    let mut reader = bamutil::get_reader(&input);
    let header = bamutil::get_header(&reader);
    let mut res = LPMDResult::new(header);

    let bar = progressbar::ProgressBar::new();

    for r in reader.records().map(|r| r.unwrap()) {
        res.inc_n_read(1);
        if r.mapq() < min_qual { continue; }

        let br = readutil::BismarkRead::new(&r);
        let (c, d, pair2concordance) = br.compute_pairwise_cpg_concordance_discordance(min_distance, max_distance);
        

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

    // bar.finish_with_message(res.done_string(HumanDuration(bar.elapsed())));

    res
}

fn get_all_cpgs_in_read(r: &Record, xm: &str) -> Vec<(i32, char)> {
    
    let mut cpgs: Vec<(i32, char)> = Vec::new();
    for (pos, c) in r.reference_positions_full().zip(xm.chars()) {
        
        if (c != 'z') && (c != 'Z') { continue; }

        match pos {
            Some(pos) => {
                if (r.flags() == 99) || (r.flags() == 147) { // Forward
                    cpgs.push((pos as i32, c));
                } else {
                    cpgs.push(((pos - 1) as i32, c));
                }
            },
            None => {}
        }
    }

    cpgs
}

fn get_subset_cpgs_in_read(r: &Record, xm: &str, chrom: &str, target_cpgs: &HashSet<(&str, i32)>) -> Vec<(i32, char)> {
    
    let mut cpgs: Vec<(i32, char)> = Vec::new();
    for (pos, c) in r.reference_positions_full().zip(xm.chars()) {
        
        if (c != 'z') && (c != 'Z') { continue; }

        match pos {
            Some(pos) => {
                if (r.flags() == 99) || (r.flags() == 147) { // Forward
                    let this_cpg = (chrom, pos as i32);
                    if target_cpgs.contains(&this_cpg) {
                        cpgs.push((pos as i32, c));
                    }
                } else {
                    let this_cpg = (chrom, (pos - 1) as i32);
                    if target_cpgs.contains(&this_cpg) {
                        cpgs.push(((pos - 1) as i32, c));
                    }
                }
            },
            None => {}
        }
    }

    cpgs
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
    
        let res = compute_all(input, min_distance, max_distance, min_qual);
    
        assert_eq!(res.n_concordant, 48);
        assert_eq!(res.n_discordant, 48);
    }

    #[test]
    fn test2() {

        let input = "tests/test2.bam";
        let min_distance = 2;
        let max_distance = 16;
        let min_qual = 10;
    
        let res = compute_all(input, min_distance, max_distance, min_qual);
    
        assert_eq!(res.n_concordant, 96);
        assert_eq!(res.n_discordant, 0);
    }
}
