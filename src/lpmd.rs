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

    fn print_pair_statistics(&self, output: &str) {
        let mut pairs: Vec<&(readutil::CpGPosition, readutil::CpGPosition)> = self.pair2n_concordant.keys().collect::<Vec<&(readutil::CpGPosition, readutil::CpGPosition)>>();
        pairs.sort();

        let mut out = fs::OpenOptions::new().create(true).read(true).write(true).open(output).unwrap();

        writeln!(out, "chrom\tcpg1\tcpg2\tlpmd\tn_concordant\tn_discordant")
            .ok()
            .expect("Error writing to output file.");

        for (cpg1, cpg2) in pairs {

            let k = (*cpg1, *cpg2);
            let n_concordant = self.pair2n_concordant[&k];
            let n_discordant = self.pair2n_discordant[&k];
            let lpmd = (n_discordant as f32) / (n_concordant as f32 + n_discordant as f32);

            let chrom = bamutil::tid2chrom(cpg1.tid, &self.header);

            writeln!(out, "{}\t{}\t{}\t{}\t{}\t{}", chrom, cpg1.pos, cpg2.pos, lpmd, n_concordant, n_discordant)
                .ok()
                .expect("Error writing to output file.");
        }
    }
}

pub fn compute(input: &str, output: &str, min_distance: i32, max_distance: i32, min_qual: u8, cpg_set: &Option<String>, pairs: &Option<String>) -> LPMDResult {

    match cpg_set {
        Some(cpg_set) => run_subset(input, output, min_distance, max_distance, min_qual, &cpg_set, pairs),
        None => run_all(input, output, min_distance, max_distance, min_qual, pairs),
    }

    // if cpg_set.is_empty() {
        // run_all(input, output, min_distance, max_distance, min_qual)
    // } else {
        // run_subset(input, output, min_distance, max_distance, cpg_set, min_qual)
    // }
}

fn run_all(input: &str, output: &str, min_distance: i32, max_distance: i32, min_qual: u8, pairs: &Option<String>) -> LPMDResult {
    eprintln!("Computing LPMD with parameters input={}, output={}, min_distance={}, max_distance={}", input, output, min_distance, max_distance);

    let res = compute_all(input, min_distance, max_distance, min_qual);

    let lpmd: f32 = res.compute_lpmd();
    let mut out = fs::OpenOptions::new().create(true).read(true).write(true).open(output).unwrap();

    writeln!(out, "name\tlpmd")
        .ok()
        .expect("Error writing to output file.");

    writeln!(out, "{}\t{}", input, lpmd)
        .ok()
        .expect("Error writing to output file.");

    match pairs {
        Some(f) => res.print_pair_statistics(f),
        None => (),
    }

    res
}

fn run_subset(input: &str, output: &str, min_distance: i32, max_distance: i32, min_qual: u8, cpg_set: &str, pairs: &Option<String>) -> LPMDResult {
    eprintln!("Computing subset-LPMD with parameters input={}, output={}, min_distance={}, max_distance={}", input, output, min_distance, max_distance);
    let mut reader = bamutil::get_reader(&input);
    let header = bamutil::get_header(&reader);

    eprint!("Processing target CpG set... ");
    let mut target_cpgs: HashSet<readutil::CpGPosition> = HashSet::new();

    let contents = fs::read_to_string(cpg_set)
                    .expect("Could not read target CpG file.");

    for line in contents.lines() {
        let tokens: Vec<&str> = line.split("\t").collect();

        let chrom = tokens[0];
        let pos = tokens[1].parse::<i32>().unwrap();
        
        target_cpgs.insert(readutil::CpGPosition{ tid: bamutil::chrom2tid(chrom.as_bytes(), &header) as i32, pos: pos });
    }

    eprintln!("Analyzing {} CpGs in total.", target_cpgs.len());

    let res = compute_subset(input, min_distance, max_distance, &target_cpgs, min_qual);
    let lpmd = res.compute_lpmd();
    let mut out = fs::OpenOptions::new().create(true).read(true).write(true).open(output).unwrap();

    writeln!(out, "name\tlpmd")
        .ok()
        .expect("Error writing to output file.");

    writeln!(out, "{}\t{}", input, lpmd)
        .ok()
        .expect("Error writing to output file.");

    match pairs {
        Some(f) => res.print_pair_statistics(f),
        None => (),
    }

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

fn compute_subset(input: &str, min_distance: i32, max_distance: i32, target_cpgs: &HashSet<readutil::CpGPosition>, min_qual: u8) -> LPMDResult {
    let mut reader = bamutil::get_reader(&input);
    let header = bamutil::get_header(&reader);
    let mut res = LPMDResult::new(header);
    // let _flag_counter: HashMap<u16, i32> = HashMap::new();

    let bar = progressbar::ProgressBar::new();

    // Iterate over reads and compute LPMD.
    for r in reader.records().map(|r| r.unwrap()) {
        res.inc_n_read(1);
        if r.mapq() < min_qual { continue; }

        let mut br = readutil::BismarkRead::new(&r);
        br.filter_isin(&target_cpgs);

        let (c, d, pair2concordance) = br.compute_pairwise_cpg_concordance_discordance(min_distance, max_distance);

        res.inc_n_valid_read(1);
        res.inc_n_concordant(c);
        res.inc_n_discordant(d);
        for (cpg1, cpg2, concordance) in &pair2concordance {
            res.add_pair_concordance(cpg1, cpg2, concordance);
        }
    
        if res.n_read % 10000 == 0 { bar.update_lpmd(res.progress_string()); }
    }

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test1() {
        // Maximum heterogeneity.
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
        // Half fully methylated, the other half fully unmethylated.
        let input = "tests/test2.bam";
        let min_distance = 2;
        let max_distance = 16;
        let min_qual = 10;
    
        let res = compute_all(input, min_distance, max_distance, min_qual);
    
        assert_eq!(res.n_concordant, 96);
        assert_eq!(res.n_discordant, 0);
    }

    #[test]
    fn test3() {
        // Same as test2, but there are only two reads with quality >= 10.
        let input = "tests/test3.bam";
        let min_distance = 2;
        let max_distance = 16;
        let min_qual = 10;

        let res = compute_all(input, min_distance, max_distance, min_qual);

        assert_eq!(res.n_concordant, 12);
        assert_eq!(res.n_discordant, 0);
    }

    #[test]
    fn test4_all() {
        // Similar to test1 (maximally heterogeneous methylation status),
        // but there are two independent bunches of reads.
        // Mimicking RRBS read clusters!
        let input = "tests/test4.bam";
        let min_distance = 2;
        let max_distance = 16;
        let min_qual = 10;

        let res = compute_all(input, min_distance, max_distance, min_qual);

        assert_eq!(res.n_concordant, 96);
        assert_eq!(res.n_discordant, 96);
    }

    #[test]
    fn test4_for_each_subset() {
        // Similar to test4_all, but compute only for each of the first and second bunch of reads.
        let input = "tests/test4.bam";
        let min_distance = 2;
        let max_distance = 16;
        let min_qual = 10;

        let mut target_cpgs: HashSet<readutil::CpGPosition> = HashSet::new();
        let pos = Vec::from([0, 2, 4, 6]);
        for p in pos.iter() {
            target_cpgs.insert(readutil::CpGPosition{ tid: 0, pos: *p });
        }

        let res = compute_subset(input, min_distance, max_distance, &target_cpgs, min_qual);

        assert_eq!(res.n_concordant, 48);
        assert_eq!(res.n_discordant, 48);

        let mut target_cpgs: HashSet<readutil::CpGPosition> = HashSet::new();
        let pos = Vec::from([13, 15, 17, 19]);
        for p in pos.iter() {
            target_cpgs.insert(readutil::CpGPosition{ tid: 0, pos: *p });
        }

        let res = compute_subset(input, min_distance, max_distance, &target_cpgs, min_qual);

        assert_eq!(res.n_concordant, 48);
        assert_eq!(res.n_discordant, 48);
    }

    fn test5() {
        // No valid reads.
        let input = "tests/test5.bam";
        let min_distance = 2;
        let max_distance = 16;
        let min_qual = 10;

        let res = compute_all(input, min_distance, max_distance, min_qual);

        assert_eq!(res.n_concordant, 0);
        assert_eq!(res.n_discordant, 0);
    }
}
