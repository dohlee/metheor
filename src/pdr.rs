use rust_htslib::{bam, bam::Read, bam::record::{Record}};
use std::vec::Vec;
use std::fs;
use std::io::Write;
use std::str;
use std::fmt;
use std::cmp::Ordering;
use std::collections::{HashMap};

use crate::{readutil, bamutil, progressbar};

#[derive(Eq)]
struct PDRResult {
    pos: readutil::CpGPosition,
    n_concordant: u32,
    n_discordant: u32,
}

impl PDRResult {
    fn new(pos: readutil::CpGPosition) -> Self {
        Self{ pos: pos, n_concordant: 0, n_discordant: 0 }
    }

    fn inc_concordant(&mut self) { self.n_concordant += 1; }

    fn inc_discordant(&mut self) { self.n_discordant += 1; }

    fn get_coverage(&self) -> u32 { self.n_concordant + self.n_discordant }

    fn to_bedgraph_field(&self, header: &bam::HeaderView) -> String {
        let chrom = bamutil::tid2chrom(self.pos.tid, header);
        let pdr = (self.n_discordant as f32) / (self.n_concordant as f32 + self.n_discordant as f32);
        format!("{}\t{}\t{}\t{}\t{}\t{}", chrom, self.pos.pos, self.pos.pos + 2, pdr, self.n_concordant, self.n_discordant)
    }
}

impl fmt::Display for PDRResult {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let pdr = (self.n_discordant as f32) / (self.n_concordant as f32 + self.n_discordant as f32);
        write!(f, "{}\t{}\t{}\t{}", self.pos, pdr, self.n_concordant, self.n_discordant)
    }
}

impl PartialEq for PDRResult {
    fn eq(&self, other: &Self) -> bool { self.pos == other.pos }
}

impl PartialOrd for PDRResult {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> { Some(self.cmp(&other)) }
}

impl Ord for PDRResult {
    fn cmp(&self, other: &Self) -> Ordering { self.pos.cmp(&other.pos) }
}

pub fn compute(input: &str, output: &str, min_depth: u32, min_cpgs: usize, min_qual: u8, cpg_set: &Option<String>) {
    match cpg_set {
        Some(cpg_set) => compute_subset(input, output, min_depth, min_cpgs, min_qual, cpg_set),
        None => compute_all(input, output, min_depth, min_cpgs, min_qual),
    }
}

pub fn compute_subset(input: &str, output: &str, min_depth: u32, min_cpgs: usize, min_qual: u8, cpg_set: &str) {
    let mut reader = bamutil::get_reader(&input);
    let header = bamutil::get_header(&reader);

    let target_cpgs = readutil::get_target_cpgs(cpg_set, &header);

    let mut pdr_result: HashMap<readutil::CpGPosition, PDRResult> = HashMap::new();

    let mut readcount = 0;
    let mut valid_readcount = 0;

    let bar = progressbar::ProgressBar::new();

    for r in reader.records().map(|r| r.unwrap()) {
        let mut br = readutil::BismarkRead::new(&r);
        br.filter_isin(&target_cpgs);

        if br.get_num_cpgs() < min_cpgs { continue; }

        readcount += 1;
        if r.mapq() < min_qual { continue; } // Read filtering: Minimum quality should be >= min_qual.

        let mut cpg_positions = br.get_cpg_positions();
        if cpg_positions.len() == 0 { continue; } // Read filtering: Ignore reads without CpGs.

        for cpg_position in cpg_positions.iter_mut() {
            let r = pdr_result.entry(*cpg_position)
                            .or_insert(PDRResult::new(*cpg_position));
            
            match br.get_concordance_state() {
                readutil::ReadConcordanceState::Concordant => r.inc_concordant(),
                readutil::ReadConcordanceState::Discordant => r.inc_discordant(),
            }
        }
        
        valid_readcount += 1;
        if readcount % 10000 == 0 { bar.update_pdr(readcount, valid_readcount) };
    }
    
    let mut results: Vec<&PDRResult> = pdr_result.values().collect::<Vec<&PDRResult>>();
    results.sort();

    let mut out = fs::OpenOptions::new().create(true).read(true).write(true).open(output).unwrap();

    for r in results {
        if r.get_coverage() < min_depth { continue; } // Output filtering: CpG coverage should be >= min_depth.
        writeln!(out, "{}", r.to_bedgraph_field(&header))
            .ok()
            .expect("Error writing to output file.");
    }
}

pub fn compute_all(input: &str, output: &str, min_depth: u32, min_cpgs: usize, min_qual: u8) {
    let mut reader = bamutil::get_reader(&input);
    let header = bamutil::get_header(&reader);

    let mut pdr_result: HashMap<readutil::CpGPosition, PDRResult> = HashMap::new();

    let mut readcount = 0;
    let mut valid_readcount = 0;

    let bar = progressbar::ProgressBar::new();

    for r in reader.records().map(|r| r.unwrap()) {
        let br = readutil::BismarkRead::new(&r);

        readcount += 1;
        if br.get_num_cpgs() < min_cpgs { continue; } // Read filtering: Ignore reads with number of CpGs less than min_cpgs.
        if r.mapq() < min_qual { continue; } // Read filtering: Minimum quality should be >= min_qual.

        let mut cpg_positions = br.get_cpg_positions();

        for cpg_position in cpg_positions.iter_mut() {
            let r = pdr_result.entry(*cpg_position)
                            .or_insert(PDRResult::new(*cpg_position));
            
            match br.get_concordance_state() {
                readutil::ReadConcordanceState::Concordant => r.inc_concordant(),
                readutil::ReadConcordanceState::Discordant => r.inc_discordant(),
            }
        }
        
        valid_readcount += 1;
        if readcount % 10000 == 0 { bar.update_pdr(readcount, valid_readcount) };
    }
    
    let mut results: Vec<&PDRResult> = pdr_result.values().collect::<Vec<&PDRResult>>();
    results.sort();

    let mut out = fs::OpenOptions::new().create(true).read(true).write(true).open(output).unwrap();

    for r in results {
        if r.get_coverage() < min_depth { continue; } // Output filtering: CpG coverage should be >= min_depth.
        writeln!(out, "{}", r.to_bedgraph_field(&header))
            .ok()
            .expect("Error writing to output file.");
    }
}