use rust_htslib::{bam, bam::Read};
use std::vec::Vec;
use std::fs;
use std::io::Write;
use std::str;
use std::cmp::Ordering;
use std::collections::{HashMap};

use crate::{readutil, bamutil, progressbar};

#[derive(Eq)]
struct AssociatedReads {
    pos: readutil::CpGPosition,
    stretch_info: HashMap<i32, i32>,
    num_cpgs: Vec<i32>,
}

impl AssociatedReads {
    fn new(pos: readutil::CpGPosition) -> Self {
        let stretch_info: HashMap<i32, i32> = HashMap::new();
        let num_cpgs: Vec<i32> = Vec::new();
        Self{ pos, stretch_info, num_cpgs }
    }

    fn get_coverage(&self) -> u32 { self.num_cpgs.len() as u32 }

    fn add_stretch_info(&mut self, stretch_info: HashMap<i32, i32>) {
        for (l, count) in stretch_info.iter() {
            let curr_count = self.stretch_info.entry(*l).or_insert(0);
            *curr_count += count;
        }
    }

    fn add_num_cpgs(&mut self, num_cpgs: usize) {
        self.num_cpgs.push(num_cpgs as i32);
    }

    fn to_bedgraph_field(&self, header: &bam::HeaderView) -> String {
        let chrom = bamutil::tid2chrom(self.pos.tid, header);
        let mut mhl = 0.0;
        let mut l_sum = 0.0;

        for ((l, count), num_cpgs) in self.stretch_info.iter().zip(self.num_cpgs.iter()) {
            mhl += (l * count) as f32 / (num_cpgs - l + 1) as f32;
            l_sum += *l as f32;
        }

        mhl /= l_sum;

        format!("{}\t{}\t{}\t{}", chrom, self.pos.pos, self.pos.pos + 2, mhl)
    }
}

impl PartialEq for AssociatedReads {
    fn eq(&self, other: &Self) -> bool { self.pos == other.pos }
}

impl PartialOrd for AssociatedReads {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> { Some(self.cmp(&other)) }
}

impl Ord for AssociatedReads {
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

    let mut pdr_result: HashMap<readutil::CpGPosition, AssociatedReads> = HashMap::new();

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
                            .or_insert(AssociatedReads::new(*cpg_position));
            
            r.add_num_cpgs(br.get_num_cpgs());
            r.add_stretch_info(br.get_stretch_info());
        }
        
        valid_readcount += 1;
        if readcount % 10000 == 0 { bar.update(readcount, valid_readcount) };
    }
    
    let mut results: Vec<&AssociatedReads> = pdr_result.values().collect::<Vec<&AssociatedReads>>();
    results.sort();

    let mut out = fs::OpenOptions::new().create(true).read(true).write(true).truncate(true).open(output).unwrap();

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

    let mut pdr_result: HashMap<readutil::CpGPosition, AssociatedReads> = HashMap::new();

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
                            .or_insert(AssociatedReads::new(*cpg_position));
            
            r.add_num_cpgs(br.get_num_cpgs());
            r.add_stretch_info(br.get_stretch_info());
        }
        
        valid_readcount += 1;
        if readcount % 10000 == 0 { bar.update(readcount, valid_readcount) };
    }
    
    let mut results: Vec<&AssociatedReads> = pdr_result.values().collect::<Vec<&AssociatedReads>>();
    results.sort();

    let mut out = fs::OpenOptions::new().create(true).read(true).write(true).truncate(true).open(output).unwrap();

    for r in results {
        if r.get_coverage() < min_depth { continue; } // Output filtering: CpG coverage should be >= min_depth.
        writeln!(out, "{}", r.to_bedgraph_field(&header))
            .ok()
            .expect("Error writing to output file.");
    }
}