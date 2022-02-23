use rust_htslib::{bam, bam::Read, bam::ext::BamRecordExtensions, bam::record::{Record}};
use std::fs;
use std::fs::OpenOptions;
use std::vec::Vec;
use std::str;
use std::io::Write;
use std::collections::{HashMap};

use crate::{readutil, bamutil, progressbar};

struct MEResult {
    pos1: readutil::CpGPosition,
    pos2: readutil::CpGPosition,
    pos3: readutil::CpGPosition,
    pos4: readutil::CpGPosition,
    quartet_pattern_counts: [u32; 16],
}

impl MEResult {

    fn new(q: readutil::Quartet) -> Self {
        let pos1 = q.pos1;
        let pos2 = q.pos2;
        let pos3 = q.pos3;
        let pos4 = q.pos4;

        let quartet_pattern_counts = [0; 16];
        Self{ pos1, pos2, pos3, pos4, quartet_pattern_counts }
    }

    fn get_read_depth(&self) -> u32 {
        self.quartet_pattern_counts.iter().sum()
    }

    fn add_quartet_pattern(&mut self, p: readutil::QuartetPattern) {
        self.quartet_pattern_counts[p] += 1;
    }

    fn to_bedgraph_field(&self, header: &bam::HeaderView) -> String {
        let chrom = bamutil::tid2chrom(self.pos1.tid, header);
        let mut me: f32 = 0.0;

        let total: u32 = self.quartet_pattern_counts.iter().sum();
        for count in self.quartet_pattern_counts.iter() {
            let p: f32 = (*count as f32) / (total as f32);
            if *count > 0 { me += p * p.log2(); }
        }
        me *= -0.25;

        format!("{}\t{}\t{}\t{}\t{}\t{}", chrom, self.pos1.pos, self.pos2.pos, self.pos3.pos, self.pos4.pos, me)
    }
}

pub fn compute(input: &str, output: &str, min_depth: u32, min_qual: u8, cpg_set: &Option<String>) {
    match cpg_set {
        Some(cpg_set) => compute_subset(input, output, min_depth, min_qual, cpg_set),
        None => compute_all(input, output, min_depth, min_qual),
    }
}

pub fn compute_all(input: &str, output: &str, min_depth: u32, min_qual: u8) {
    let mut reader = bamutil::get_reader(&input);
    let header = bamutil::get_header(&reader);

    let mut quartet2stat: HashMap<readutil::Quartet, MEResult> = HashMap::new();

    let mut readcount = 0;
    let mut valid_readcount = 0;

    let bar = progressbar::ProgressBar::new();

    for r in reader.records().map(|r| r.unwrap()) {
        let br = readutil::BismarkRead::new(&r);

        readcount += 1;
        
        if r.mapq() < min_qual { continue; }
        valid_readcount += 1;

        let (quartets, patterns) = br.get_cpg_quartets_and_patterns();
        for (q, p) in quartets.iter().zip(patterns.iter()) {
            let stat = quartet2stat.entry(*q)
                        .or_insert(MEResult::new(*q));

            stat.add_quartet_pattern(*p);
        }
    }

    let mut out = fs::OpenOptions::new().create(true).read(true).write(true).truncate(true).open(output).unwrap();
    for stat in quartet2stat.values() {
        if stat.get_read_depth() < min_depth { continue; }
        writeln!(out, "{}", stat.to_bedgraph_field(&header))
            .ok()
            .expect("Error writing to output file.");
    }
}

pub fn compute_subset(input: &str, output: &str, min_depth: u32, min_qual: u8, cpg_set: &str) {
    let mut reader = bamutil::get_reader(&input);
    let header = bamutil::get_header(&reader);
    
    let target_cpgs = readutil::get_target_cpgs(cpg_set, &header);
    let mut quartet2stat: HashMap<readutil::Quartet, MEResult> = HashMap::new();

    let mut readcount = 0;
    let mut valid_readcount = 0;

    let bar = progressbar::ProgressBar::new();

    for r in reader.records().map(|r| r.unwrap()) {
        let mut br = readutil::BismarkRead::new(&r);
        br.filter_isin(&target_cpgs);

        readcount += 1;
        
        if r.mapq() < min_qual { continue; }
        valid_readcount += 1;

        let (quartets, patterns) = br.get_cpg_quartets_and_patterns();
        for (q, p) in quartets.iter().zip(patterns.iter()) {
            let stat = quartet2stat.entry(*q)
                        .or_insert(MEResult::new(*q));

            stat.add_quartet_pattern(*p);
        }
    }

    let mut out = fs::OpenOptions::new().create(true).read(true).write(true).truncate(true).open(output).unwrap();
    for stat in quartet2stat.values() {
        if stat.get_read_depth() < min_depth { continue; }
        writeln!(out, "{}", stat.to_bedgraph_field(&header))
            .ok()
            .expect("Error writing to output file.");
    }
}
