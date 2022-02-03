use rust_htslib::{bam, bam::Read, bam::ext::BamRecordExtensions, bam::record::{Record}};
use std::fs::OpenOptions;
use std::vec::Vec;
use std::str;
use std::io::Write;
use std::collections::{HashMap};

use crate::{readutil, bamutil, progressbar};

struct PMResult {
    pos1: readutil::CpGPosition,
    pos2: readutil::CpGPosition,
    pos3: readutil::CpGPosition,
    pos4: readutil::CpGPosition,
    quartet_pattern_counts: [u32; 16],
}

impl PMResult {

    fn new(q: readutil::Quartet) -> Self {
        let pos1 = q.pos1;
        let pos2 = q.pos2;
        let pos3 = q.pos3;
        let pos4 = q.pos4;

        let quartet_pattern_counts = [0; 16];
        Self{ pos1, pos2, pos3, pos4, quartet_pattern_counts }
    }

    fn add_quartet_pattern(&mut self, p: readutil::QuartetPattern) {
        self.quartet_pattern_counts[p] += 1;
    }

    fn to_bedgraph_field(&self, header: &bam::HeaderView) -> String {
        let chrom = bamutil::tid2chrom(self.pos1.tid, header);
        let mut pm = 1.0;
        let total: u32 = self.quartet_pattern_counts.iter().sum();
        for count in self.quartet_pattern_counts.iter() {
            pm -= ((*count as f32) / (total as f32)) * ((*count as f32) / (total as f32));
        }
        format!("{}\t{}\t{}\t{}\t{}\t{}", chrom, self.pos1.pos, self.pos2.pos, self.pos3.pos, self.pos4.pos, pm)
    }
}

pub fn compute(input: &str, _output: &str, min_depth: u32, min_qual: u8) {
    let mut reader = bamutil::get_reader(&input);
    let header = bamutil::get_header(&reader);

    let mut quartet2stat: HashMap<readutil::Quartet, PMResult> = HashMap::new();

    let mut readcount = 0;
    let mut valid_readcount = 0;

    let bar = progressbar::ProgressBar::new();

    for r in reader.records().map(|r| r.unwrap()) {
        let br = readutil::BismarkRead::new(&r);

        readcount += 1;

        // TODO: read filtering.
        valid_readcount += 1;

        let (quartets, patterns) = br.get_cpg_quartets_and_patterns();
        for (q, p) in quartets.iter().zip(patterns.iter()) {
            let stat = quartet2stat.entry(*q)
                        .or_insert(PMResult::new(*q));

            stat.add_quartet_pattern(*p);
        }
    }

    for stat in quartet2stat.values() {
        println!("{}", stat.to_bedgraph_field(&header));
    }
}