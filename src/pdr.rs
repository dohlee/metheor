use rust_htslib::{bam::Read, bam::ext::BamRecordExtensions, bam::record::{Aux}};
use std::fs::OpenOptions;
use std::vec::Vec;
use std::str;
use std::fmt;
use std::io::Write;
use std::collections::{BTreeSet, HashMap};
use indicatif::{ProgressBar, ProgressStyle};
use interval_tree::IntervalTree;

use crate::{readutil, bamutil};

struct PDRResult {
    pos: readutil::CpGPosition,
    n_concordant: u32,
    n_discordant: u32,
}

impl PDRResult {
    fn new(pos: readutil::CpGPosition) -> Self {
        Self{ pos: pos, n_concordant: 0, n_discordant: 0 }
    }

    fn inc_concordant(&mut self) {
        self.n_concordant += 1;
    }

    fn inc_discordant(&mut self) {
        self.n_discordant += 1;
    }
}

impl fmt::Display for PDRResult {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let pdr = (self.n_discordant as f32) / (self.n_concordant as f32 + self.n_discordant as f32);
        write!(f, "{}\t{}\t{}\t{}", self.pos, pdr, self.n_concordant, self.n_discordant)
    }
}

pub fn compute(input: &str, output: &str, min_cpgs: i32, min_depth: i32) {
    let mut reader = bamutil::get_reader(&input);
    let header = bamutil::get_header(&reader);

    let mut pdr_result: HashMap<readutil::CpGPosition, PDRResult> = HashMap::new();

    for r in reader.records() {
        let mut r = r.unwrap();
        let br = readutil::BismarkRead::new(&r);

        let mut cpg_positions = br.get_cpg_positions();
        if cpg_positions.len() == 0 { continue; }

        for cpg_position in cpg_positions.iter_mut() {
            let mut r = pdr_result.entry(*cpg_position)
                            .or_insert(PDRResult{ pos:*cpg_position, n_concordant: 0, n_discordant: 0 });
            
            match br.get_concordance_state() {
                readutil::ReadConcordanceState::Concordant => r.inc_concordant(),
                readutil::ReadConcordanceState::Discordant => r.inc_discordant(),
            }
        }
    }

    for (pos, r) in &pdr_result {
        println!("{}", r);
    }
}

pub fn compute_old(input: &str, output: &str, min_cpgs: i32, min_depth: i32) {
    println!("Computing PDR with parameters input={}, output={}, min_depth={}", input, output, min_depth);

    let mut reader = bamutil::get_reader(&input);
    let header = bamutil::get_header(&reader);

    let mut trees: Vec<IntervalTree<i64, i64>> = Vec::new();
    let mut tree = IntervalTree::new();
    let mut cpg_set = BTreeSet::new();  // BTreeSet ensures that the set entries are sorted.
    let mut curr_tid = 0;
    let mut readcount = 0;
    let mut valid_readcount = 0;

    let bar = ProgressBar::new(1);
    bar.set_style(ProgressStyle::default_bar()
        .template("{spinner} {elapsed_precise} {msg}")
    );

    for r in reader.records() {
        let r = r.unwrap();
        let tid = r.tid();
        let start = r.reference_start();
        let end = r.reference_end();

        // Count CpGs and filter.
        match r.aux(b"XM") {
            Ok(value) => {
                if let Aux::String(xm) = value {
                    let n_z = readutil::count_z(xm);

                    if n_z > min_cpgs { // Skip reads with few CpGs.
                        let mut xm_char_set = BTreeSet::new();

                        valid_readcount += 1;
                        for (ref_pos, xm_char) in r.reference_positions_full().zip(xm.chars()) {
                            if let Some(ref_pos) = ref_pos {
                                if xm_char == 'z' || xm_char == 'Z' {
                                    cpg_set.insert((tid, ref_pos));
                                    xm_char_set.insert(xm_char);
                                }
                            }
                        }
                        
                        // Determine if the read is concordant or discordant.
                        if xm_char_set.len() == 1 { // Concordant.
                            tree.insert(start..=end, 0);
                        } else { // Discordant.
                            tree.insert(start..=end, 1);
                        }
                    }
                }
            }
            Err(_) => {
                panic!("Error reading XM tag in BAM record. Make sure the reads are aligned using Bismark!");
            }
        }
        
        if tid != curr_tid {
            curr_tid = tid;
            trees.push(tree);
            tree = IntervalTree::new();
        }

        readcount += 1;
        if readcount % 10000 == 0 {
            bar.inc_length(10000);
            bar.inc(10000);
            bar.set_message(format!("Processed {} reads, found {} valid reads and {} CpGs were covered.", readcount, valid_readcount, cpg_set.len()));
        }
    }
    trees.push(tree);

    bar.finish_with_message(format!("Processed {} reads, found {} valid reads and {} CpGs were covered. Done in {:?}.", readcount, valid_readcount, cpg_set.len(), bar.elapsed()));

    /*
    Done processing.
    */

    // Output file.
    let mut out = OpenOptions::new().create(true).read(true).write(true).open(output).unwrap();

    // let mut csv_reader = csv::Reader::from_path("../data/cpg_quartet_catalog.csv").unwrap();
    let bar = ProgressBar::new(cpg_set.len() as u64);
    bar.set_style(ProgressStyle::default_bar()
        .template("{spinner} {elapsed_precise}>{eta_precise} {bar:30} {msg} ({pos}/{len})")
    );

    // for row in csv_reader.records() {
    for (tid, ref_pos) in cpg_set {
        let tree = trees.get(tid as usize).unwrap();

        let start = ref_pos;
        let end = ref_pos + 1;
        
        let mut n_c: f64 = 0.0;
        let mut n_d: f64 = 0.0;
        for r in tree.find(start..=end) {
            if *r.get() == 1 {
                n_d += 1.0;
            } else {
                n_c += 1.0;
            }
        }

        if n_c + n_d >= (min_depth as f64) { // Skip CpGs covered by insufficient reads.
            let pdr: f64 = n_d / (n_c + n_d);

            let chrom = str::from_utf8(header.tid2name(tid as u32))
                .ok()
                .expect("Error parsing chromosome name.");

            write!(out, "{}\t{}\t{}\t{}\n", chrom, start, end, pdr)
                .ok()
                .expect("Error writing to output file.");
        }

        bar.inc(1);
    }

    bar.finish_with_message(format!("Done in {:?}.", bar.elapsed()));
}