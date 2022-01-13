use rust_htslib::{bam::Read, bam::ext::BamRecordExtensions, bam::record::{Aux}};
use std::fs::OpenOptions;
use std::vec::Vec;
use std::str;
use std::io::Write;
use std::collections::{BTreeMap};
use indicatif::{ProgressBar, ProgressStyle};
use interval_tree::IntervalTree;

use crate::{readutil, bamutil};

pub fn compute(input: &str, output: &str, min_depth: i32) {
    println!("Running compute_pm with parameters input={}, output={}, min_depth={}", input, output, min_depth);

    let mut reader = bamutil::get_reader(&input);
    let header = bamutil::get_header(&reader);

    let mut trees: Vec<IntervalTree<i64, i64>> = Vec::new();
    let mut tree = IntervalTree::new();
    let mut quartet_map: BTreeMap<(i32, i64, i64, i64, i64), Vec<i64>> = BTreeMap::new();
    let mut curr_tid = 0;
    let mut readcount = 0;
    let mut valid_readcount = 0;

    let bar = ProgressBar::new(1);
    bar.set_style(ProgressStyle::default_bar()
        .template("{spinner} {elapsed_precise} {msg}")
    );

    for r in reader.records() {
        let r = r.unwrap();
        let _tid: i32 = r.tid();
        let _start: i64 = r.reference_start();
        let _end: i64 = r.reference_end();

        let tid = r.tid();
        
    
        match r.aux(b"XM") {
            Ok(value) => {
                if let Aux::String(xm) = value {
                    let mut ref_positions: Vec<i64> = Vec::new();
                    let mut meth_states: Vec<i64> = Vec::new();
                    let _n_z = readutil::count_z(xm);
    
                    for (ref_pos, xm_char) in r.reference_positions_full().zip(xm.chars()) {
                        if let Some(ref_pos) = ref_pos {
                            if xm_char == 'z' || xm_char == 'Z' {
                                ref_positions.push(ref_pos);

                                if xm_char == 'Z' {
                                    meth_states.push(1);
                                } else {
                                    meth_states.push(0);
                                }
                            }
                        }
                    }

                    if ref_positions.len() >= 4 {
                        valid_readcount += 1;

                        for i in 0..ref_positions.len() - 3 {
                            let c1 = ref_positions[i];
                            let c2 = ref_positions[i + 1];
                            let c3 = ref_positions[i + 2];
                            let c4 = ref_positions[i + 3];

                            let m1 = meth_states[i];
                            let m2 = meth_states[i + 1];
                            let m3 = meth_states[i + 2];
                            let m4 = meth_states[i + 3];

                            let quartet = (tid, c1, c2, c3, c4);
                            let this_pat = 8 * m1 + 4 * m2 + 2 * m3 + m4;

                            if quartet_map.contains_key(&quartet) {
                                let mut vec = quartet_map.get(&quartet).unwrap().clone();
                                vec.push(this_pat);
                                quartet_map.insert(quartet, vec);
                            } else {
                                let mut vec: Vec<i64> = Vec::new();
                                vec.push(this_pat);
                                quartet_map.insert(quartet, vec);
                            }
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
            bar.set_message(format!("Processed {} reads, found {} valid reads and {} quartets were covered.", readcount, valid_readcount, quartet_map.len()));
        }
    }
    trees.push(tree);
    bar.finish_with_message(format!("Processed {} reads, found {} valid reads and {} quartets were covered. Done in {:?}.", readcount, valid_readcount, quartet_map.len(), bar.elapsed()));

    /*
    Done processing.
    */

    // Output file.
    let mut out = OpenOptions::new().create(true).read(true).write(true).open(output).unwrap();

    // let mut csv_reader = csv::Reader::from_path("../data/cpg_quartet_catalog.csv").unwrap();
    let bar = ProgressBar::new(quartet_map.len() as u64);
    bar.set_style(ProgressStyle::default_bar()
        .template("{spinner} {elapsed_precise}>{eta_precise} {bar:30} {msg} ({pos}/{len})")
    );

    // for row in csv_reader.records() {
    for ((tid, c1, c2, c3, c4), pats) in quartet_map {
        let _tree = trees.get(tid as usize).unwrap();
        let _quartet = (tid, c1, c2, c3, c4);
        let n_pat_total: f64 = pats.len() as f64;

        if n_pat_total >= min_depth as f64 { // Skip CpGs covered by insufficient reads.
            let mut pm: f64 = 1.0;
            // Count patterns.
            let mut c: [f64; 16] = [0.0; 16];
            for pat in pats {
                c[pat as usize] += 1.0;
            }

            for n_pat in c {
                pm -= (n_pat / n_pat_total).powi(2);
            }

            let chrom = str::from_utf8(header.tid2name(tid as u32))
                .ok()
                .expect("Error parsing chromosome name.");

            write!(out, "{}\t{}\t{}\t{}\t{}\t{}\n", chrom, c1, c2, c3, c4, pm)
                .ok()
                .expect("Error writing to output file.");
        }

        bar.inc(1);
    }

    bar.finish_with_message(format!("Done in {:?}.", bar.elapsed()));
}