use rust_htslib::{bam::Read, bam::ext::BamRecordExtensions, bam::record::{Aux}};
use bio_types::genome::AbstractInterval;
use std::fs;
use std::vec::Vec;
use std::str;

use std::collections::{HashSet, HashMap};
use indicatif::{ProgressBar, ProgressStyle};

use crate::{readutil, bamutil};


pub fn compute(input: &str, output: &str, min_distance: i32, max_distance: i32, cpg_set: &str, min_qual: u8) -> (i32, i32) {

    if cpg_set.is_empty() {
        let (n_concordant, n_discordant) = compute_all(input, output, min_distance, max_distance, min_qual);
        (n_concordant, n_discordant)
    } else {
        let (n_concordant, n_discordant) = compute_subset(input, output, min_distance, max_distance, cpg_set, min_qual);
        (n_concordant, n_discordant)
    }
}

fn compute_all(input: &str, output: &str, min_distance: i32, max_distance: i32, min_qual: u8) -> (i32, i32) {
    println!("Computing LPMD with parameters input={}, output={}, min_distance={}, max_distance={}", input, output, min_distance, max_distance);

    let mut reader = bamutil::get_reader(&input);
    
    let mut readcount = 0;
    let mut valid_readcount = 0;
    let mut n_concordant = 0;
    let mut n_discordant = 0;

    let bar = ProgressBar::new(1);
    bar.set_style(ProgressStyle::default_bar()
        .template("{spinner} {elapsed_precise} {msg}")
    );

    for r in reader.records() {
        let r = r.unwrap();

        readcount += 1;
        if r.mapq() < min_qual { continue; }
        
        match r.aux(b"XM") {
            Ok(value) => {
                if let Aux::String(xm) = value {
                    valid_readcount += 1;


                    let mut cpgs: Vec<(i32, char)> = Vec::new();
                    for (pos, c) in r.range().zip(xm.chars()) {
                        
                        if (c != 'z') && (c != 'Z') {
                            continue;
                        }
                        
                        if (r.flags() == 99) || (r.flags() == 147) { // Forward
                            cpgs.push((pos as i32, c));
                        } else {
                            cpgs.push(((pos - 1)as i32, c));
                        }
                    }

                    let (c, d) = readutil::compute_pairwise_concordance_discordance(cpgs, min_distance, max_distance);

                    n_concordant += c;
                    n_discordant += d;
                }
            }
            Err(_) => {
                panic!("Error reading XM tag in BAM record. Make sure the reads are aligned using Bismark!");
            }
        }
    
        readcount += 1;
        if readcount % 10000 == 0 {
            bar.inc_length(10000);
            bar.inc(10000);
            bar.set_message(format!("Processed {} reads, found {} valid reads.", readcount, valid_readcount));
        }
    }

    println!("{} {} {}", n_concordant, n_discordant, n_discordant / (n_concordant + n_discordant));
    println!("{} {} {}", n_concordant, n_discordant, n_discordant / (n_concordant + n_discordant));
    println!("{} valid reads.", valid_readcount);

    (n_concordant, n_discordant)
}

fn compute_subset(input: &str, output: &str, min_distance: i32, max_distance: i32, cpg_set: &str, min_qual: u8) -> (i32, i32) {
    println!("Computing subset-LPMD with parameters input={}, output={}, min_distance={}, max_distance={}", input, output, min_distance, max_distance);

    println!("Processing target CpG set...");
    let mut target_cpgs: HashSet<(&str, i32)> = HashSet::new();

    let contents = fs::read_to_string(cpg_set)
                    .expect("Could not read target CpG file.");

    println!("{}", cpg_set);

    for line in contents.lines() {
        let tokens: Vec<&str> = line.split("\t").collect();

        let chrom = tokens[0];
        let pos = tokens[1].parse::<i32>().unwrap();
        
        target_cpgs.insert((chrom, pos));
    }

    println!("Processed {} CpGs in total.", target_cpgs.len());


    let mut reader = bamutil::get_reader(&input);
    let header = bamutil::get_header(&reader);
    
    let mut readcount = 0;
    let mut valid_readcount = 0;
    let mut n_forward = 0;
    let mut n_reverse = 0;
    let mut n_concordant = 0;
    let mut n_discordant = 0;

    let mut flag_counter: HashMap<u16, i32> = HashMap::new();

    let bar = ProgressBar::new(1);
    bar.set_style(ProgressStyle::default_bar()
        .template("{spinner} {elapsed_precise} {msg}")
    );

    for r in reader.records() {
        let mut r = r.unwrap();

        readcount += 1;
        if r.mapq() < min_qual { continue; }


        r.cache_cigar();

        let tid: i32 = r.tid();
        let chrom = str::from_utf8(header.tid2name(tid as u32))
            .ok()
            .expect("Error parsing chromosome name.");

        if (r.flags() == 99) || (r.flags() == 147) { // Forward
            n_forward += 1;
        } else {
            n_reverse += 1;
        }

        let v = flag_counter.entry(r.flags()).or_insert(0);
        *v += 1;
        
        match r.aux(b"XM") {
            Ok(value) => {
                if let Aux::String(xm) = value {

                    // println!("{}", r.flags());

                    let mut cpgs: Vec<(i32, char)> = Vec::new();
                    for (pos, c) in r.reference_positions_full().zip(xm.chars()) {
                        
                        if (c != 'z') && (c != 'Z') {
                            continue;
                        }

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

                    valid_readcount += 1;

                    let (c, d) = readutil::compute_pairwise_concordance_discordance(cpgs, min_distance, max_distance);

                    n_concordant += c;
                    n_discordant += d;
                }
            }
            Err(_) => {
                panic!("Error reading XM tag in BAM record. Make sure the reads are aligned using Bismark!");
            }
        }
    
        if readcount % 10000 == 0 {
            bar.inc_length(10000);
            bar.inc(10000);
            bar.set_message(format!("Processed {} reads, found {} valid reads. ({}, {})", readcount, valid_readcount, n_concordant, n_discordant));
        }
    }

    println!("{} {}", n_concordant, n_discordant);
    println!("{} {}", n_concordant, n_discordant);
    println!("{} valid reads.", valid_readcount);
    println!("{} {} (F, R)", n_forward, n_reverse);
    println!("{} valid reads.", valid_readcount);

    println!("{:?}", flag_counter);
    println!("{:?}", flag_counter);

    (n_concordant, n_discordant)
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
