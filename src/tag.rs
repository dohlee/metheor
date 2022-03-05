use rust_htslib::{bam, bam::Read, bam::record::{Aux, Cigar, Record}, bam::ext::BamRecordExtensions};
use rust_htslib::faidx;
use std::collections::HashMap;
use std::str;

use crate::{bamutil};

fn need_reverse_complement(read: &Record) -> bool {
    if (!read.is_reverse() && read.is_first_in_template()) || (read.is_reverse() && read.is_last_in_template()) {
        return false;
    } else {
        return true;
    }
}

fn reverse_complement(seq: &str, switch_base: &HashMap<char, char>) -> String {
    let mut res: String = String::with_capacity(seq.len());
    for base in seq.chars().rev() {
        res.push(switch_base[&base]);
    }

    res
}

fn substring(seq: &str, i: usize, j: usize) -> String {
    seq.chars().skip(i).take(j - i - 1).collect::<String>()
}

fn char_at(seq: &str, i: usize) -> char {
    seq.chars().nth(i).unwrap()
}

fn is_chg_context(seq: &str) -> bool {
    match seq {
        "CAG" | "CTG" | "CCG" => true,
        _ => false
    }
}

fn is_chh_context(seq: &str) -> bool {
    match seq {
        "CAA" | "CAT" | "CAC" | "CTA" | "CTT" | "CTC" | "CCA" | "CCT" | "CCC" => true,
        _ => false
    }
}

fn is_unknown_context(seq: &str) -> bool {
    let mut flag = false;
    for base in seq.chars() {
        if base == '-' || base == 'N' {
            flag = true;
            break;
        }
    }

    flag
}

pub fn run(input: &str, output: &str, genome: &str) {
    let mut reader = bamutil::get_reader(&input);
    let header = bamutil::get_header(&reader);

    let bam = bam::Reader::from_path(&input).unwrap();
    let header_ = bam::Header::from_template(bam.header());

    let mut tid2size: HashMap<usize, usize> = HashMap::new();

    for (key, records) in header_.to_hashmap() {
        if key == "SQ" {
            for (tid, record) in records.iter().enumerate() {
                tid2size.insert(tid, record["LN"].parse().unwrap());
            }
        }
    }

    let mut switch_base: HashMap<char, char> = HashMap::new();
    switch_base.insert('A', 'T');
    switch_base.insert('C', 'G');
    switch_base.insert('G', 'C');
    switch_base.insert('T', 'A');
    switch_base.insert('N', 'N');
    switch_base.insert('-', '-');

    let ref_genome = faidx::Reader::from_path(&genome)
                .ok().expect("Error opening reference genome file");
    let flag_paired_end = bamutil::is_paired_end(&input);

    println!("Parsing reference genome...");
    let mut my_ref_genome: HashMap<usize, &[u8]> = HashMap::new();
    for (tid, _size) in tid2size.iter() {
        let ref_array = ref_genome.fetch_seq(bamutil::tid2chrom(*tid as i32, &header), 0, tid2size[tid])
            .ok().expect("Error fetching reference genome sequence.");
        my_ref_genome.insert(*tid, ref_array);
    }
    println!("Done!");

    let mut writer = bam::Writer::from_path(&output, &header_, bam::Format::Bam).ok().expect("Error opening BAM file to write.");

    for mut r in reader.records().map(|r| r.unwrap()) {
        let tid = r.tid();
        let start = r.reference_start() as usize;
        let end = r.reference_end() as usize;

        let flag_reverse_complement = match flag_paired_end {
            true => need_reverse_complement(&r),
            false => r.is_reverse(),
        };

        let read_seq = str::from_utf8(&r.seq().as_bytes())
            .ok().expect("Error parsing read sequence.").to_string().to_uppercase();

        let ref_seq = str::from_utf8(&my_ref_genome[&(tid as usize)][start - 2..end + 2])
            .ok().expect("Error parsing reference genome.").to_string().to_uppercase();

        let mut tmp_read_seq: Vec<char> = Vec::new();
        let mut tmp_ref_seq: Vec<char> = Vec::new();

        tmp_read_seq.push('-');
        tmp_read_seq.push('-');

        tmp_ref_seq.push(ref_seq.chars().nth(0).unwrap());
        tmp_ref_seq.push(ref_seq.chars().nth(1).unwrap());

        let mut used_read_len: usize = 0;
        let mut used_ref_len: usize = 2;

        for cigar in r.cigar().iter() {
            match cigar {
                Cigar::Match(length) => {

                    tmp_read_seq.append(&mut read_seq.chars().skip(used_read_len).take(*length as usize).collect());
                    tmp_ref_seq.append(&mut ref_seq.chars().skip(used_ref_len).take(*length as usize).collect());

                    used_read_len += *length as usize;
                    used_ref_len += *length as usize;
                },
                Cigar::Ins(length) => {
                    tmp_read_seq.append(&mut read_seq.chars().skip(used_read_len).take(*length as usize).collect());
                    for _ in 0..*length {
                        tmp_ref_seq.push('-');
                    }

                    used_read_len += *length as usize;
                },
                Cigar::Del(length) => {
                    for _ in 0..*length {
                        tmp_read_seq.push('-');
                    }
                    tmp_ref_seq.append(&mut ref_seq.chars().skip(used_ref_len).take(*length as usize).collect());

                    used_ref_len += *length as usize;
                },
                _ => {}
            }
        }

        tmp_read_seq.push('-');
        tmp_read_seq.push('-');

        tmp_ref_seq.push(ref_seq.chars().nth(ref_seq.len() - 2).unwrap());
        tmp_ref_seq.push(ref_seq.chars().nth(ref_seq.len() - 1).unwrap());

        let mut target_read_seq = String::new();
        let mut target_ref_seq = String::new();

        if flag_reverse_complement {
            target_read_seq = tmp_read_seq.iter().take(tmp_read_seq.len() - 2).collect::<String>();
            target_ref_seq = tmp_ref_seq.iter().take(tmp_ref_seq.len() - 2).collect::<String>();

            target_read_seq = reverse_complement(&target_read_seq, &switch_base);
            target_ref_seq = reverse_complement(&target_ref_seq, &switch_base);
        }
        else {
            target_read_seq = tmp_read_seq.iter().skip(2).collect::<String>();
            target_ref_seq = tmp_ref_seq.iter().skip(2).collect::<String>();
        }

        let mut xm_tag: Vec<char> = Vec::new();
        for idx in 0..target_read_seq.len() - 2 {

            if char_at(&target_read_seq, idx) == '-' {
                continue;
            }
            else if char_at(&target_read_seq, idx) == 'N' {
                xm_tag.push('.');
            } else if char_at(&target_ref_seq, idx) == 'C' {
                if (char_at(&target_read_seq, idx+1) == '-' || char_at(&target_read_seq, idx+2) == '-') && ((idx != target_read_seq.len() - 3) && (idx != target_read_seq.len() - 4)) {

                    let mut tmp_target_read_seq: Vec<char> = Vec::new();
                    let mut tmp_target_ref_seq: Vec<char> = Vec::new();

                    tmp_target_read_seq.push(char_at(&target_read_seq, idx));
                    tmp_target_ref_seq.push(char_at(&target_ref_seq, idx));
    
                    let mut flag_tmp = 0;
                    let mut tmp_count = 1;
    
                    while flag_tmp != 2 {
                        if char_at(&target_read_seq, idx + tmp_count) != '-' {
                            tmp_target_read_seq.push(char_at(&target_read_seq, idx + tmp_count));
                            tmp_target_ref_seq.push(char_at(&target_ref_seq, idx + tmp_count));
                            flag_tmp += 1;
                        }
    
                        tmp_count += 1;
                    }
    
                    let ref_context = tmp_target_ref_seq.iter().collect::<String>();

                    if (tmp_target_ref_seq[0] == 'C') && (tmp_target_ref_seq[1] == 'G') {
                        if tmp_target_read_seq[0] == 'C' { xm_tag.push('Z'); }
                        else if tmp_target_read_seq[0] == 'T' { xm_tag.push('z'); }
                        else { xm_tag.push('.'); }
                    }

                    else if is_chg_context(&ref_context) {
                        if tmp_target_read_seq[0] == 'C' { xm_tag.push('X'); }
                        else if tmp_target_read_seq[0] == 'T' { xm_tag.push('x'); }
                        else { xm_tag.push('.'); }
                    }

                    else if is_chh_context(&ref_context) {
                        if tmp_target_read_seq[0] == 'C' { xm_tag.push('H'); }
                        else if tmp_target_read_seq[0] == 'T' { xm_tag.push('h'); }
                        else { xm_tag.push('.'); }
                    }

                    else if is_unknown_context(&ref_context) {
                        if tmp_target_read_seq[0] == 'C' { xm_tag.push('U'); }
                        else if tmp_target_read_seq[0] == 'T' { xm_tag.push('u'); }
                        else { xm_tag.push('.'); }
                    }
                }
                // No deletion
                else {
                    let ref_context = target_ref_seq.chars().skip(idx).take(3).collect::<String>();

                    if (char_at(&target_ref_seq, idx) == 'C') && (char_at(&target_ref_seq, idx+1) == 'G') {
                        if char_at(&target_read_seq, idx) == 'C' { xm_tag.push('Z'); }
                        else if char_at(&target_read_seq, idx) == 'T' { xm_tag.push('z'); }
                        else { xm_tag.push('.'); }
                    }

                    else if is_chg_context(&ref_context) {
                        if char_at(&target_read_seq, idx) == 'C' { xm_tag.push('X'); }
                        else if char_at(&target_read_seq, idx) == 'T' { xm_tag.push('x'); }
                        else { xm_tag.push('.'); }
                    }

                    else if is_chh_context(&ref_context) {
                        if char_at(&target_read_seq, idx) == 'C' { xm_tag.push('H'); }
                        else if char_at(&target_read_seq, idx) == 'T' { xm_tag.push('h'); }
                        else { xm_tag.push('.'); }
                    }

                    else if is_unknown_context(&ref_context) {
                        if char_at(&target_read_seq, idx) == 'C' { xm_tag.push('U'); }
                        else if char_at(&target_read_seq, idx) == 'T' { xm_tag.push('u'); }
                        else { xm_tag.push('.'); }
                    }
                }
            } else {
                xm_tag.push('.');
            }
        }

        let xm_tag_string = match flag_reverse_complement {
            true => xm_tag.iter().rev().collect::<String>(),
            false => xm_tag.iter().collect::<String>()
        };

        r.push_aux("XM".as_bytes(), Aux::String(&xm_tag_string));
        writer.write(&r).ok().expect("Error writing to output file.");
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use super::super::bamutil;
    use rust_htslib::bam::Read;

    #[test]
    fn test_reverse_complement() {
        let mut switch_base: HashMap<char, char> = HashMap::new();
        switch_base.insert('A', 'T');
        switch_base.insert('C', 'G');
        switch_base.insert('G', 'C');
        switch_base.insert('T', 'A');
        switch_base.insert('N', 'N');
        switch_base.insert('-', '-');

        assert_eq!("ATGC", reverse_complement("GCAT", &switch_base));
    }
}