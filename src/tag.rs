use rust_htslib::faidx;
use rust_htslib::{
    bam,
    bam::ext::BamRecordExtensions,
    bam::record::{Aux, Cigar, Record},
    bam::Read,
};
use std::cmp::{max, min};
use std::collections::HashMap;
use std::path::PathBuf;
use std::str;

use crate::bamutil;

fn need_reverse_complement(read: &Record) -> bool {
    if (!read.is_reverse() && read.is_first_in_template())
        || (read.is_reverse() && read.is_last_in_template())
    {
        return false;
    } else {
        return true;
    }
}

fn reverse_complement(seq: &str, rcmapping: &HashMap<char, char>) -> String {
    let mut res: String = String::with_capacity(seq.len());
    for base in seq.chars().rev() {
        res.push(rcmapping[&base]);
    }
    res
}

fn char_at(seq: &str, i: usize) -> char {
    seq.chars().nth(i).unwrap()
}

fn is_chg_context(seq: &str) -> bool {
    match seq {
        "CAG" | "CTG" | "CCG" => true,
        _ => false,
    }
}

fn is_chh_context(seq: &str) -> bool {
    match seq {
        "CAA" | "CAT" | "CAC" | "CTA" | "CTT" | "CTC" | "CCA" | "CCT" | "CCC" => true,
        _ => false,
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

pub fn get_header_template_from_bam(input: &str) -> bam::Header {
    let bam = bam::Reader::from_path(&input).unwrap();
    bam::Header::from_template(bam.header())
}

pub fn get_tid2size_from_bam(input: &str) -> HashMap<usize, usize> {
    let header_ = get_header_template_from_bam(input);

    let mut tid2size: HashMap<usize, usize> = HashMap::new();
    for (key, records) in header_.to_hashmap() {
        // SQ = Reference sequence dictionary.
        // The order of @SQ lines defines the alignment sorting order.
        // Each @SQ entry has "SN" (Reference sequence name) and "LN" (Reference sequence length) field.
        // https://samtools.github.io/hts-specs/SAMv1.pdf
        if key != "SQ" {
            continue;
        }

        for (tid, record) in records.iter().enumerate() {
            tid2size.insert(tid, record["LN"].parse().unwrap());
        }
    }
    tid2size
}

pub fn get_rcmapping() -> HashMap<char, char> {
    // reverse-complement mapping
    let mut rcmapping: HashMap<char, char> = HashMap::new();
    rcmapping.insert('A', 'T');
    rcmapping.insert('C', 'G');
    rcmapping.insert('G', 'C');
    rcmapping.insert('T', 'A');
    rcmapping.insert('N', 'N');
    rcmapping.insert('M', 'K');
    rcmapping.insert('R', 'Y');
    rcmapping.insert('W', 'W');
    rcmapping.insert('S', 'S');
    rcmapping.insert('Y', 'R');
    rcmapping.insert('K', 'M');
    rcmapping.insert('V', 'B');
    rcmapping.insert('H', 'D');
    rcmapping.insert('D', 'H');
    rcmapping.insert('B', 'V');
    rcmapping.insert('-', '-');

    rcmapping
}

// fn process_match(
//     tmp_read_seq: Vec<char>,
//     tmp_ref_seq: Vec<char>,
//     read_seq: String,
//     ref_seq: String,
//     used_read_len: &mut usize,
//     used_ref_len: &mut usize,
//     length: &mut u32,
// ) -> () {
//     tmp_read_seq.append(
//         &mut read_seq
//             .chars()
//             .skip(used_read_len)
//             .take(*length as usize)
//             .collect(),
//     );
//     tmp_ref_seq.append(
//         &mut ref_seq
//             .chars()
//             .skip(used_ref_len)
//             .take(*length as usize)
//             .collect(),
//     );

//     *used_read_len += *length as usize;
//     *used_ref_len += *length as usize;
// }

pub fn determine_xm_tag_string(
    r: &Record,
    refgenome: &HashMap<usize, &[u8]>,
    tid2size: &HashMap<usize, usize>,
    rcmapping: &HashMap<char, char>,
    is_paired_end: bool,
) -> String {
    let tid = r.tid();
    let start = r.reference_start();
    let end = r.reference_end();

    let flag_reverse_complement = match is_paired_end {
        true => need_reverse_complement(&r),
        false => r.is_reverse(),
    };

    // Extract read sequence from alignment record.
    let read_seq = match str::from_utf8(&r.seq().as_bytes()) {
        Ok(read_seq) => read_seq.to_string().to_uppercase(),
        Err(error) => panic!("Error parsing alignment record: {}", error),
    };
    // For reference sequence,
    // we should additionally consider upstream & downstream 2-bp positions,
    // to determine the cytosine context near the left & right edge of the alignment.
    let chromsize = tid2size[&(tid as usize)] as i64;
    let clipped_start = max(start - 2, 0) as usize;
    let clipped_end = min(end + 2, chromsize) as usize;

    let ref_seq_result = str::from_utf8(&refgenome[&(tid as usize)][clipped_start..clipped_end]);
    let ref_seq = match ref_seq_result {
        Ok(ref_seq) => ref_seq.to_string().to_uppercase(),
        Err(error) => panic!("Error extracting reference sequence: {}", error),
    };
    // For reads aligned at the edge of the reference genome,
    // we may not be able to extract flanking 2bp. In that case, just pad with N as much as needed.
    let padding = ["", "N", "NN"];
    let pad_nbases_start = max(2 - start, 0) as usize;
    let pad_nbases_end = max(end - chromsize + 2, 0) as usize;
    let ref_seq = format!(
        "{}{}{}",
        padding[pad_nbases_start], ref_seq, padding[pad_nbases_end]
    );

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
                tmp_read_seq.append(
                    &mut read_seq
                        .chars()
                        .skip(used_read_len)
                        .take(*length as usize)
                        .collect(),
                );
                tmp_ref_seq.append(
                    &mut ref_seq
                        .chars()
                        .skip(used_ref_len)
                        .take(*length as usize)
                        .collect(),
                );

                used_read_len += *length as usize;
                used_ref_len += *length as usize;
            }
            Cigar::Ins(length) => {
                tmp_read_seq.append(
                    &mut read_seq
                        .chars()
                        .skip(used_read_len)
                        .take(*length as usize)
                        .collect(),
                );
                for _ in 0..*length {
                    tmp_ref_seq.push('-');
                }

                used_read_len += *length as usize;
            }
            Cigar::Del(length) => {
                for _ in 0..*length {
                    tmp_read_seq.push('-');
                }
                tmp_ref_seq.append(
                    &mut ref_seq
                        .chars()
                        .skip(used_ref_len)
                        .take(*length as usize)
                        .collect(),
                );

                used_ref_len += *length as usize;
            }
            _ => {}
        }
    }

    tmp_read_seq.push('-');
    tmp_read_seq.push('-');

    tmp_ref_seq.push(ref_seq.chars().nth(ref_seq.len() - 2).unwrap());
    tmp_ref_seq.push(ref_seq.chars().nth(ref_seq.len() - 1).unwrap());

    let target_read_seq;
    let target_ref_seq;

    if flag_reverse_complement {
        let read = tmp_read_seq
            .iter()
            .take(tmp_read_seq.len() - 2)
            .collect::<String>();
        let reference = tmp_ref_seq
            .iter()
            .take(tmp_ref_seq.len() - 2)
            .collect::<String>();

        target_read_seq = reverse_complement(&read, &rcmapping);
        target_ref_seq = reverse_complement(&reference, &rcmapping);
    } else {
        target_read_seq = tmp_read_seq.iter().skip(2).collect::<String>();
        target_ref_seq = tmp_ref_seq.iter().skip(2).collect::<String>();
    }

    let mut xm_tag: Vec<char> = Vec::new();
    for idx in 0..target_read_seq.len() - 2 {
        if char_at(&target_read_seq, idx) == '-' {
            continue;
        } else if char_at(&target_read_seq, idx) == 'N' {
            xm_tag.push('.');
        } else if char_at(&target_ref_seq, idx) == 'C' {
            if (char_at(&target_read_seq, idx + 1) == '-'
                || char_at(&target_read_seq, idx + 2) == '-')
                && ((idx != target_read_seq.len() - 3) && (idx != target_read_seq.len() - 4))
            {
                let mut tmp_target_read_seq: Vec<char> = Vec::new();
                let mut tmp_target_ref_seq: Vec<char> = Vec::new();

                tmp_target_read_seq.push(char_at(&target_read_seq, idx));
                tmp_target_ref_seq.push(char_at(&target_ref_seq, idx));

                let mut flag_tmp = 0;
                let mut tmp_count = 1;

                while flag_tmp != 2 {
                    if idx + tmp_count > target_read_seq.len() - 1 {
                        break;
                    }
                    if char_at(&target_read_seq, idx + tmp_count) != '-' {
                        tmp_target_read_seq.push(char_at(&target_read_seq, idx + tmp_count));
                        tmp_target_ref_seq.push(char_at(&target_ref_seq, idx + tmp_count));
                        flag_tmp += 1;
                    }

                    tmp_count += 1;
                }

                let ref_context = tmp_target_ref_seq.iter().collect::<String>();

                if (tmp_target_ref_seq[0] == 'C') && (tmp_target_ref_seq[1] == 'G') {
                    if tmp_target_read_seq[0] == 'C' {
                        xm_tag.push('Z');
                    } else if tmp_target_read_seq[0] == 'T' {
                        xm_tag.push('z');
                    } else {
                        xm_tag.push('.');
                    }
                } else if is_chg_context(&ref_context) {
                    if tmp_target_read_seq[0] == 'C' {
                        xm_tag.push('X');
                    } else if tmp_target_read_seq[0] == 'T' {
                        xm_tag.push('x');
                    } else {
                        xm_tag.push('.');
                    }
                } else if is_chh_context(&ref_context) {
                    if tmp_target_read_seq[0] == 'C' {
                        xm_tag.push('H');
                    } else if tmp_target_read_seq[0] == 'T' {
                        xm_tag.push('h');
                    } else {
                        xm_tag.push('.');
                    }
                } else if is_unknown_context(&ref_context) {
                    if tmp_target_read_seq[0] == 'C' {
                        xm_tag.push('U');
                    } else if tmp_target_read_seq[0] == 'T' {
                        xm_tag.push('u');
                    } else {
                        xm_tag.push('.');
                    }
                }
            }
            // No deletion
            else {
                let ref_context = target_ref_seq.chars().skip(idx).take(3).collect::<String>();

                if (char_at(&target_ref_seq, idx) == 'C')
                    && (char_at(&target_ref_seq, idx + 1) == 'G')
                {
                    // Reference context is CG, read 'C' -> 'methylated in CG context (Z)'
                    if char_at(&target_read_seq, idx) == 'C' {
                        xm_tag.push('Z');
                    }
                    // Reference context is CG, read 'T' -> 'unmethylated in CG context (z)'
                    else if char_at(&target_read_seq, idx) == 'T' {
                        xm_tag.push('z');
                    }
                    // Reference context is CG, read 'A or G' -> Nothing.
                    else {
                        xm_tag.push('.');
                    }
                } else if is_chg_context(&ref_context) {
                    if char_at(&target_read_seq, idx) == 'C' {
                        xm_tag.push('X');
                    } else if char_at(&target_read_seq, idx) == 'T' {
                        xm_tag.push('x');
                    } else {
                        xm_tag.push('.');
                    }
                } else if is_chh_context(&ref_context) {
                    if char_at(&target_read_seq, idx) == 'C' {
                        xm_tag.push('H');
                    } else if char_at(&target_read_seq, idx) == 'T' {
                        xm_tag.push('h');
                    } else {
                        xm_tag.push('.');
                    }
                } else if is_unknown_context(&ref_context) {
                    if char_at(&target_read_seq, idx) == 'C' {
                        xm_tag.push('U');
                    } else if char_at(&target_read_seq, idx) == 'T' {
                        xm_tag.push('u');
                    } else {
                        xm_tag.push('.');
                    }
                }
            }
        } else {
            xm_tag.push('.');
        }
    }

    match flag_reverse_complement {
        true => xm_tag.iter().rev().collect::<String>(),
        false => xm_tag.iter().collect::<String>(),
    }
}

pub fn run(input: &str, output: &str, genome: &str) {
    let mut reader = bamutil::get_reader(&input);
    let is_paired_end = bamutil::is_paired_end(&input);
    let header = bamutil::get_header(&reader);
    let tid2size: HashMap<usize, usize> = get_tid2size_from_bam(&input);

    let rcmapping = get_rcmapping();

    // Assert if the output directory exists.
    let path = PathBuf::from(&output);
    let dir = path.parent().unwrap();

    if !dir.is_dir() {
        panic!(
            "No such directory for output alignment file: {}",
            dir.to_str().unwrap()
        )
    }
    // Prepare output writer.
    let header_tmpl = get_header_template_from_bam(&input);
    let mut writer = match bam::Writer::from_path(&output, &header_tmpl, bam::Format::Sam) {
        Ok(writer) => writer,
        Err(error) => panic!("Error opening alignment file to write: {}", error),
    };
    // Prepare reference genome.
    let refgenome_reader = match faidx::Reader::from_path(&genome) {
        Ok(refgenome_reader) => refgenome_reader,
        Err(error) => {
            panic!("Error opening reference genome file: {}", error);
        }
    };
    println!("Parsing reference genome...");
    let mut refgenome: HashMap<usize, &[u8]> = HashMap::new();
    for (tid, _size) in tid2size.iter() {
        let ref_array = refgenome_reader
            .fetch_seq(bamutil::tid2chrom(*tid as i32, &header), 0, tid2size[tid])
            .expect("Error fetching reference genome sequence.");

        refgenome.insert(*tid, ref_array);
    }
    println!("Done!");

    // Main loop
    // Iterate aligned reads and determine xm tag string.
    for mut r in reader.records().map(|r| r.unwrap()) {
        // Determine XM tag string by comparing read sequence and reference sequence.
        let xm_tag_string =
            determine_xm_tag_string(&r, &refgenome, &tid2size, &rcmapping, is_paired_end);
        // Attach XM tag to the record.
        let add_result = r.push_aux("XM".as_bytes(), Aux::String(&xm_tag_string));
        match add_result {
            Ok(_) => (),
            Err(e) => panic!("Error adding XM tag to alignment record. {}", e),
        }
        // Write record to output.
        writer
            .write(&r)
            .ok()
            .expect("Error writing to output file.");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        let rcmapping = get_rcmapping();
        assert_eq!("ATGC", reverse_complement("GCAT", &rcmapping));
    }
    #[test]
    #[should_panic]
    fn error_when_output_directory_is_not_found() {
        run(
            "tests/test1.bam",
            "tests/no_such_directory/out.bam",
            "tests/tinyref.fa",
        )
    }
    #[test]
    #[should_panic]
    fn error_when_reference_genome_is_not_found() {
        run(
            "tests/test1.bam",
            "tests/out.tagged.bam",
            "tests/there_is_no_such.fa",
        )
    }
}
