mod common;

use common::Workspace;
use rust_htslib::{bam, bam::Read};
use std::{
    convert::TryFrom,
    path::{Path, PathBuf},
};

fn bam_copy(w: &Workspace, input: &Path) -> PathBuf {
    let mut reader = bam::Reader::from_path(input).unwrap();
    let output = w.path("input.bam");
    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(&output, &header, bam::Format::Bam).unwrap();
    for record in reader.records() {
        writer.write(&record.unwrap()).unwrap();
    }
    output
}

#[test]
fn reference_skips_advance_context_without_adding_xm_bases() {
    for (flags, expected) in [(0, "Z.Z."), (16, ".Z.Z")] {
        let w = Workspace::new();
        let reference = w.write("reference.fa", ">chr1\nCAGCGTCC\n");
        rust_htslib::faidx::build(&reference).unwrap();
        let sam = w.write(
            "input.sam",
            &format!("@SQ\tSN:chr1\tLN:8\nr\t{flags}\tchr1\t1\t60\t1M1N3M\t*\t0\t0\tCGCG\t*\n"),
        );
        for input in [sam.clone(), bam_copy(&w, &sam)] {
            w.command("tag", &input, &["-g", reference.to_str().unwrap()])
                .assert()
                .success();
            let mut reader = bam::Reader::from_path(w.path("output.tsv")).unwrap();
            let records: Vec<_> = reader.records().map(Result::unwrap).collect();
            assert_eq!(records.len(), 1);
            assert_eq!(records[0].flags(), flags);
            assert_eq!(records[0].cigar().to_string(), "1M1N3M");
            assert_eq!(records[0].seq().as_bytes(), b"CGCG");
            assert_eq!(
                records[0].aux(b"XM").unwrap(),
                bam::record::Aux::String(expected)
            );
        }
    }
}

#[test]
fn reported_complex_cigars_match_bismark_calls_and_preserve_alignment_data() {
    use bam::record::{Aux, Cigar, CigarString};

    let w = Workspace::new();
    let genome = "ACGTCAGTCACTGAAG".repeat(1250);
    let reference = w.write("reference.fa", &format!(">chr1\n{genome}\n"));
    rust_htslib::faidx::build(&reference).unwrap();
    let cases: Vec<_> = include_str!("tag_cigar_bismark.tsv")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .collect();
    assert_eq!(cases.len(), 10);
    let header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20000\n@RG\tID:sample\tSM:test\n";
    let mut source = header.to_owned();
    for case in &cases {
        assert_eq!(case.len(), 5);
        let name = case[0];
        let flags: u16 = case[1].parse().unwrap();
        let pos: usize = case[2].parse().unwrap();
        let cigar = CigarString::try_from(case[3]).unwrap();
        let mut reference_pos = pos - 1;
        let mut sequence = Vec::new();
        for operation in &cigar {
            match operation {
                Cigar::Match(length) | Cigar::Equal(length) | Cigar::Diff(length) => {
                    let end = reference_pos + *length as usize;
                    sequence.extend_from_slice(&genome.as_bytes()[reference_pos..end]);
                    reference_pos = end;
                }
                Cigar::Ins(length) | Cigar::SoftClip(length) => {
                    sequence.extend(std::iter::repeat_n(b'A', *length as usize));
                }
                Cigar::Del(length) | Cigar::RefSkip(length) => reference_pos += *length as usize,
                _ => {}
            }
        }
        let (original, converted) = if flags == 16 {
            (b'G', b'A')
        } else {
            (b'C', b'T')
        };
        for (index, base) in sequence.iter_mut().enumerate() {
            if index % 3 == 0 && *base == original {
                *base = converted;
            }
        }
        let sequence = std::str::from_utf8(&sequence).unwrap();
        source.push_str(&format!(
            "{name}\t{flags}\tchr1\t{pos}\t60\t{cigar}\t*\t0\t0\t{sequence}\t*\tRG:Z:sample\n"
        ));
    }
    let sam = w.write("input.sam", &source);
    let mut reader = bam::Reader::from_path(&sam).unwrap();
    let originals: Vec<_> = reader.records().map(Result::unwrap).collect();
    for input in [sam.clone(), bam_copy(&w, &sam)] {
        w.command("tag", &input, &["-g", reference.to_str().unwrap()])
            .assert()
            .success();
        let mut reader = bam::Reader::from_path(w.path("output.tsv")).unwrap();
        assert_eq!(reader.header().as_bytes(), header.as_bytes());
        let records: Vec<_> = reader.records().map(Result::unwrap).collect();
        assert_eq!(records.len(), cases.len());
        for ((record, original), case) in records.iter().zip(&originals).zip(&cases) {
            let name = case[0];
            assert_eq!(record.qname(), original.qname());
            assert_eq!(record.flags(), original.flags());
            assert_eq!(record.tid(), original.tid());
            assert_eq!(record.pos(), original.pos());
            assert_eq!(record.mapq(), original.mapq());
            assert_eq!(record.cigar().to_string(), original.cigar().to_string());
            assert!(
                record.seq().as_bytes() == original.seq().as_bytes(),
                "{}: sequence changed",
                name
            );
            assert!(
                record.qual() == original.qual(),
                "{}: quality changed",
                name
            );
            assert_eq!(record.aux(b"RG").unwrap(), Aux::String("sample"));
            let actual = match record.aux(b"XM").unwrap() {
                Aux::String(xm) => xm,
                _ => panic!("{}: XM is not a string", name),
            };
            assert_eq!(actual.len(), original.seq_len(), "{name}");
            assert_eq!(actual.len(), case[4].len(), "{name}");
            for (index, (actual, expected)) in actual.bytes().zip(case[4].bytes()).enumerate() {
                assert_eq!(actual, expected, "{name}: XM at query offset {index}");
            }
        }
    }
}
