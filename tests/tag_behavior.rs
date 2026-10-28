mod common;

use common::Workspace;
use rust_htslib::{bam, bam::Read};
use std::fs;

fn reference(w: &Workspace) -> std::path::PathBuf {
    let path = w.write("reference.fa", ">chr1\nCGCGCGCGCGCGCGCGCGCG\n");
    rust_htslib::faidx::build(&path).unwrap();
    path
}

fn input(w: &Workspace, flags: u16, cigar: &str, seq: &str) -> std::path::PathBuf {
    w.write("input.sam", &format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\nr\t{flags}\tchr1\t1\t60\t{cigar}\t*\t0\t0\t{seq}\t*\n"))
}

#[test]
fn missing_reference_preserves_existing_output() {
    let w = Workspace::new();
    let input = input(&w, 0, "4M", "CGCG");
    let output = w.write("output.tsv", "previous successful result\n");
    w.command(
        "tag",
        &input,
        &["-g", w.path("missing.fa").to_str().unwrap()],
    )
    .assert()
    .failure();
    assert_eq!(
        fs::read_to_string(output).unwrap(),
        "previous successful result\n"
    );
}

#[test]
fn xm_keeps_soft_clipped_bases_and_reference_edge_context() {
    for (cigar, seq, expected) in [
        ("4M", "CGCG", "Z.Z."),
        ("2S4M", "AACGCG", "..Z.Z."),
        ("4M2S", "CGCGAA", "Z.Z..."),
        ("1M2S", "CAA", "Z.."),
        ("2H4M1H", "CGCG", "Z.Z."),
        ("2M2I2M", "CGAACG", "Z...Z."),
        ("2M2D2M", "CGCG", "Z.Z."),
    ] {
        let w = Workspace::new();
        let input = input(&w, 0, cigar, seq);
        let reference = reference(&w);
        w.command("tag", &input, &["-g", reference.to_str().unwrap()])
            .assert()
            .success();
        let mut reader = bam::Reader::from_path(w.path("output.tsv")).unwrap();
        let records: Vec<_> = reader.records().map(Result::unwrap).collect();
        assert_eq!(records.len(), 1);
        assert_eq!(
            records[0].aux(b"XM").unwrap(),
            bam::record::Aux::String(expected)
        );
        assert_eq!(records[0].seq_len(), expected.len());
    }
}

#[test]
fn xm_supports_equal_and_mismatch_cigars() {
    for (cigar, seq, expected) in [("4=", "CGCG", "Z.Z."), ("2=1X1=", "CGTG", "Z.z.")] {
        let w = Workspace::new();
        let input = input(&w, 0, cigar, seq);
        let reference = reference(&w);
        w.command("tag", &input, &["-g", reference.to_str().unwrap()])
            .assert()
            .success();
        let mut reader = bam::Reader::from_path(w.path("output.tsv")).unwrap();
        let r = reader.records().next().unwrap().unwrap();
        assert_eq!(r.aux(b"XM").unwrap(), bam::record::Aux::String(expected));
    }
}

#[test]
fn bam_output_is_compressed_bam_with_the_expected_calls() {
    let w = Workspace::new();
    let input = input(&w, 0, "4M", "CGCG");
    let reference = reference(&w);
    let output = w.path("tagged.bam");
    assert_cmd::Command::cargo_bin("metheor")
        .unwrap()
        .arg("tag")
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(&output)
        .arg("-g")
        .arg(reference)
        .assert()
        .success();
    assert_eq!(&fs::read(&output).unwrap()[..2], &[0x1f, 0x8b]);
    let mut reader = bam::Reader::from_path(&output).unwrap();
    let r = reader.records().next().unwrap().unwrap();
    assert_eq!(r.aux(b"XM").unwrap(), bam::record::Aux::String("Z.Z."));
}

#[test]
fn xm_strand_and_clipped_padding_follow_single_and_paired_end_orientation() {
    for (flags, expected) in [
        (0, "..Z.Z..."),
        (16, "...Z.Z.."),
        (99, "..Z.Z..."),
        (147, "..Z.Z..."),
        (83, "...Z.Z.."),
        (163, "...Z.Z.."),
    ] {
        let w = Workspace::new();
        let input = input(&w, flags, "2S4M2S", "AACGCGAA");
        let reference = reference(&w);
        w.command("tag", &input, &["-g", reference.to_str().unwrap()])
            .assert()
            .success();
        let mut reader = bam::Reader::from_path(w.path("output.tsv")).unwrap();
        let r = reader.records().next().unwrap().unwrap();
        assert_eq!(
            r.aux(b"XM").unwrap(),
            bam::record::Aux::String(expected),
            "flags={flags}"
        );
    }
}

#[test]
fn xm_pads_missing_flanks_and_keeps_records_across_chromosomes() {
    let cases = [
        ("right_one_flank", 0, "chr1", 7, "C", "U"),
        ("right_no_flanks", 0, "chr1", 8, "T", "u"),
        ("left_no_flanks", 16, "chr17", 1, "G", "U"),
        ("left_one_flank", 16, "chr17", 2, "A", "u"),
        ("interior", 0, "chr17", 4, "C", "Z"),
        ("final", 0, "chr17", 8, "C", "U"),
    ];
    for format in ["sam", "bam"] {
        let w = Workspace::new();
        let reference = w.write("reference.fa", ">chr1\nGGAAAACC\n>chr17\nGGACGTCC\n");
        rust_htslib::faidx::build(&reference).unwrap();
        let mut sam =
            String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:8\n@SQ\tSN:chr17\tLN:8\n");
        for (name, flags, chrom, pos, seq, _) in cases {
            sam.push_str(&format!(
                "{name}\t{flags}\t{chrom}\t{pos}\t60\t1M\t*\t0\t0\t{seq}\t*\n"
            ));
        }
        let mut input = w.write("input.sam", &sam);
        if format == "bam" {
            let mut reader = bam::Reader::from_path(&input).unwrap();
            input = w.path("input.bam");
            let header = bam::Header::from_template(reader.header());
            let mut writer = bam::Writer::from_path(&input, &header, bam::Format::Bam).unwrap();
            for record in reader.records() {
                writer.write(&record.unwrap()).unwrap();
            }
        }
        w.command("tag", &input, &["-g", reference.to_str().unwrap()])
            .assert()
            .success();
        let mut reader = bam::Reader::from_path(w.path("output.tsv")).unwrap();
        let header = reader.header().to_owned();
        let records: Vec<_> = reader.records().map(Result::unwrap).collect();
        assert_eq!(records.len(), cases.len(), "{format}");
        for (record, (name, flags, chrom, pos, seq, expected)) in records.iter().zip(cases) {
            assert_eq!(record.qname(), name.as_bytes(), "{format}");
            assert_eq!(record.flags(), flags);
            assert_eq!(header.tid2name(record.tid() as u32), chrom.as_bytes());
            assert_eq!(record.pos(), pos - 1);
            assert_eq!(record.seq().as_bytes(), seq.as_bytes());
            assert_eq!(
                record.aux(b"XM").unwrap(),
                bam::record::Aux::String(expected)
            );
        }
    }
}
