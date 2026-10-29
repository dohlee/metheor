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
fn thread_counts_preserve_all_records_and_headers_across_sam_bam_and_batches() {
    let w = Workspace::new();
    let reference = reference(&w);
    let header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\n@RG\tID:sample\tSM:test\n@CO\tkeep this comment\n";
    let mut source = header.to_owned();
    let mut expected = header.to_owned();
    for i in 0..4103 {
        let (flags, cigar, seq, xm) = match i % 4 {
            0 => (0, "4M", "CGCG", "Z.Z."),
            1 => (0, "2S4M", "AACGCG", "..Z.Z."),
            2 => (0, "2M2I2M", "CGAACG", "Z...Z."),
            _ => (16, "2S4M2S", "AACGCGAA", "...Z.Z.."),
        };
        let record = format!("r{i}\t{flags}\tchr1\t1\t60\t{cigar}\t*\t0\t0\t{seq}\t*\tRG:Z:sample");
        source.push_str(&format!("{record}\n"));
        expected.push_str(&format!("{record}\tXM:Z:{xm}\n"));
    }
    let input_sam = w.write("input.sam", &source);
    let input_bam = w.path("input.bam");
    let mut reader = bam::Reader::from_path(&input_sam).unwrap();
    let mut writer = bam::Writer::from_path(
        &input_bam,
        &bam::Header::from_template(reader.header()),
        bam::Format::Bam,
    )
    .unwrap();
    for record in reader.records() {
        writer.write(&record.unwrap()).unwrap();
    }
    drop(writer);
    let mut expected_reader = bam::Reader::from_path(w.write("expected.sam", &expected)).unwrap();
    let expected_records: Vec<_> = expected_reader.records().map(Result::unwrap).collect();

    for input in [input_sam, input_bam] {
        for extension in ["sam", "bam"] {
            for threads in [None, Some("1"), Some("3")] {
                let output = w.path(&format!("output.{extension}"));
                let mut cmd = assert_cmd::Command::cargo_bin("metheor").unwrap();
                cmd.arg("tag")
                    .arg("-i")
                    .arg(&input)
                    .arg("-o")
                    .arg(&output)
                    .arg("-g")
                    .arg(&reference);
                if let Some(count) = threads {
                    cmd.args(["--threads", count]);
                }
                cmd.assert().success();
                let mut reader = bam::Reader::from_path(&output).unwrap();
                assert_eq!(reader.header().as_bytes(), header.as_bytes());
                let records: Vec<_> = reader.records().map(Result::unwrap).collect();
                assert_eq!(records, expected_records);
                if extension == "sam" {
                    assert_eq!(fs::read_to_string(&output).unwrap(), expected);
                }
            }
        }
    }
}

#[test]
fn invalid_thread_counts_are_rejected_before_touching_output() {
    let w = Workspace::new();
    let input = input(&w, 0, "4M", "CGCG");
    let reference = reference(&w);
    let output = w.write("output.tsv", "previous successful result\n");
    for count in ["0", "-1", "many", "18446744073709551616"] {
        w.command("tag", &input, &["-g", reference.to_str().unwrap()])
            .arg(format!("--threads={count}"))
            .assert()
            .code(2)
            .stderr(predicates::str::contains("--threads"));
        assert_eq!(
            fs::read_to_string(&output).unwrap(),
            "previous successful result\n"
        );
    }
}

#[test]
fn parallel_tagging_handles_empty_input_and_more_workers_than_records() {
    for count in [0, 1, 2] {
        let w = Workspace::new();
        let reference = reference(&w);
        let header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\n";
        let mut source = header.to_owned();
        let mut expected = header.to_owned();
        for i in 0..count {
            let record = format!("r{i}\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCGCG\t*");
            source.push_str(&format!("{record}\n"));
            expected.push_str(&format!("{record}\tXM:Z:Z.Z.\n"));
        }
        let input = w.write("input.sam", &source);
        w.command(
            "tag",
            &input,
            &["-g", reference.to_str().unwrap(), "--threads", "4"],
        )
        .assert()
        .success();
        assert_eq!(fs::read_to_string(w.path("output.tsv")).unwrap(), expected);
    }
}
