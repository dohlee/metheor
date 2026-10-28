mod common;

use assert_cmd::Command;
use common::Workspace;
use rust_htslib::{bam, bam::Read};

#[test]
fn tag_stdout_is_compressed_bam_without_progress_text() {
    let w = Workspace::new();
    let reference = w.write("reference.fa", ">chr1\nCCGCGTCC\n");
    rust_htslib::faidx::build(&reference).unwrap();
    let input = w.write(
        "input.sam",
        "@SQ\tSN:chr1\tLN:8\nread\t0\tchr1\t1\t60\t8M\t*\t0\t0\tCCGCGTCC\t*\n",
    );
    let stdout = Command::cargo_bin("metheor")
        .unwrap()
        .args([
            "tag",
            "-i",
            input.to_str().unwrap(),
            "-g",
            reference.to_str().unwrap(),
            "-o",
            "-",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    assert!(stdout.starts_with(&[0x1f, 0x8b]));
    let output = w.path("stream.bam");
    std::fs::write(&output, stdout).unwrap();
    let mut reader = bam::Reader::from_path(output).unwrap();
    assert_eq!(reader.header().target_names(), [b"chr1".as_slice()]);
    let records: Vec<_> = reader.records().map(Result::unwrap).collect();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].qname(), b"read");
    assert_eq!(records[0].seq().as_bytes(), b"CCGCGTCC");
    assert_eq!(
        records[0].aux(b"XM").unwrap(),
        bam::record::Aux::String("XZ.Z..UU")
    );
}

#[test]
fn tag_accepts_bare_relative_bam_and_sam_output_filenames() {
    for extension in ["bam", "sam"] {
        let w = Workspace::new();
        let reference = w.write("reference.fa", ">chr1\nCCGCGTCC\n");
        rust_htslib::faidx::build(&reference).unwrap();
        let input = w.write(
            "input.sam",
            "@SQ\tSN:chr1\tLN:8\nread\t0\tchr1\t1\t60\t8M\t*\t0\t0\tCCGCGTCC\t*\n",
        );
        let name = format!("tagged.{extension}");
        Command::cargo_bin("metheor")
            .unwrap()
            .current_dir(w.path(""))
            .args([
                "tag",
                "-i",
                input.to_str().unwrap(),
                "-g",
                reference.to_str().unwrap(),
                "-o",
                &name,
            ])
            .assert()
            .success();
        let bytes = std::fs::read(w.path(&name)).unwrap();
        assert_eq!(bytes.starts_with(&[0x1f, 0x8b]), extension == "bam");
        let mut reader = bam::Reader::from_path(w.path(&name)).unwrap();
        let records: Vec<_> = reader.records().map(Result::unwrap).collect();
        assert_eq!(records.len(), 1);
        assert_eq!(
            records[0].aux(b"XM").unwrap(),
            bam::record::Aux::String("XZ.Z..UU")
        );
    }
}

#[test]
fn tag_rejects_reusing_the_input_file_without_changing_it() {
    let w = Workspace::new();
    let reference = w.write("reference.fa", ">chr1\nCCGCGTCC\n");
    rust_htslib::faidx::build(&reference).unwrap();
    let input = w.write(
        "input.sam",
        "@SQ\tSN:chr1\tLN:8\nread\t0\tchr1\t1\t60\t8M\t*\t0\t0\tCCGCGTCC\t*\n",
    );
    let before = std::fs::read(&input).unwrap();
    Command::cargo_bin("metheor")
        .unwrap()
        .timeout(std::time::Duration::from_secs(5))
        .args([
            "tag",
            "-i",
            input.to_str().unwrap(),
            "-g",
            reference.to_str().unwrap(),
            "-o",
            input.to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains("different"));
    assert_eq!(std::fs::read(&input).unwrap(), before);
}

#[test]
#[cfg(unix)]
fn tag_rejects_source_aliases_and_reference_overwrites() {
    let w = Workspace::new();
    let reference = w.write("reference.fa", ">chr1\nCCGCGTCC\n");
    rust_htslib::faidx::build(&reference).unwrap();
    let input = w.write(
        "input.sam",
        "@SQ\tSN:chr1\tLN:8\nread\t0\tchr1\t1\t60\t8M\t*\t0\t0\tCCGCGTCC\t*\n",
    );
    let hardlink = w.path("hardlink.bam");
    std::fs::hard_link(&input, &hardlink).unwrap();
    let symlink = w.path("symlink.bam");
    std::os::unix::fs::symlink(&input, &symlink).unwrap();
    let input_before = std::fs::read(&input).unwrap();
    let reference_before = std::fs::read(&reference).unwrap();
    for output in [&hardlink, &symlink, &reference] {
        Command::cargo_bin("metheor")
            .unwrap()
            .args([
                "tag",
                "-i",
                input.to_str().unwrap(),
                "-g",
                reference.to_str().unwrap(),
                "-o",
                output.to_str().unwrap(),
            ])
            .assert()
            .failure()
            .stderr(predicates::str::contains("different"));
        assert_eq!(std::fs::read(&input).unwrap(), input_before);
        assert_eq!(std::fs::read(&reference).unwrap(), reference_before);
    }
}

#[test]
#[cfg(unix)]
fn tag_stdout_reports_a_failure_when_the_downstream_is_already_closed() {
    use std::{
        os::{fd::OwnedFd, unix::net::UnixStream},
        process::{Command, Stdio},
    };
    let w = Workspace::new();
    let reference = w.write("reference.fa", ">chr1\nCCGCGTCC\n");
    rust_htslib::faidx::build(&reference).unwrap();
    let input = w.write(
        "input.sam",
        "@SQ\tSN:chr1\tLN:8\nread\t0\tchr1\t1\t60\t8M\t*\t0\t0\tCCGCGTCC\t*\n",
    );
    let (downstream, upstream) = UnixStream::pair().unwrap();
    drop(downstream);
    let descriptor: OwnedFd = upstream.into();
    let output = Command::new(env!("CARGO_BIN_EXE_metheor"))
        .args([
            "tag",
            "-i",
            input.to_str().unwrap(),
            "-g",
            reference.to_str().unwrap(),
            "-o",
            "-",
        ])
        .stdout(Stdio::from(descriptor))
        .stderr(Stdio::piped())
        .output()
        .unwrap();
    assert!(
        !output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}
