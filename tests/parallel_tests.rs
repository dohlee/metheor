use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;

#[test]
fn test_parallel_cli_options() {
    let mut cmd = Command::cargo_bin("metheor").unwrap();
    cmd.args(["--help"]);

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("--threads"))
        .stdout(predicate::str::contains("--parallel-threshold"));
}

#[test]
fn test_fdrp_with_custom_threads() {
    let mut cmd = Command::cargo_bin("metheor").unwrap();
    cmd.args([
        "--threads",
        "2",
        "--parallel-threshold",
        "50",
        "fdrp",
        "--input",
        "tests/test1.bam",
        "--output",
        "/tmp/test_parallel_fdrp.tsv",
        "--min-qual",
        "10",
        "--min-depth",
        "2",
        "--max-depth",
        "40",
        "--min-overlap",
        "35",
    ]);

    cmd.assert().success();
}

#[test]
fn test_qfdrp_with_custom_threads() {
    let mut cmd = Command::cargo_bin("metheor").unwrap();
    cmd.args([
        "--threads",
        "2",
        "--parallel-threshold",
        "50",
        "qfdrp",
        "--input",
        "tests/test1.bam",
        "--output",
        "/tmp/test_parallel_qfdrp.tsv",
        "--min-qual",
        "10",
        "--min-depth",
        "2",
        "--max-depth",
        "40",
        "--min-overlap",
        "35",
    ]);

    cmd.assert().success();
}

#[test]
fn test_parallel_threshold_validation() {
    // Test that parallel threshold of 0 works (forces parallel processing)
    let mut cmd = Command::cargo_bin("metheor").unwrap();
    cmd.args([
        "--parallel-threshold",
        "0",
        "fdrp",
        "--input",
        "tests/test1.bam",
        "--output",
        "/tmp/test_threshold_0.tsv",
        "--min-qual",
        "10",
        "--min-depth",
        "2",
        "--max-depth",
        "40",
        "--min-overlap",
        "35",
    ]);

    cmd.assert().success();
}

#[test]
fn test_auto_thread_detection() {
    // Test with threads=0 (auto-detect)
    let mut cmd = Command::cargo_bin("metheor").unwrap();
    cmd.args([
        "--threads",
        "0", // Auto-detect
        "fdrp",
        "--input",
        "tests/test1.bam",
        "--output",
        "/tmp/test_auto_threads.tsv",
        "--min-qual",
        "10",
        "--min-depth",
        "2",
        "--max-depth",
        "40",
        "--min-overlap",
        "35",
    ]);

    cmd.assert().success();
}
