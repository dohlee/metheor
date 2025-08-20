use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs

#[test]
fn simple_run() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("metheor")?;

    cmd.arg("me")
        .arg("-i")
        .arg("tests/test1.bam")
        .arg("-o")
        .arg("tests/test1.me.tsv")
        .assert()
        .success();

    Ok(())
}
#[test]
fn input_bam_doesnt_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("metheor")?;

    cmd.arg("me")
        .arg("-i")
        .arg("tests/no_such.bam")
        .arg("-o")
        .arg("tests/test1.me.tsv")
        .assert()
        .failure()
        .stderr(predicate::str::contains("file not found"))
        .stderr(predicate::str::contains("no_such.bam"));

    Ok(())
}
