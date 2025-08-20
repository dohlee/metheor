use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::fs;
use std::process::Command; // Run programs

#[test]
fn input_file_doesnt_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("metheor")?;

    cmd.arg("tag")
        .arg("-i")
        .arg("tests/no_such_input.bam")
        .arg("-o")
        .arg("tests/out.bam")
        .arg("-g")
        .arg("tests/hg38.chr19.fa")
        .assert()
        .failure()
        .stderr(predicate::str::contains("file not found"))
        .stderr(predicate::str::contains("no_such_input.bam"));

    Ok(())
}
#[test]
fn output_directory_doesnt_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("metheor")?;

    cmd.arg("tag")
        .arg("-i")
        .arg("tests/test.chr19.noXM.sam")
        .arg("-o")
        .arg("no_such_dir/out.bam")
        .arg("-g")
        .arg("tests/hg38.chr19.fa")
        .assert()
        .failure()
        .stderr(predicate::str::contains("No such directory"))
        .stderr(predicate::str::contains("no_such_dir"));

    Ok(())
}
#[test]
fn reference_genome_doesnt_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("metheor")?;

    cmd.arg("tag")
        .arg("-i")
        .arg("tests/test.chr19.noXM.sam")
        .arg("-o")
        .arg("tests/out.bam")
        .arg("-g")
        .arg("tests/no_such.fa")
        .assert()
        .failure()
        .stderr(predicate::str::contains("file not found"))
        .stderr(predicate::str::contains("no_such.fa"));

    Ok(())
}
#[test]
fn test_whether_xmtag_generated_correctly_for_1000_chr19_reads(
) -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("metheor")?;

    cmd.arg("tag")
        .arg("-i")
        .arg("tests/test.chr19.noXM.sam")
        .arg("-o")
        .arg("tests/test.chr19.metheor_tag_out.sam")
        .arg("-g")
        .arg("tests/hg38.chr19.fa")
        .assert()
        .success();

    let original = fs::read_to_string("tests/test.chr19.XM.sam")?;
    let generated = fs::read_to_string("tests/test.chr19.metheor_tag_out.sam")?;
    assert_eq!(original, generated);

    Ok(())
}
