mod common;

use assert_cmd::Command;
use common::Workspace;
use predicates::prelude::*;

#[cfg(test)]
mod cli_error_tests {
    use super::*;

    // Test core CLI infrastructure
    #[test]
    fn test_help_output() -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("--help");

        cmd.assert()
            .success()
            .stdout(predicate::str::contains("Usage:"))
            .stdout(predicate::str::contains("Commands:"))
            .stdout(predicate::str::contains("pdr"))
            .stdout(predicate::str::contains("fdrp"))
            .stdout(predicate::str::contains("tag"));

        Ok(())
    }

    #[test]
    fn test_version_output() -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("--version");

        cmd.assert()
            .success()
            .stdout(predicate::str::contains("metheor"));

        Ok(())
    }

    #[test]
    fn test_no_arguments_shows_help() -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin("metheor")?;

        cmd.assert()
            .failure()
            .stderr(predicate::str::contains("Usage:"));

        Ok(())
    }

    #[test]
    fn test_invalid_subcommand() -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("invalid_command");

        cmd.assert()
            .failure()
            .stderr(predicate::str::contains("error:"))
            .stderr(predicate::str::contains("subcommand"));

        Ok(())
    }

    // Test argument validation for each subcommand
    fn test_subcommand_missing_input(subcommand: &str) -> Result<(), Box<dyn std::error::Error>> {
        let w = Workspace::new();
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg(subcommand)
            .arg("--output")
            .arg(w.path("output.tsv"));

        cmd.assert()
            .failure()
            .stderr(predicate::str::contains("required"));

        Ok(())
    }

    fn test_subcommand_missing_output(subcommand: &str) -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg(subcommand).arg("--input").arg("tests/test1.bam");

        cmd.assert()
            .failure()
            .stderr(predicate::str::contains("required"));

        Ok(())
    }

    #[test]
    fn test_pdr_missing_required_args() -> Result<(), Box<dyn std::error::Error>> {
        test_subcommand_missing_input("pdr")?;
        test_subcommand_missing_output("pdr")?;
        Ok(())
    }

    #[test]
    fn test_lpmd_missing_required_args() -> Result<(), Box<dyn std::error::Error>> {
        test_subcommand_missing_input("lpmd")?;
        test_subcommand_missing_output("lpmd")?;
        Ok(())
    }

    #[test]
    fn test_mhl_missing_required_args() -> Result<(), Box<dyn std::error::Error>> {
        test_subcommand_missing_input("mhl")?;
        test_subcommand_missing_output("mhl")?;
        Ok(())
    }

    #[test]
    fn test_pm_missing_required_args() -> Result<(), Box<dyn std::error::Error>> {
        test_subcommand_missing_input("pm")?;
        test_subcommand_missing_output("pm")?;
        Ok(())
    }

    #[test]
    fn test_me_missing_required_args() -> Result<(), Box<dyn std::error::Error>> {
        test_subcommand_missing_input("me")?;
        test_subcommand_missing_output("me")?;
        Ok(())
    }

    #[test]
    fn test_fdrp_missing_required_args() -> Result<(), Box<dyn std::error::Error>> {
        test_subcommand_missing_input("fdrp")?;
        test_subcommand_missing_output("fdrp")?;
        Ok(())
    }

    #[test]
    fn test_qfdrp_missing_required_args() -> Result<(), Box<dyn std::error::Error>> {
        test_subcommand_missing_input("qfdrp")?;
        test_subcommand_missing_output("qfdrp")?;
        Ok(())
    }

    // Test parameter validation
    #[test]
    fn test_invalid_min_depth_negative() -> Result<(), Box<dyn std::error::Error>> {
        let w = Workspace::new();
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("pdr")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(w.path("output.tsv"))
            .arg("--min-depth")
            .arg("-5");

        cmd.assert()
            .failure()
            .stderr(predicate::str::contains("invalid").or(predicate::str::contains("error")));

        Ok(())
    }

    #[test]
    fn test_invalid_quality_threshold_out_of_range() -> Result<(), Box<dyn std::error::Error>> {
        let w = Workspace::new();
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("pdr")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(w.path("output.tsv"))
            .arg("--min-qual")
            .arg("300"); // Quality scores are typically 0-60

        // This might not fail at argument parsing but at runtime
        // depending on implementation
        cmd.assert().failure();

        Ok(())
    }

    #[test]
    fn test_lpmd_invalid_distance_parameters() -> Result<(), Box<dyn std::error::Error>> {
        let w = Workspace::new();
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("lpmd")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(w.path("output.tsv"))
            .arg("--min-distance")
            .arg("100")
            .arg("--max-distance")
            .arg("50"); // min > max

        cmd.assert()
            .failure()
            .stderr(predicate::str::contains("min-distance <= max-distance"));

        Ok(())
    }

    // Test file system error conditions
    #[test]
    fn test_output_directory_readonly() -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("pdr")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg("/nonexistent_directory/readonly_output.tsv");

        // This should fail, but the specific error might vary by system
        cmd.assert().failure();

        Ok(())
    }

    #[test]
    fn test_invalid_bam_file_format() -> Result<(), Box<dyn std::error::Error>> {
        let w = Workspace::new();
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("pdr")
            .arg("--input")
            .arg("Cargo.toml") // Not a BAM file
            .arg("--output")
            .arg(w.path("output.tsv"));

        cmd.assert()
            .failure()
            .stderr(predicate::str::contains("Error opening BAM file"));

        Ok(())
    }

    #[test]
    fn test_nonexistent_cpg_set_file() -> Result<(), Box<dyn std::error::Error>> {
        let w = Workspace::new();
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("pdr")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(w.path("output.tsv"))
            .arg("--cpg-set")
            .arg("nonexistent.bed");

        cmd.assert().failure();

        Ok(())
    }

    // Test tag-specific error conditions
    #[test]
    fn test_tag_missing_reference() -> Result<(), Box<dyn std::error::Error>> {
        let w = Workspace::new();
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("tag")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(w.path("output.sam"));
        // Missing required --genome argument

        cmd.assert()
            .failure()
            .stderr(predicate::str::contains("required"));

        Ok(())
    }

    #[test]
    fn test_tag_nonexistent_reference() -> Result<(), Box<dyn std::error::Error>> {
        let w = Workspace::new();
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("tag")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(w.path("output.sam"))
            .arg("--genome")
            .arg("nonexistent.fa");

        cmd.assert().failure();

        Ok(())
    }

    // Test subcommand help messages
    #[test]
    fn test_subcommand_help_pdr() -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("pdr").arg("--help");

        cmd.assert()
            .success()
            .stdout(predicate::str::contains("PDR"))
            .stdout(predicate::str::contains("--input"))
            .stdout(predicate::str::contains("--output"));

        Ok(())
    }

    #[test]
    fn test_subcommand_help_lpmd() -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("lpmd").arg("--help");

        cmd.assert()
            .success()
            .stdout(predicate::str::contains("LPMD"))
            .stdout(predicate::str::contains("--min-distance"))
            .stdout(predicate::str::contains("--max-distance"));

        Ok(())
    }

    #[test]
    fn test_subcommand_help_tag() -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("tag").arg("--help");

        cmd.assert()
            .success()
            .stdout(predicate::str::contains("Add bismark XM tag"))
            .stdout(predicate::str::contains("--genome")); // It's --genome, not --reference

        Ok(())
    }

    // Test edge cases with empty or minimal files
    #[test]
    fn test_zero_byte_input() -> Result<(), Box<dyn std::error::Error>> {
        let w = Workspace::new();
        // Create a temporary empty file
        let empty_file = w.write("empty.bam", "");

        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("pdr")
            .arg("--input")
            .arg(empty_file)
            .arg("--output")
            .arg(w.path("output.tsv"));

        let _result = cmd.assert().failure();

        Ok(())
    }

    // Test parameter boundary conditions
    #[test]
    fn test_extreme_parameter_values() -> Result<(), Box<dyn std::error::Error>> {
        let w = Workspace::new();
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("pdr")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(w.path("output.tsv"))
            .arg("--min-depth")
            .arg("999999");

        // Should run but likely produce no results due to high threshold
        cmd.assert().success();

        Ok(())
    }

    #[test]
    fn test_zero_min_qual_threshold() -> Result<(), Box<dyn std::error::Error>> {
        let w = Workspace::new();
        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("pdr")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(w.path("output.tsv"))
            .arg("--min-qual")
            .arg("0");

        cmd.assert().success();

        Ok(())
    }
}
