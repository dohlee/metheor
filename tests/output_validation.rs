use assert_cmd::Command;
use std::fs;
use std::path::Path;

#[cfg(test)]
mod output_validation_tests {
    use super::*;

    fn clean_output_file(path: &str) {
        if Path::new(path).exists() {
            let _ = fs::remove_file(path);
        }
    }

    // Helper function to validate output file structure and basic content
    fn validate_output_structure(
        output_file: &str,
        expected_columns: usize,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let content = fs::read_to_string(output_file)?;
        let lines: Vec<&str> = content.lines().collect();

        // Should have at least one line of content
        assert!(!lines.is_empty(), "Output file should not be empty");

        // Check that each line has consistent number of fields
        for (i, line) in lines.iter().enumerate() {
            if line.trim().is_empty() {
                continue; // Skip empty lines
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if expected_columns > 0 {
                assert_eq!(
                    fields.len(),
                    expected_columns,
                    "Line {} has {} fields, expected {}",
                    i + 1,
                    fields.len(),
                    expected_columns
                );
            }

            // Ensure no completely empty fields (except last field which might be empty)
            for (j, field) in fields.iter().enumerate() {
                if j < fields.len() - 1 {
                    assert!(
                        !field.trim().is_empty(),
                        "Empty field at line {}, column {}",
                        i + 1,
                        j + 1
                    );
                }
            }
        }

        Ok(())
    }

    #[test]
    fn test_pdr_output_correctness() -> Result<(), Box<dyn std::error::Error>> {
        let output_file = "test_pdr_output.tsv";
        clean_output_file(output_file);

        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("pdr")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(output_file)
            .arg("--min-depth")
            .arg("1")
            .arg("--min-cpgs")
            .arg("1")
            .arg("--min-qual")
            .arg("10");

        cmd.assert().success();

        // Verify output file exists and has content
        assert!(Path::new(output_file).exists());
        let content = fs::read_to_string(output_file)?;
        assert!(!content.is_empty());

        // Validate output structure - PDR should have 6 columns (chr, start, end, pdr_value, concordant, discordant)
        validate_output_structure(output_file, 6)?;

        // Validate that we have numeric PDR values in the last column
        let lines: Vec<&str> = content.lines().collect();
        for line in lines.iter() {
            // Skip header if present
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 6 {
                // Try to parse the PDR value (4th column, 0-indexed as 3)
                let pdr_str = fields[3];
                if let Ok(pdr_val) = pdr_str.parse::<f64>() {
                    assert!(
                        (0.0..=1.0).contains(&pdr_val),
                        "PDR value {} out of bounds [0,1]",
                        pdr_val
                    );
                }
            }
        }

        clean_output_file(output_file);
        Ok(())
    }

    #[test]
    fn test_lpmd_output_correctness() -> Result<(), Box<dyn std::error::Error>> {
        let output_file = "test_lpmd_output.tsv";
        clean_output_file(output_file);

        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("lpmd")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(output_file)
            .arg("--min-distance")
            .arg("1")
            .arg("--max-distance")
            .arg("1000")
            .arg("--min-qual")
            .arg("10");

        cmd.assert().success();

        // Verify output
        assert!(Path::new(output_file).exists());
        let content = fs::read_to_string(output_file)?;
        assert!(!content.is_empty());

        // LPMD outputs summary format with name and lpmd columns
        validate_output_structure(output_file, 2)?;

        // Check that the output contains expected structure
        let lines: Vec<&str> = content.lines().collect();
        assert!(lines.len() >= 2, "LPMD output should have header and data");

        // First line should be header
        let header = lines[0];
        assert!(
            header.contains("lpmd") || header.contains("name"),
            "Header should contain lpmd or name: {}",
            header
        );

        clean_output_file(output_file);
        Ok(())
    }

    #[test]
    fn test_mhl_output_correctness() -> Result<(), Box<dyn std::error::Error>> {
        let output_file = "test_mhl_output.tsv";
        clean_output_file(output_file);

        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("mhl")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(output_file)
            .arg("--min-depth")
            .arg("1")
            .arg("--min-cpgs")
            .arg("1")
            .arg("--min-qual")
            .arg("10");

        cmd.assert().success();

        assert!(Path::new(output_file).exists());
        let content = fs::read_to_string(output_file)?;
        assert!(!content.is_empty());

        // MHL typically outputs 4 columns (chr, start, end, mhl_value)
        validate_output_structure(output_file, 4)?;

        // Validate MHL values are in valid range
        let lines: Vec<&str> = content.lines().collect();
        for line in lines.iter() {
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 4 {
                let mhl_str = fields[3];
                if let Ok(mhl_val) = mhl_str.parse::<f64>() {
                    assert!(
                        (0.0..=1.0).contains(&mhl_val),
                        "MHL value {} out of bounds [0,1]",
                        mhl_val
                    );
                }
            }
        }

        clean_output_file(output_file);
        Ok(())
    }

    #[test]
    fn test_pm_output_correctness() -> Result<(), Box<dyn std::error::Error>> {
        let output_file = "test_pm_output.tsv";
        clean_output_file(output_file);

        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("pm")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(output_file)
            .arg("--min-depth")
            .arg("1")
            .arg("--min-qual")
            .arg("10");

        cmd.assert().success();

        assert!(Path::new(output_file).exists());
        let content = fs::read_to_string(output_file)?;
        assert!(!content.is_empty());

        // PM typically outputs 6 columns (chr, start, end, pos1, pos2, pm_value)
        validate_output_structure(output_file, 6)?;

        // Validate PM values are in valid range (last column)
        let lines: Vec<&str> = content.lines().collect();
        for line in lines.iter() {
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 6 {
                let pm_str = fields[5]; // PM value is in the last column
                if let Ok(pm_val) = pm_str.parse::<f64>() {
                    assert!(
                        (0.0..=1.0).contains(&pm_val),
                        "PM value {} out of bounds [0,1]",
                        pm_val
                    );
                }
            }
        }

        clean_output_file(output_file);
        Ok(())
    }

    #[test]
    fn test_me_output_correctness() -> Result<(), Box<dyn std::error::Error>> {
        let output_file = "test_me_output.tsv";
        clean_output_file(output_file);

        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("me")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(output_file)
            .arg("--min-depth")
            .arg("1")
            .arg("--min-qual")
            .arg("10");

        cmd.assert().success();

        assert!(Path::new(output_file).exists());
        let content = fs::read_to_string(output_file)?;
        assert!(!content.is_empty());

        // ME typically outputs 6 columns (chr, start, end, pos1, pos2, me_value)
        validate_output_structure(output_file, 6)?;

        // Validate ME values are in valid range (last column)
        let lines: Vec<&str> = content.lines().collect();
        for line in lines.iter() {
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 6 {
                let me_str = fields[5]; // ME value is in the last column
                if let Ok(me_val) = me_str.parse::<f64>() {
                    assert!(
                        (0.0..=1.0).contains(&me_val),
                        "ME value {} out of bounds [0,1]",
                        me_val
                    );
                }
            }
        }

        clean_output_file(output_file);
        Ok(())
    }

    #[test]
    fn test_fdrp_output_correctness() -> Result<(), Box<dyn std::error::Error>> {
        let output_file = "test_fdrp_output.tsv";
        clean_output_file(output_file);

        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("fdrp")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(output_file)
            .arg("--min-depth")
            .arg("1")
            .arg("--max-depth")
            .arg("100")
            .arg("--min-overlap")
            .arg("1")
            .arg("--min-qual")
            .arg("10");

        cmd.assert().success();

        assert!(Path::new(output_file).exists());
        let content = fs::read_to_string(output_file)?;
        assert!(!content.is_empty());

        // FDRP typically outputs 4 columns (chr, start, end, fdrp_value)
        validate_output_structure(output_file, 4)?;

        // Validate FDRP values are in valid range
        let lines: Vec<&str> = content.lines().collect();
        for line in lines.iter() {
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 4 {
                let fdrp_str = fields[3];
                if let Ok(fdrp_val) = fdrp_str.parse::<f64>() {
                    assert!(
                        (0.0..=1.0).contains(&fdrp_val),
                        "FDRP value {} out of bounds [0,1]",
                        fdrp_val
                    );
                }
            }
        }

        clean_output_file(output_file);
        Ok(())
    }

    #[test]
    fn test_qfdrp_output_correctness() -> Result<(), Box<dyn std::error::Error>> {
        let output_file = "test_qfdrp_output.tsv";
        clean_output_file(output_file);

        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("qfdrp")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(output_file)
            .arg("--min-depth")
            .arg("1")
            .arg("--max-depth")
            .arg("100")
            .arg("--min-overlap")
            .arg("1")
            .arg("--min-qual")
            .arg("10");

        cmd.assert().success();

        assert!(Path::new(output_file).exists());
        let content = fs::read_to_string(output_file)?;
        assert!(!content.is_empty());

        // qFDRP typically outputs 4 columns (chr, start, end, qfdrp_value)
        validate_output_structure(output_file, 4)?;

        // Validate qFDRP values are in valid range
        let lines: Vec<&str> = content.lines().collect();
        for line in lines.iter() {
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 4 {
                let qfdrp_str = fields[3];
                if let Ok(qfdrp_val) = qfdrp_str.parse::<f64>() {
                    assert!(
                        (0.0..=1.0).contains(&qfdrp_val),
                        "qFDRP value {} out of bounds [0,1]",
                        qfdrp_val
                    );
                }
            }
        }

        clean_output_file(output_file);
        Ok(())
    }

    // Test with multiple BAM files to ensure consistency
    #[test]
    fn test_pdr_multiple_bam_files() -> Result<(), Box<dyn std::error::Error>> {
        for i in 1..=6 {
            let input_file = format!("tests/test{}.bam", i);
            let output_file = format!("test_pdr_test{}.tsv", i);

            // Skip if input file doesn't exist
            if !Path::new(&input_file).exists() {
                continue;
            }

            clean_output_file(&output_file);

            let mut cmd = Command::cargo_bin("metheor")?;
            cmd.arg("pdr")
                .arg("--input")
                .arg(&input_file)
                .arg("--output")
                .arg(&output_file)
                .arg("--min-depth")
                .arg("1")
                .arg("--min-cpgs")
                .arg("1")
                .arg("--min-qual")
                .arg("10");

            cmd.assert().success();

            // Verify output was created
            assert!(
                Path::new(&output_file).exists(),
                "Output file not created for {}",
                input_file
            );

            let _content = fs::read_to_string(&output_file)?;
            // Content can be empty if no CpGs pass filters, but file should exist

            clean_output_file(&output_file);
        }

        Ok(())
    }

    // Test parameter variation effects
    #[test]
    fn test_min_depth_parameter_effects() -> Result<(), Box<dyn std::error::Error>> {
        let output_low = "test_depth_low.tsv";
        let output_high = "test_depth_high.tsv";

        clean_output_file(output_low);
        clean_output_file(output_high);

        // Run with low min-depth
        let mut cmd_low = Command::cargo_bin("metheor")?;
        cmd_low
            .arg("pdr")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(output_low)
            .arg("--min-depth")
            .arg("1");
        cmd_low.assert().success();

        // Run with high min-depth
        let mut cmd_high = Command::cargo_bin("metheor")?;
        cmd_high
            .arg("pdr")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(output_high)
            .arg("--min-depth")
            .arg("20");
        cmd_high.assert().success();

        // Compare output sizes - higher threshold should produce fewer results
        let content_low = fs::read_to_string(output_low)?;
        let content_high = fs::read_to_string(output_high)?;

        let lines_low = content_low.lines().count();
        let lines_high = content_high.lines().count();

        // Higher min-depth should result in same or fewer lines
        assert!(
            lines_high <= lines_low,
            "Higher min-depth should not produce more results"
        );

        clean_output_file(output_low);
        clean_output_file(output_high);
        Ok(())
    }

    // Test output file format validation
    #[test]
    fn test_output_file_format_structure() -> Result<(), Box<dyn std::error::Error>> {
        let output_file = "test_format.tsv";
        clean_output_file(output_file);

        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("pdr")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(output_file)
            .arg("--min-depth")
            .arg("1");

        cmd.assert().success();

        let content = fs::read_to_string(output_file)?;
        let lines: Vec<&str> = content.lines().collect();

        if !lines.is_empty() {
            // Check that each line has consistent number of tab-separated fields
            let first_line_fields = lines[0].split('\t').count();

            for (i, line) in lines.iter().enumerate() {
                let fields = line.split('\t').count();
                assert_eq!(
                    fields,
                    first_line_fields,
                    "Inconsistent field count at line {}",
                    i + 1
                );

                // Ensure no empty fields (except possibly last field)
                let field_vec: Vec<&str> = line.split('\t').collect();
                for (j, field) in field_vec.iter().enumerate() {
                    if j < field_vec.len() - 1 {
                        // Not last field
                        assert!(
                            !field.is_empty(),
                            "Empty field at line {}, position {}",
                            i + 1,
                            j + 1
                        );
                    }
                }
            }
        }

        clean_output_file(output_file);
        Ok(())
    }

    // Test output file permissions and creation
    #[test]
    fn test_output_file_permissions() -> Result<(), Box<dyn std::error::Error>> {
        let output_file = "test_permissions.tsv";
        clean_output_file(output_file);

        let mut cmd = Command::cargo_bin("metheor")?;
        cmd.arg("pdr")
            .arg("--input")
            .arg("tests/test1.bam")
            .arg("--output")
            .arg(output_file);

        cmd.assert().success();

        // Check file was created and is readable
        assert!(Path::new(output_file).exists());
        let _content = fs::read_to_string(output_file)?; // Should not fail

        // Check file metadata
        let metadata = fs::metadata(output_file)?;
        assert!(metadata.is_file());
        // File exists and is a file (length check is redundant since len() returns u64)

        clean_output_file(output_file);
        Ok(())
    }
}
