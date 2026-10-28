mod common;

use common::{assert_close, record, Workspace};

// S1 Text, Algorithms 1–11. Literals are hand-worked expectations, not snapshots
// of Metheor output. Four CpGs at 100, 102, 104, 106; all read pairs overlap.
fn check_metric(metric: &str, expected: [f64; 4]) {
    let cases = [
        vec!["z.z.z.z.", "z.z.z.z."],
        vec!["Z.Z.Z.Z.", "Z.Z.Z.Z."],
        vec!["z.z.z.z.", "Z.Z.Z.Z."],
        vec![
            "z.z.z.z.", "z.z.z.Z.", "z.z.Z.z.", "z.z.Z.Z.", "z.Z.z.z.", "z.Z.z.Z.", "z.Z.Z.z.",
            "z.Z.Z.Z.", "Z.z.z.z.", "Z.z.z.Z.", "Z.z.Z.z.", "Z.z.Z.Z.", "Z.Z.z.z.", "Z.Z.z.Z.",
            "Z.Z.Z.z.", "Z.Z.Z.Z.",
        ],
    ];
    for (patterns, expected) in cases.iter().zip(expected) {
        let w = Workspace::new();
        let records: String = patterns
            .iter()
            .enumerate()
            .map(|(i, xm)| record(&format!("r{i}"), 0, 101, 60, "8M", xm))
            .collect();
        let input = w.sam(&records);
        let mut args = vec!["--min-qual", "0"];
        if metric != "lpmd" {
            args.extend(["--min-depth", "2"]);
        }
        if matches!(metric, "fdrp" | "qfdrp") {
            args.extend(["--min-overlap", "1", "--max-depth", "40"]);
        }
        w.command(metric, &input, &args).assert().success();
        let rows = w.rows();
        let value_column = match metric {
            "pdr" | "mhl" | "fdrp" | "qfdrp" => 3,
            "lpmd" => 1,
            _ => 5,
        };
        let data = if metric == "lpmd" {
            &rows[1..]
        } else {
            &rows[..]
        };
        assert_eq!(
            data.len(),
            if matches!(metric, "lpmd" | "pm" | "me") {
                1
            } else {
                4
            }
        );
        for (i, row) in data.iter().enumerate() {
            assert_eq!(
                row.len(),
                match metric {
                    "pdr" | "pm" | "me" => 6,
                    "lpmd" => 2,
                    _ => 4,
                }
            );
            match metric {
                "lpmd" => {
                    assert_eq!(rows[0], ["name", "lpmd"]);
                    assert_eq!(row[0], input.to_str().unwrap());
                }
                "pm" | "me" => assert_eq!(&row[..5], &["chr1", "100", "102", "104", "106"]),
                _ => {
                    assert_eq!(row[0], "chr1");
                    assert_eq!(row[1].parse::<u32>().unwrap(), [100, 102, 104, 106][i]);
                    assert_eq!(row[2].parse::<u32>().unwrap(), [102, 104, 106, 108][i]);
                }
            }
            assert_close(row[value_column].parse().unwrap(), expected);
            if metric == "pdr" {
                let counts = if patterns.len() == 16 {
                    (2, 14)
                } else {
                    (2, 0)
                };
                assert_eq!(
                    (
                        row[4].parse::<u32>().unwrap(),
                        row[5].parse::<u32>().unwrap()
                    ),
                    counts
                );
            }
        }
    }
}

#[test]
fn pdr_matches_hand_worked_patterns() {
    check_metric("pdr", [0.0, 0.0, 0.0, 0.875]);
}

#[test]
fn lpmd_matches_hand_worked_patterns() {
    check_metric("lpmd", [0.0, 0.0, 0.0, 0.5]);
}

#[test]
fn mhl_matches_hand_worked_patterns() {
    check_metric("mhl", [0.0, 1.0, 0.5, 0.1625]);
}

#[test]
fn pm_matches_hand_worked_patterns() {
    check_metric("pm", [0.0, 0.0, 0.5, 0.9375]);
}

#[test]
fn me_matches_hand_worked_patterns() {
    check_metric("me", [0.0, 0.0, 0.25, 1.0]);
}

#[test]
fn fdrp_matches_hand_worked_patterns() {
    check_metric("fdrp", [0.0, 0.0, 1.0, 1.0]);
}

#[test]
fn qfdrp_matches_hand_worked_patterns() {
    check_metric("qfdrp", [0.0, 0.0, 1.0, 0.5333333333333333]);
}
