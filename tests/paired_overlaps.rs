mod common;

use common::{assert_close, record, Workspace};

fn paired_record(
    name: &str,
    flags: u16,
    pos: u32,
    mate_pos: u32,
    mapq: u8,
    cigar: &str,
    xm: &str,
) -> String {
    record(name, flags, pos, mapq, cigar, xm).replacen(
        "\t*\t0\t0\t",
        &format!("\t=\t{mate_pos}\t0\t"),
        1,
    )
}

fn overlapping_fragments() -> String {
    paired_record("fragment-a", 99, 101, 101, 60, "8M", "Z.Z.Z.Z.")
        + &paired_record("fragment-a", 147, 101, 101, 60, "8M", "Z.Z.Z.Z.")
        + &paired_record("fragment-b", 99, 101, 201, 60, "8M", "z.z.z.z.")
        + &paired_record("fragment-b", 147, 201, 101, 60, "8M", "........")
}

#[test]
fn pm_overlap_correction_counts_a_shared_quartet_once_per_fragment() {
    let w = Workspace::new();
    let input = w.sam(&overlapping_fragments());
    for (extra, expected) in [
        (vec![], 0.4444444444444444),
        (vec!["--ignore-overlaps"], 0.5),
    ] {
        let mut args = vec!["-d", "2"];
        args.extend(extra);
        w.command("pm", &input, &args).assert().success();
        let rows = w.rows();
        assert_eq!(rows.len(), 1);
        assert_close(rows[0][5].parse().unwrap(), expected);
    }
}

#[test]
fn me_overlap_correction_counts_a_shared_quartet_once_per_fragment() {
    let w = Workspace::new();
    let input = w.sam(&overlapping_fragments());
    w.command("me", &input, &["-d", "2", "--ignore-overlaps"])
        .assert()
        .success();
    let rows = w.rows();
    assert_eq!(rows.len(), 1);
    assert_close(rows[0][5].parse().unwrap(), 0.25);
}

#[test]
fn pdr_overlap_correction_uses_fragment_depth_at_each_cpg() {
    let w = Workspace::new();
    let input = w.sam(&overlapping_fragments());
    w.command("pdr", &input, &["-d", "2", "--ignore-overlaps"])
        .assert()
        .success();
    let rows = w.rows();
    assert_eq!(rows.len(), 4);
    for row in rows {
        assert_eq!(&row[3..], &["0", "2", "0"]);
    }
}

#[test]
fn mhl_overlap_correction_does_not_overweight_a_fully_methylated_fragment() {
    let w = Workspace::new();
    let input = w.sam(&overlapping_fragments());
    w.command("mhl", &input, &["-d", "2", "--ignore-overlaps"])
        .assert()
        .success();
    let rows = w.rows();
    assert_eq!(rows.len(), 4);
    for row in rows {
        assert_close(row[3].parse().unwrap(), 0.5);
    }
}

#[test]
fn lpmd_overlap_correction_counts_each_shared_pair_once_per_fragment() {
    let w = Workspace::new();
    let input = w.sam(&overlapping_fragments());
    let pairs = w.path("pairs.tsv");
    w.command(
        "lpmd",
        &input,
        &["--ignore-overlaps", "--pairs", pairs.to_str().unwrap()],
    )
    .assert()
    .success();
    let text = std::fs::read_to_string(pairs).unwrap();
    let rows: Vec<_> = text.lines().skip(1).collect();
    assert_eq!(rows.len(), 6);
    for row in rows {
        let fields: Vec<_> = row.split('\t').collect();
        assert_eq!(&fields[3..], &["0", "2", "0"]);
    }
}

#[test]
fn fdrp_overlap_correction_compares_distinct_fragments() {
    let w = Workspace::new();
    let input = w.sam(&overlapping_fragments());
    w.command("fdrp", &input, &["-d", "2", "-l", "1", "--ignore-overlaps"])
        .assert()
        .success();
    let rows = w.rows();
    assert_eq!(rows.len(), 4);
    for row in rows {
        assert_close(row[3].parse().unwrap(), 1.0);
    }
}

#[test]
fn qfdrp_overlap_correction_compares_distinct_fragments() {
    let w = Workspace::new();
    let input = w.sam(&overlapping_fragments());
    w.command(
        "qfdrp",
        &input,
        &["-d", "2", "-l", "1", "--ignore-overlaps"],
    )
    .assert()
    .success();
    let rows = w.rows();
    assert_eq!(rows.len(), 4);
    for row in rows {
        assert_close(row[3].parse().unwrap(), 1.0);
    }
}

#[test]
fn partial_overlap_keeps_unique_quartets_with_their_full_read_context() {
    let w = Workspace::new();
    let records = paired_record("pair", 99, 101, 103, 60, "10M", "Z.Z.Z.Z.Z.")
        + &record("other", 0, 101, 60, "12M", "z.z.z.z.z.z.")
        + &paired_record("pair", 147, 103, 101, 60, "10M", "z.z.z.z.z.");
    let input = w.sam(&records);
    for (metric, expected) in [("pm", 0.5), ("me", 0.25)] {
        w.command(metric, &input, &["-d", "2", "--ignore-overlaps"])
            .assert()
            .success();
        let mut rows = w.rows();
        rows.sort();
        assert_eq!(rows.len(), 3);
        assert_eq!(
            rows.iter().map(|r| r[1].as_str()).collect::<Vec<_>>(),
            ["100", "102", "104"]
        );
        for row in &rows[..2] {
            assert_close(row[5].parse().unwrap(), expected);
        }
        assert_close(rows[2][5].parse().unwrap(), 0.0);
    }
}

#[test]
fn pdr_keeps_full_read_concordance_for_nonoverlapping_cpgs() {
    let w = Workspace::new();
    let input = w.sam(
        &(paired_record("pair", 99, 101, 103, 60, "10M", "Z.Z.Z.Z.Z.")
            + &record("other", 0, 101, 60, "12M", "z.z.z.z.z.z.")
            + &paired_record("pair", 147, 103, 101, 60, "10M", "Z.Z.Z.z.Z.")),
    );
    w.command("pdr", &input, &["-d", "2", "--ignore-overlaps"])
        .assert()
        .success();
    let rows = w.rows();
    let unique = rows.iter().find(|r| r[1] == "110").unwrap();
    assert_eq!(&unique[3..], &["0.5", "1", "1"]);
    let shared = rows.iter().find(|r| r[1] == "108").unwrap();
    assert_eq!(&shared[3..], &["0", "2", "0"]);
}

#[test]
fn filtered_mates_do_not_suppress_a_qualifying_observation() {
    for metric in ["pdr", "mhl", "pm", "me", "fdrp", "qfdrp"] {
        let w = Workspace::new();
        let input = w.sam(
            &(paired_record("pair", 99, 101, 101, 9, "8M", "Z.Z.Z.Z.")
                + &paired_record("pair", 147, 101, 101, 60, "8M", "Z.Z.Z.Z.")
                + &record("other", 0, 101, 60, "8M", "z.z.z.z.")),
        );
        let mut args = vec!["-d", "2", "--ignore-overlaps"];
        if matches!(metric, "fdrp" | "qfdrp") {
            args.extend(["-l", "1"]);
        }
        w.command(metric, &input, &args).assert().success();
        let rows = w.rows();
        assert_eq!(
            rows.len(),
            if matches!(metric, "pm" | "me") { 1 } else { 4 }
        );
        let (column, expected) = match metric {
            "pm" => (5, 0.5),
            "me" => (5, 0.25),
            "pdr" => (3, 0.0),
            "mhl" => (3, 0.5),
            _ => (3, 1.0),
        };
        for row in rows {
            assert_close(row[column].parse().unwrap(), expected);
        }
    }
}

#[test]
fn read_groups_and_reciprocal_mate_coordinates_prevent_false_matches() {
    let first = paired_record("pair", 99, 101, 101, 60, "8M", "Z.Z.Z.Z.");
    let second = paired_record("pair", 147, 101, 101, 60, "8M", "Z.Z.Z.Z.");
    let other = record("other", 0, 101, 60, "8M", "z.z.z.z.");
    let different_groups =
        first.trim_end().to_owned() + "\tRG:Z:a\n" + second.trim_end() + "\tRG:Z:b\n" + &other;
    let mismatched_positions =
        paired_record("pair", 99, 101, 103, 60, "8M", "Z.Z.Z.Z.") + &second + &other;
    let same_read_end = first.clone() + &first + &other;
    for records in [different_groups, mismatched_positions, same_read_end] {
        let w = Workspace::new();
        let input = w.sam(&records);
        w.command("pdr", &input, &["-d", "1", "--ignore-overlaps"])
            .assert()
            .success();
        let rows = w.rows();
        assert_eq!(rows.len(), 4);
        for row in rows {
            assert_eq!(&row[4..], &["3", "0"]);
        }
    }
}

#[test]
fn coordinate_first_mate_is_retained_even_when_it_is_read_two() {
    let w = Workspace::new();
    let input = w.sam(
        &(paired_record("pair", 147, 101, 103, 60, "10M", "Z.Z.Z.Z.Z.")
            + &record("other", 0, 101, 60, "12M", "z.z.z.z.z.z.")
            + &paired_record("pair", 99, 103, 101, 60, "10M", "z.z.z.z.z.")),
    );
    w.command("pm", &input, &["-d", "2", "--ignore-overlaps"])
        .assert()
        .success();
    let rows = w.rows();
    assert_eq!(rows.len(), 3);
    let shared = rows.iter().find(|r| r[1] == "102").unwrap();
    assert_close(shared[5].parse().unwrap(), 0.5);
}

#[test]
fn shared_cpgs_do_not_imply_that_the_same_quartet_was_observed() {
    let w = Workspace::new();
    let input = w.sam(
        &(paired_record("pair", 99, 101, 101, 60, "10M", "Z.Z.Z.Z.Z.")
            + &paired_record("pair", 147, 101, 101, 60, "10M", "z...z.z.z.")
            + &record("other", 0, 101, 60, "10M", "z...z.z.z.")),
    );
    w.command("pm", &input, &["-d", "2", "--ignore-overlaps"])
        .assert()
        .success();
    let rows = w.rows();
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0], ["chr1", "100", "104", "106", "108", "0"]);
}

#[test]
fn rejected_long_read_observations_do_not_mask_a_usable_mate() {
    for metric in ["fdrp", "qfdrp"] {
        let w = Workspace::new();
        let long_xm = "Z.Z.Z.Z.".to_owned() + &".".repeat(292);
        let input = w.sam(
            &(paired_record("pair", 99, 101, 101, 60, "300M", &long_xm)
                + &paired_record("pair", 147, 101, 101, 60, "8M", "Z.Z.Z.Z.")
                + &record("other", 0, 101, 60, "8M", "z.z.z.z.")),
        );
        w.command(metric, &input, &["-d", "2", "-l", "1", "--ignore-overlaps"])
            .assert()
            .success();
        let rows = w.rows();
        assert_eq!(rows.len(), 4);
        for row in rows {
            assert_close(row[3].parse().unwrap(), 1.0);
        }
    }
}

#[test]
fn reads_below_the_cpg_threshold_do_not_mask_eligible_mates() {
    for metric in ["pdr", "mhl"] {
        let w = Workspace::new();
        let input = w.sam(
            &(paired_record("pair", 99, 101, 101, 60, "6M", "Z.Z.Z.")
                + &paired_record("pair", 147, 101, 101, 60, "8M", "Z.Z.Z.Z.")
                + &record("other", 0, 101, 60, "8M", "z.z.z.z.")),
        );
        w.command(metric, &input, &["-d", "2", "-p", "4", "--ignore-overlaps"])
            .assert()
            .success();
        let rows = w.rows();
        assert_eq!(rows.len(), 4);
        for row in rows {
            assert_close(
                row[3].parse().unwrap(),
                if metric == "mhl" { 0.5 } else { 0.0 },
            );
        }
    }
}

#[test]
fn single_end_and_nonoverlapping_pairs_keep_their_existing_results() {
    let paired = paired_record("pair", 99, 101, 201, 60, "8M", "Z.Z.Z.Z.")
        + &record("other-a", 0, 101, 60, "8M", "z.z.z.z.")
        + &paired_record("pair", 147, 201, 101, 60, "8M", "Z.Z.Z.Z.")
        + &record("other-b", 0, 201, 60, "8M", "z.z.z.z.");
    let single_end =
        record("a", 0, 101, 60, "8M", "Z.Z.Z.Z.") + &record("b", 0, 101, 60, "8M", "z.z.z.z.");
    for records in [paired, single_end] {
        for metric in ["pdr", "mhl", "lpmd", "pm", "me", "fdrp", "qfdrp"] {
            let w = Workspace::new();
            let input = w.sam(&records);
            let mut args = if metric == "lpmd" {
                vec![]
            } else {
                vec!["-d", "2"]
            };
            if matches!(metric, "fdrp" | "qfdrp") {
                args.extend(["-l", "1"]);
            }
            w.command(metric, &input, &args).assert().success();
            let mut original = w.rows();
            original.sort();
            args.push("--ignore-overlaps");
            w.command(metric, &input, &args).assert().success();
            let mut corrected = w.rows();
            corrected.sort();
            assert!(!original.is_empty());
            assert_eq!(original, corrected, "{metric}");
        }
    }
}

#[test]
fn overlap_correction_rejects_unsorted_input_without_overwriting_output() {
    let w = Workspace::new();
    let input = w.sam(
        &(record("later", 0, 201, 60, "8M", "Z.Z.Z.Z.")
            + &record("earlier", 0, 101, 60, "8M", "z.z.z.z.")),
    );
    let output = w.write("output.tsv", "previous successful result\n");
    w.command("pm", &input, &["-d", "1", "--ignore-overlaps"])
        .assert()
        .failure()
        .stderr(predicates::str::contains("coordinate-sorted"));
    assert_eq!(
        std::fs::read_to_string(output).unwrap(),
        "previous successful result\n"
    );
}

#[test]
fn minimum_depth_is_applied_after_overlap_correction() {
    for metric in ["pdr", "mhl", "pm", "me", "fdrp", "qfdrp"] {
        let w = Workspace::new();
        let input = w.sam(&overlapping_fragments());
        let mut args = vec!["-d", "3"];
        if matches!(metric, "fdrp" | "qfdrp") {
            args.extend(["-l", "1"]);
        }
        w.command(metric, &input, &args).assert().success();
        assert!(!w.rows().is_empty());
        args.push("--ignore-overlaps");
        w.command(metric, &input, &args).assert().success();
        assert!(w.rows().is_empty(), "{}", metric);
    }
}

#[test]
fn unmatched_and_nonprimary_mates_keep_their_existing_contributions() {
    let first = paired_record("pair", 99, 101, 101, 60, "8M", "Z.Z.Z.Z.");
    let second = paired_record("pair", 147, 101, 101, 60, "8M", "Z.Z.Z.Z.");
    let other = record("other", 0, 101, 60, "8M", "z.z.z.z.");
    let cases = [
        (first.clone() + &other, 2),
        (
            paired_record("pair", 99 | 8, 101, 101, 60, "8M", "Z.Z.Z.Z.") + &second + &other,
            3,
        ),
        (
            first.clone()
                + &paired_record("pair", 147 | 256, 101, 101, 60, "8M", "Z.Z.Z.Z.")
                + &second
                + &other,
            3,
        ),
        (
            first
                + &paired_record("pair", 147 | 2048, 101, 101, 60, "8M", "Z.Z.Z.Z.")
                + &second
                + &other,
            3,
        ),
    ];
    for (records, depth) in cases {
        let w = Workspace::new();
        let input = w.sam(&records);
        w.command("pdr", &input, &["-d", "1", "--ignore-overlaps"])
            .assert()
            .success();
        let rows = w.rows();
        assert_eq!(rows.len(), 4);
        for row in rows {
            assert_eq!(row[4].parse::<u32>().unwrap(), depth);
        }
    }
}

#[test]
fn an_unmatched_mate_does_not_interfere_with_later_fragments() {
    let w = Workspace::new();
    let input = w.sam(
        &(paired_record("missing-mate", 99, 101, 103, 60, "8M", "Z.Z.Z.Z.")
            + &paired_record("later", 99, 111, 113, 60, "8M", "Z.Z.Z.Z.")
            + &record("other", 0, 111, 60, "8M", "z.z.z.z.")
            + &paired_record("later", 147, 113, 111, 60, "8M", "Z.Z.Z.Z.")),
    );
    w.command("pdr", &input, &["-d", "2", "--ignore-overlaps"])
        .assert()
        .success();
    let rows = w.rows();
    assert_eq!(rows.len(), 4);
    assert_eq!(rows[0][1], "110");
    for row in rows {
        assert_eq!(&row[4..], &["2", "0"]);
    }
}

#[test]
fn many_fragments_at_one_coordinate_remain_independent() {
    let w = Workspace::new();
    let records: String = (0..1024)
        .map(|i| {
            paired_record(&format!("fragment-{i}"), 99, 101, 101, 60, "8M", "Z.Z.Z.Z.")
                + &paired_record(
                    &format!("fragment-{i}"),
                    147,
                    101,
                    101,
                    60,
                    "8M",
                    "Z.Z.Z.Z.",
                )
        })
        .collect();
    let input = w.sam(&records);
    w.command("pdr", &input, &["-d", "1", "--ignore-overlaps"])
        .assert()
        .success();
    let rows = w.rows();
    assert_eq!(rows.len(), 4);
    for row in rows {
        assert_eq!(&row[4..], &["1024", "0"]);
    }
}
