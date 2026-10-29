mod common;

use common::{assert_close, record, Workspace};

const METRICS: [&str; 7] = ["pdr", "lpmd", "mhl", "pm", "me", "fdrp", "qfdrp"];

fn args(metric: &str) -> Vec<&str> {
    match metric {
        "lpmd" => vec![],
        "fdrp" | "qfdrp" => vec!["-d", "2", "-l", "1"],
        _ => vec!["-d", "2"],
    }
}

#[test]
fn mapping_quality_includes_the_threshold_and_excludes_the_read_below_it() {
    let records = record("below", 0, 101, 9, "8M", "Z.Z.Z.Z.")
        + &record("equal", 0, 101, 10, "8M", "z.z.z.z.")
        + &record("above", 0, 101, 11, "8M", "z.z.z.z.");
    for metric in METRICS {
        let w = Workspace::new();
        let input = w.sam(&records);
        let mut args = args(metric);
        args.extend(["-q", "10"]);
        w.command(metric, &input, &args).assert().success();
        let rows = w.rows();
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
        for row in data {
            assert_close(row.last().unwrap().parse().unwrap(), 0.0);
            if metric == "pdr" {
                assert_eq!(&row[3..], &["0", "2", "0"]);
            }
        }
    }
}

#[test]
fn minimum_depth_is_inclusive_and_a_higher_threshold_removes_the_locus() {
    for metric in ["pdr", "mhl", "pm", "me", "fdrp", "qfdrp"] {
        let w = Workspace::new();
        let input = w.sam(
            &(record("a", 0, 101, 60, "8M", "Z.Z.Z.Z.")
                + &record("b", 0, 101, 60, "8M", "Z.Z.Z.Z.")),
        );
        for (depth, count) in [
            ("2", if matches!(metric, "pm" | "me") { 1 } else { 4 }),
            ("3", 0),
        ] {
            let mut args = vec!["--min-depth", depth];
            if matches!(metric, "fdrp" | "qfdrp") {
                args.extend(["-l", "1"]);
            }
            w.command(metric, &input, &args).assert().success();
            assert_eq!(w.rows().len(), count, "metric={metric}, depth={depth}");
        }
    }
}

#[test]
fn minimum_cpg_count_is_inclusive() {
    for metric in ["pdr", "mhl"] {
        let w = Workspace::new();
        let input = w.sam(&record("a", 0, 101, 60, "8M", "Z.Z.Z.Z."));
        for (min_cpgs, expected_rows) in [("4", 4), ("5", 0)] {
            w.command(metric, &input, &["-d", "1", "-p", min_cpgs])
                .assert()
                .success();
            assert_eq!(w.rows().len(), expected_rows);
        }
    }
}

#[test]
fn chromosome_transition_and_eof_both_keep_the_last_locus() {
    for metric in ["pdr", "mhl", "pm", "me", "fdrp", "qfdrp"] {
        let w = Workspace::new();
        let chr1 =
            record("a", 0, 101, 60, "8M", "Z.Z.Z.Z.") + &record("b", 0, 101, 60, "8M", "Z.Z.Z.Z.");
        let chr2 = chr1.replace("chr1", "chr2");
        let input = w.sam(&(chr1 + &chr2));
        w.command(metric, &input, &args(metric)).assert().success();
        let rows = w.rows();
        let count = if matches!(metric, "pm" | "me") { 1 } else { 4 };
        for chrom in ["chr1", "chr2"] {
            assert_eq!(
                rows.iter().filter(|r| r[0] == chrom).count(),
                count,
                "{metric}"
            );
        }
    }
}

#[test]
fn valid_empty_alignment_has_no_loci_and_lpmd_is_undefined() {
    for metric in METRICS {
        let w = Workspace::new();
        let input = w.sam("");
        w.command(metric, &input, &args(metric)).assert().success();
        let rows = w.rows();
        if metric == "lpmd" {
            assert_eq!(rows.len(), 2);
            assert_eq!(rows[0], ["name", "lpmd"]);
            assert!(rows[1][1].parse::<f64>().unwrap().is_nan());
        } else {
            assert!(rows.is_empty(), "{}", metric);
        }
    }
}

#[test]
fn reservoir_sampling_preserves_homogeneous_read_results() {
    for metric in ["fdrp", "qfdrp"] {
        let w = Workspace::new();
        let records: String = (0..80)
            .map(|i| record(&format!("r{i}"), 0, 101, 60, "8M", "Z.Z.Z.Z."))
            .collect();
        let input = w.sam(&records);
        w.command(metric, &input, &["-d", "2", "-D", "3", "-l", "1"])
            .assert()
            .success();
        let rows = w.rows();
        assert_eq!(rows.len(), 4);
        for row in rows {
            assert_close(row[3].parse().unwrap(), 0.0);
        }
    }
}

#[test]
fn target_bed_filters_before_computing_read_patterns() {
    for metric in METRICS {
        let w = Workspace::new();
        let input = w.sam(
            &(record("a", 0, 101, 60, "8M", "Z.z.Z.z.")
                + &record("b", 0, 101, 60, "8M", "z.Z.z.Z.")),
        );
        let bed = w.write(
            "targets.bed",
            "chr1\t100\t102\nchr1\t104\t106\nchr2\t102\t104\n",
        );
        let mut args = args(metric);
        args.extend(["-c", bed.to_str().unwrap()]);
        if matches!(metric, "pdr" | "mhl") {
            args.extend(["-p", "2"]);
        }
        w.command(metric, &input, &args).assert().success();
        let rows = w.rows();
        match metric {
            "pm" | "me" => assert!(rows.is_empty()),
            "lpmd" => {
                assert_eq!(rows.len(), 2);
                assert_close(rows[1][1].parse().unwrap(), 0.0);
            }
            _ => {
                assert_eq!(rows.len(), 2);
                assert_eq!(
                    rows.iter().map(|r| r[1].as_str()).collect::<Vec<_>>(),
                    vec!["100", "104"]
                );
                let expected = match metric {
                    "mhl" => 0.5,
                    "fdrp" | "qfdrp" => 1.0,
                    _ => 0.0,
                };
                for row in rows {
                    assert_close(row[3].parse().unwrap(), expected);
                }
            }
        }
    }
}

#[test]
fn lpmd_summary_weights_pair_observations_and_exports_exact_pair_counts() {
    let w = Workspace::new();
    let input =
        w.sam(&(record("a", 0, 101, 60, "6M", "Z.Z.z.") + &record("b", 0, 101, 60, "4M", "Z.Z.")));
    let pairs = w.path("pairs.tsv");
    w.command(
        "lpmd",
        &input,
        &["-m", "2", "-M", "4", "--pairs", pairs.to_str().unwrap()],
    )
    .assert()
    .success();
    assert_close(w.rows()[1][1].parse().unwrap(), 0.5);
    assert_eq!(std::fs::read_to_string(pairs).unwrap(),
        "chrom\tcpg1\tcpg2\tlpmd\tn_concordant\tn_discordant\nchr1\t100\t102\t0\t2\t0\nchr1\t100\t104\t1\t0\t1\nchr1\t102\t104\t1\t0\t1\n");
}

#[test]
fn bam_and_sam_inputs_produce_identical_metric_results() {
    use rust_htslib::{bam, bam::Read};
    for metric in METRICS {
        let w = Workspace::new();
        let sam = w.sam(
            &(record("a", 0, 101, 60, "8M", "Z.Z.Z.Z.")
                + &record("b", 0, 101, 60, "8M", "z.z.z.z.")),
        );
        let bam = w.path("input.bam");
        let mut reader = bam::Reader::from_path(&sam).unwrap();
        let mut writer = bam::Writer::from_path(
            &bam,
            &bam::Header::from_template(reader.header()),
            bam::Format::Bam,
        )
        .unwrap();
        for r in reader.records() {
            writer.write(&r.unwrap()).unwrap();
        }
        drop(writer);
        w.command(metric, &sam, &args(metric)).assert().success();
        let mut sam_rows = w.rows();
        w.command(metric, &bam, &args(metric)).assert().success();
        let mut bam_rows = w.rows();
        // LPMD names its input; quartet row ordering is not part of the contract.
        if metric == "lpmd" {
            sam_rows[1][0].clear();
            bam_rows[1][0].clear();
        }
        sam_rows.sort();
        bam_rows.sort();
        assert_eq!(sam_rows, bam_rows, "{metric}");
    }
}

#[test]
fn invalid_xm_is_rejected_before_replacing_an_existing_result() {
    for metric in METRICS {
        for tag in ["\tXM:i:3", "\tXM:Z:Z."] {
            let w = Workspace::new();
            let input = w.sam(&format!("r\t0\tchr1\t101\t60\t4M\t*\t0\t0\tCGCG\t*{tag}\n"));
            let output = w.write("output.tsv", "previous result\n");
            w.command(metric, &input, &args(metric)).assert().failure();
            assert_eq!(
                std::fs::read_to_string(output).unwrap(),
                "previous result\n"
            );
        }
    }
}

#[test]
fn invalid_target_bed_is_rejected_before_replacing_an_existing_result() {
    for metric in METRICS {
        for content in ["chr1\n", "chr1\tnot-a-position\t102\n", "chr3\t100\t102\n"] {
            let w = Workspace::new();
            let input = w.sam(&record("a", 0, 101, 60, "8M", "Z.Z.Z.Z."));
            let bed = w.write("invalid.bed", content);
            let output = w.write("output.tsv", "previous result\n");
            let mut args = args(metric);
            args.extend(["-c", bed.to_str().unwrap()]);
            w.command(metric, &input, &args).assert().failure();
            assert_eq!(
                std::fs::read_to_string(output).unwrap(),
                "previous result\n"
            );
        }
    }
}
