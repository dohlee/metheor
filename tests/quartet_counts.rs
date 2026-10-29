mod common;

use common::{assert_close, record, Workspace};

#[test]
fn pm_can_append_counts_in_binary_epiallele_order_without_changing_default_output() {
    let patterns = [
        "z.z.z.z.", "z.z.z.Z.", "z.z.Z.z.", "z.z.Z.Z.", "z.Z.z.z.", "z.Z.z.Z.", "z.Z.Z.z.",
        "z.Z.Z.Z.", "Z.z.z.z.", "Z.z.z.Z.", "Z.z.Z.z.", "Z.z.Z.Z.", "Z.Z.z.z.", "Z.Z.z.Z.",
        "Z.Z.Z.z.", "Z.Z.Z.Z.",
    ];
    let mut records: String = patterns
        .iter()
        .enumerate()
        .map(|(i, pattern)| record(&format!("p{i}"), 0, 101, 60, "8M", pattern))
        .collect();
    records += &record("zero1", 0, 101, 60, "8M", patterns[0]);
    records += &record("zero2", 0, 101, 60, "8M", patterns[0]);
    for i in 0..3 {
        records += &record(&format!("full{i}"), 0, 101, 60, "8M", patterns[15]);
    }
    let w = Workspace::new();
    let input = w.sam(&records);
    w.command("pm", &input, &["-d", "1"]).assert().success();
    let default = w.rows();
    assert_eq!(default.len(), 1);
    assert_eq!(default[0].len(), 6);
    assert_close(default[0][5].parse().unwrap(), 0.9115646258503401);
    w.command("pm", &input, &["-d", "1", "--include-counts"])
        .assert()
        .success();
    let rows = w.rows();
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].len(), 22);
    assert_eq!(&rows[0][..6], default[0]);
    let counts: Vec<u32> = rows[0][6..]
        .iter()
        .map(|count| count.parse().unwrap())
        .collect();
    assert_eq!(counts, [3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4]);
}

#[test]
fn exported_counts_respect_mapping_quality_targets_and_quartet_depth() {
    let w = Workspace::new();
    let records = record("one", 0, 101, 10, "10M", "Z.z.Z.z.Z.")
        + &record("two", 0, 101, 60, "10M", "Z.z.Z.z.Z.")
        + &record("low", 0, 101, 9, "10M", "z.Z.z.Z.z.");
    let input = w.sam(&records);
    let bed = w.write(
        "targets.bed",
        "chr1\t100\t101\nchr1\t104\t105\nchr1\t106\t107\nchr1\t108\t109\n",
    );
    w.command(
        "pm",
        &input,
        &[
            "--include-counts",
            "-q",
            "10",
            "-d",
            "2",
            "-c",
            bed.to_str().unwrap(),
        ],
    )
    .assert()
    .success();
    let rows = w.rows();
    assert_eq!(rows.len(), 1);
    assert_eq!(&rows[0][..6], &["chr1", "100", "104", "106", "108", "0"]);
    let counts: Vec<u32> = rows[0][6..]
        .iter()
        .map(|value| value.parse().unwrap())
        .collect();
    assert_eq!(counts, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0]);
    w.command(
        "pm",
        &input,
        &[
            "--include-counts",
            "-q",
            "10",
            "-d",
            "3",
            "-c",
            bed.to_str().unwrap(),
        ],
    )
    .assert()
    .success();
    assert!(w.rows().is_empty());
}
