mod common;

use common::{assert_close, record, Workspace};

#[test]
fn tied_alignment_order_does_not_lose_cpg_contributions() {
    for length in [4, 202] {
        let a = record("a", 0, 101, 60, "4M", "Z.Z.");
        let b = record(
            "b",
            0,
            102,
            60,
            &format!("{length}M"),
            &(".".repeat(length - 1) + "Z"),
        );
        let c = record("c", 0, 102, 60, "4M", ".z.Z");
        for records in [format!("{a}{b}{c}"), format!("{a}{c}{b}")] {
            for (metric, expected) in [
                ("pdr", 0.5),
                ("mhl", 0.5833333333333333),
                ("fdrp", 1.0),
                ("qfdrp", 1.0),
            ] {
                let w = Workspace::new();
                let input = w.sam(&records);
                let args = if matches!(metric, "pdr" | "mhl") {
                    vec!["-d", "2", "-p", "1"]
                } else {
                    vec!["-d", "2", "-l", "1"]
                };
                w.command(metric, &input, &args).assert().success();
                let rows = w.rows();
                let row = rows.iter().find(|r| r[1] == "102").unwrap_or_else(|| {
                    panic!(
                        "missing CpG 102: metric={}, length={}, rows={:?}",
                        metric, length, rows
                    )
                });
                assert_close(row[3].parse().unwrap(), expected);
                if metric == "pdr" {
                    assert_eq!(&row[4..], &["1", "1"]);
                }
            }
        }
    }
}
