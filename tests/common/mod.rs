#![allow(dead_code)]

use assert_cmd::Command;
use std::{
    fs,
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
};

static NEXT_ID: AtomicU64 = AtomicU64::new(0);

pub struct Workspace(PathBuf);

impl Workspace {
    pub fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "metheor-test-{}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
            NEXT_ID.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }

    pub fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }

    pub fn write(&self, name: &str, content: &str) -> PathBuf {
        let path = self.path(name);
        fs::write(&path, content).unwrap();
        path
    }

    pub fn sam(&self, records: &str) -> PathBuf {
        self.write("input.sam", &format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n{records}"))
    }

    pub fn command(&self, metric: &str, input: &Path, args: &[&str]) -> Command {
        let mut cmd = Command::cargo_bin("metheor").unwrap();
        cmd.arg(metric)
            .arg("--input")
            .arg(input)
            .arg("--output")
            .arg(self.path("output.tsv"))
            .args(args);
        cmd
    }

    pub fn rows(&self) -> Vec<Vec<String>> {
        fs::read_to_string(self.path("output.tsv"))
            .unwrap()
            .lines()
            .map(|line| line.split('\t').map(str::to_owned).collect())
            .collect()
    }
}

impl Drop for Workspace {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

pub fn record(name: &str, flags: u16, pos: u32, mapq: u8, cigar: &str, xm: &str) -> String {
    let seq: String = xm
        .chars()
        .map(|c| match c {
            'Z' => 'C',
            'z' => 'T',
            _ => 'G',
        })
        .collect();
    format!("{name}\t{flags}\tchr1\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t*\tXM:Z:{xm}\n")
}

pub fn assert_close(actual: f64, expected: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() < 1e-6,
        "expected {}, got {}",
        expected,
        actual
    );
}
