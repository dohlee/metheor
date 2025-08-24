use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::path::PathBuf;
use tempfile::NamedTempFile;

mod data_generator;
use data_generator::*;

fn setup_output_file() -> (NamedTempFile, PathBuf) {
    let output_file = NamedTempFile::new().expect("Failed to create temp output file");
    let output_path = output_file.path().to_path_buf();
    (output_file, output_path)
}

/// Configure criterion to be less verbose
fn create_benchmark_group<'a>(
    name: &str,
    c: &'a mut Criterion,
) -> criterion::BenchmarkGroup<'a, criterion::measurement::WallTime> {
    let mut group = c.benchmark_group(name);
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(6));
    group.warm_up_time(std::time::Duration::from_secs(1));
    group.noise_threshold(0.05);
    group
}

fn benchmark_pdr_clean(c: &mut Criterion) {
    let mut group = create_benchmark_group("PDR", c);

    let datasets = vec![
        (
            "small",
            generate_small_dataset().expect("Failed to generate small dataset"),
        ),
        (
            "medium",
            generate_medium_dataset().expect("Failed to generate medium dataset"),
        ),
    ];

    for (name, data) in datasets {
        group.bench_function(BenchmarkId::from_parameter(name), |b| {
            b.iter(|| {
                let (_output_file, output_path) = setup_output_file();
                let output_str = output_path.to_str().unwrap();
                metheor::pdr::compute(&data.bam_path, output_str, 10, 4, 10, &None);
            });
        });
    }

    group.finish();
}

fn benchmark_lpmd_clean(c: &mut Criterion) {
    let mut group = create_benchmark_group("LPMD", c);

    let datasets = vec![
        (
            "small",
            generate_small_dataset().expect("Failed to generate small dataset"),
        ),
        (
            "medium",
            generate_medium_dataset().expect("Failed to generate medium dataset"),
        ),
    ];

    for (name, data) in datasets {
        group.bench_function(BenchmarkId::from_parameter(name), |b| {
            b.iter(|| {
                let (_output_file, output_path) = setup_output_file();
                let output_str = output_path.to_str().unwrap();
                metheor::lpmd::compute(&data.bam_path, output_str, 2, 16, 10, &None, &None);
            });
        });
    }

    group.finish();
}

fn benchmark_mhl_clean(c: &mut Criterion) {
    let mut group = create_benchmark_group("MHL", c);

    let datasets = vec![
        (
            "small",
            generate_small_dataset().expect("Failed to generate small dataset"),
        ),
        (
            "medium",
            generate_medium_dataset().expect("Failed to generate medium dataset"),
        ),
    ];

    for (name, data) in datasets {
        group.bench_function(BenchmarkId::from_parameter(name), |b| {
            b.iter(|| {
                let (_output_file, output_path) = setup_output_file();
                let output_str = output_path.to_str().unwrap();
                metheor::mhl::compute(&data.bam_path, output_str, 10, 4, 10, &None);
            });
        });
    }

    group.finish();
}

fn benchmark_pm_clean(c: &mut Criterion) {
    let mut group = create_benchmark_group("PM", c);

    let datasets = vec![
        (
            "small",
            generate_small_dataset().expect("Failed to generate small dataset"),
        ),
        (
            "medium",
            generate_medium_dataset().expect("Failed to generate medium dataset"),
        ),
    ];

    for (name, data) in datasets {
        group.bench_function(BenchmarkId::from_parameter(name), |b| {
            b.iter(|| {
                let (_output_file, output_path) = setup_output_file();
                let output_str = output_path.to_str().unwrap();
                metheor::pm::compute(&data.bam_path, output_str, 10, 10, &None);
            });
        });
    }

    group.finish();
}

fn benchmark_me_clean(c: &mut Criterion) {
    let mut group = create_benchmark_group("ME", c);

    let datasets = vec![
        (
            "small",
            generate_small_dataset().expect("Failed to generate small dataset"),
        ),
        (
            "medium",
            generate_medium_dataset().expect("Failed to generate medium dataset"),
        ),
    ];

    for (name, data) in datasets {
        group.bench_function(BenchmarkId::from_parameter(name), |b| {
            b.iter(|| {
                let (_output_file, output_path) = setup_output_file();
                let output_str = output_path.to_str().unwrap();
                metheor::me::compute(&data.bam_path, output_str, 10, 10, &None);
            });
        });
    }

    group.finish();
}

fn benchmark_quick_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("quick_comparison");
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(5));
    group.warm_up_time(std::time::Duration::from_secs(1));

    let data = generate_small_dataset().expect("Failed to generate dataset");

    // Individual benchmarks for each measure
    group.bench_function(BenchmarkId::from_parameter("PDR"), |b| {
        b.iter(|| {
            let (_output_file, output_path) = setup_output_file();
            let output_str = output_path.to_str().unwrap();
            metheor::pdr::compute(&data.bam_path, output_str, 10, 4, 10, &None);
        });
    });

    group.bench_function(BenchmarkId::from_parameter("LPMD"), |b| {
        b.iter(|| {
            let (_output_file, output_path) = setup_output_file();
            let output_str = output_path.to_str().unwrap();
            metheor::lpmd::compute(&data.bam_path, output_str, 2, 16, 10, &None, &None);
        });
    });

    group.bench_function(BenchmarkId::from_parameter("MHL"), |b| {
        b.iter(|| {
            let (_output_file, output_path) = setup_output_file();
            let output_str = output_path.to_str().unwrap();
            metheor::mhl::compute(&data.bam_path, output_str, 10, 4, 10, &None);
        });
    });

    group.bench_function(BenchmarkId::from_parameter("PM"), |b| {
        b.iter(|| {
            let (_output_file, output_path) = setup_output_file();
            let output_str = output_path.to_str().unwrap();
            metheor::pm::compute(&data.bam_path, output_str, 10, 10, &None);
        });
    });

    group.bench_function(BenchmarkId::from_parameter("ME"), |b| {
        b.iter(|| {
            let (_output_file, output_path) = setup_output_file();
            let output_str = output_path.to_str().unwrap();
            metheor::me::compute(&data.bam_path, output_str, 10, 10, &None);
        });
    });

    group.bench_function(BenchmarkId::from_parameter("FDRP"), |b| {
        b.iter(|| {
            let (_output_file, output_path) = setup_output_file();
            let output_str = output_path.to_str().unwrap();
            metheor::fdrp::compute(&data.bam_path, output_str, 10, 10, 40, 35, &None);
        });
    });

    group.bench_function(BenchmarkId::from_parameter("qFDRP"), |b| {
        b.iter(|| {
            let (_output_file, output_path) = setup_output_file();
            let output_str = output_path.to_str().unwrap();
            metheor::qfdrp::compute(&data.bam_path, output_str, 10, 10, 40, 35, &None);
        });
    });

    group.finish();
}

criterion_group!(
    benches_clean,
    benchmark_pdr_clean,
    benchmark_lpmd_clean,
    benchmark_mhl_clean,
    benchmark_pm_clean,
    benchmark_me_clean,
    benchmark_quick_comparison
);

criterion_main!(benches_clean);
