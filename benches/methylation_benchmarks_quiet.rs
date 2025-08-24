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

/// Configure criterion to be quiet and fast
fn create_quiet_group<'a>(
    name: &str,
    c: &'a mut Criterion,
) -> criterion::BenchmarkGroup<'a, criterion::measurement::WallTime> {
    let mut group = c.benchmark_group(name);
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(5));
    group.warm_up_time(std::time::Duration::from_secs(1));
    group.noise_threshold(0.05);
    group
}

/// Benchmark the quiet measures (non-verbose output)
fn benchmark_quiet_measures(c: &mut Criterion) {
    let mut group = create_quiet_group("quiet_measures", c);

    let data = generate_small_dataset().expect("Failed to generate dataset");

    // PDR - relatively quiet
    group.bench_function(BenchmarkId::from_parameter("PDR"), |b| {
        b.iter(|| {
            let (_output_file, output_path) = setup_output_file();
            let output_str = output_path.to_str().unwrap();
            metheor::pdr::compute(&data.bam_path, output_str, 10, 4, 10, &None);
        });
    });

    // MHL - quiet
    group.bench_function(BenchmarkId::from_parameter("MHL"), |b| {
        b.iter(|| {
            let (_output_file, output_path) = setup_output_file();
            let output_str = output_path.to_str().unwrap();
            metheor::mhl::compute(&data.bam_path, output_str, 10, 4, 10, &None);
        });
    });

    // PM - quiet
    group.bench_function(BenchmarkId::from_parameter("PM"), |b| {
        b.iter(|| {
            let (_output_file, output_path) = setup_output_file();
            let output_str = output_path.to_str().unwrap();
            metheor::pm::compute(&data.bam_path, output_str, 10, 10, &None);
        });
    });

    // ME - quiet
    group.bench_function(BenchmarkId::from_parameter("ME"), |b| {
        b.iter(|| {
            let (_output_file, output_path) = setup_output_file();
            let output_str = output_path.to_str().unwrap();
            metheor::me::compute(&data.bam_path, output_str, 10, 10, &None);
        });
    });

    // FDRP - quiet
    group.bench_function(BenchmarkId::from_parameter("FDRP"), |b| {
        b.iter(|| {
            let (_output_file, output_path) = setup_output_file();
            let output_str = output_path.to_str().unwrap();
            metheor::fdrp::compute(&data.bam_path, output_str, 10, 10, 40, 35, &None);
        });
    });

    // qFDRP - quiet
    group.bench_function(BenchmarkId::from_parameter("qFDRP"), |b| {
        b.iter(|| {
            let (_output_file, output_path) = setup_output_file();
            let output_str = output_path.to_str().unwrap();
            metheor::qfdrp::compute(&data.bam_path, output_str, 10, 10, 40, 35, &None);
        });
    });

    group.finish();
}

/// Benchmark dataset size comparison for PDR only
fn benchmark_dataset_sizes(c: &mut Criterion) {
    let mut group = create_quiet_group("dataset_sizes", c);

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

/// Quick parameter variation test
fn benchmark_parameter_quick(c: &mut Criterion) {
    let mut group = create_quiet_group("parameter_quick", c);

    let data = generate_small_dataset().expect("Failed to generate dataset");

    for min_depth in &[10, 20] {
        group.bench_function(BenchmarkId::new("PDR_min_depth", min_depth), |b| {
            b.iter(|| {
                let (_output_file, output_path) = setup_output_file();
                let output_str = output_path.to_str().unwrap();
                metheor::pdr::compute(&data.bam_path, output_str, *min_depth, 4, 10, &None);
            });
        });
    }

    group.finish();
}

criterion_group!(
    benches_quiet,
    benchmark_quiet_measures,
    benchmark_dataset_sizes,
    benchmark_parameter_quick
);

criterion_main!(benches_quiet);
