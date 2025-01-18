use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::path::PathBuf;
use std::process::Command;

const QUERY_TYPES: &[(&str, &str)] = &[
    ("refgene", "repeat_masker"),
    ("refgene", "repeat_masker_autosomes"),
];

fn run_hgidx_query(regions: &str, input: &str) -> std::io::Result<std::process::Output> {
    Command::new("./target/release/hgidx")
        .arg("query")
        .arg("--regions")
        .arg(regions)
        .arg("--input")
        .arg(input)
        .output()
}

fn run_tabix_query(regions: &str, input: &str) -> std::io::Result<std::process::Output> {
    Command::new("tabix")
        .arg(input)
        .arg("--regions")
        .arg(regions)
        .output()
}

fn check_files_exist() -> Result<(), Box<dyn std::error::Error>> {
    let test_dir = PathBuf::from("tests/data");

    // Check if required files exist
    for (region_base, query_base) in QUERY_TYPES {
        // Check region files
        let region_gz = test_dir.join(format!("{}.bed.gz", region_base));
        if !region_gz.exists() {
            return Err(format!("Missing region file: {}", region_gz.display()).into());
        }

        // Check query files - both hgidx and bgz
        let query_hgidx = test_dir.join(format!("{}.hgidx", query_base));
        let query_bgz = test_dir.join(format!("{}.bed.bgz", query_base));

        if !query_hgidx.exists() {
            return Err(format!("Missing hgidx file: {}", query_hgidx.display()).into());
        }
        if !query_bgz.exists() {
            return Err(format!("Missing bgz file: {}", query_bgz.display()).into());
        }
    }

    Ok(())
}

fn bench_queries(c: &mut Criterion) {
    // Check if all required files exist
    if let Err(e) = check_files_exist() {
        eprintln!("Error: {}. Please run `make test-data` first.", e);
        return;
    }

    let mut group = c.benchmark_group("genomic_queries");
    group.sample_size(20); // Match hyperfine's min-runs
    group.warm_up_time(std::time::Duration::from_secs(10)); // Match hyperfine's warmup

    for (region_base, query_base) in QUERY_TYPES {
        let region_file = format!("tests/data/{}.bed.gz", region_base);
        let query_hgidx = format!("tests/data/{}.hgidx", query_base);
        let query_bgz = format!("tests/data/{}.bed.bgz", query_base);

        // Benchmark hgidx query
        let bench_id = BenchmarkId::new("hgidx", format!("{}_{}", region_base, query_base));
        group.bench_with_input(
            bench_id,
            &(region_file.clone(), query_hgidx),
            |b, (region, query)| {
                b.iter(|| run_hgidx_query(region, query).unwrap());
            },
        );

        // Benchmark tabix query
        let bench_id = BenchmarkId::new("tabix", format!("{}_{}", region_base, query_base));
        group.bench_with_input(bench_id, &(region_file, query_bgz), |b, (region, query)| {
            b.iter(|| run_tabix_query(region, query).unwrap());
        });
    }

    group.finish();
}

criterion_group!(benches, bench_queries);
criterion_main!(benches);
