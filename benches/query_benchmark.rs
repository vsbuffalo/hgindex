// bench/query_benchmark.rs

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::process::Command;
use std::time::Duration;

fn run_tabix_query(file: &str, region: &str) -> std::io::Result<String> {
    let output = Command::new("tabix")
        .arg("data/test.bed.bgz") // Use .bgz for tabix
        .arg(region)
        .output()?;
    Ok(String::from_utf8_lossy(&output.stdout).into_owned())
}

fn run_hgindex_query(file: &str, region: &str) -> std::io::Result<String> {
    let output = Command::new("./target/release/hgidx")
        .arg("query")
        .arg("data/test.bed") // .gz does not work! TODO
        .arg(region)
        .output()?;
    Ok(String::from_utf8_lossy(&output.stdout).into_owned())
}

fn query_benchmark(c: &mut Criterion) {
    let test_regions = vec!["chr1:1000-2000", "chr1:10000-20000", "chr1:100000-200000"];

    let mut group = c.benchmark_group("genomic_query");
    group.measurement_time(Duration::from_secs(20));
    group.sample_size(10);

    for region in test_regions {
        group.bench_with_input(BenchmarkId::new("tabix", region), &region, |b, region| {
            b.iter(|| run_tabix_query("data/test.bed.gz", region))
        });

        group.bench_with_input(BenchmarkId::new("hgindex", region), &region, |b, region| {
            b.iter(|| run_hgindex_query("data/test.bed.hgidx", region))
        });
    }
    group.finish();
}

criterion_group!(benches, query_benchmark);
criterion_main!(benches);
