// bench/query_benchmark.rs

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::process::Command;
use std::time::Duration;

fn run_tabix_query(region: &str) -> std::io::Result<String> {
    let output = Command::new("tabix")
        .arg("data/test.bed.bgz")
        .arg(region)
        .output()?;
    Ok(String::from_utf8_lossy(&output.stdout).into_owned())
}

fn run_hgindex_query(region: &str) -> std::io::Result<String> {
    let output = Command::new("./target/release/hgidx")
        .arg("query")
        .arg(region) // region goes here
        .arg("data/test.bed.hgidx") // path to hgidx file
        .output()?;
    Ok(String::from_utf8_lossy(&output.stdout).into_owned())
}

fn query_benchmark(c: &mut Criterion) {
    let test_regions = vec!["chr1:1000-2000", "chr1:10000-20000", "chr1:100000-200000"];

    // Add output validation before benchmarking
    for region in &test_regions {
        let tabix_out = run_tabix_query(region).expect("Failed to run tabix query");
        let hgidx_out = run_hgindex_query(region).expect("Failed to run hgindex query");

        let mut tabix_lines: Vec<&str> = tabix_out.lines().filter(|l| !l.is_empty()).collect();
        let mut hgidx_lines: Vec<&str> = hgidx_out.lines().filter(|l| !l.is_empty()).collect();
        tabix_lines.sort();
        hgidx_lines.sort();
        // dbg!(&hgidx_lines);

        assert_eq!(
            tabix_lines,
            hgidx_lines,
            "\nOutput mismatch for {}:\nTabix ({} lines):\n{:#?}\nHgIndex ({} lines):\n{:#?}",
            region,
            tabix_lines.len(),
            tabix_lines,
            hgidx_lines.len(),
            hgidx_lines
        );
    }

    let mut group = c.benchmark_group("genomic_query");
    group.measurement_time(Duration::from_secs(20));
    group.sample_size(10);
    for region in test_regions {
        group.bench_with_input(BenchmarkId::new("tabix", region), &region, |b, region| {
            b.iter(|| run_tabix_query(region))
        });
        group.bench_with_input(BenchmarkId::new("hgindex", region), &region, |b, region| {
            b.iter(|| run_hgindex_query(region))
        });
    }
    group.finish();
}

criterion_group!(benches, query_benchmark);
criterion_main!(benches);
