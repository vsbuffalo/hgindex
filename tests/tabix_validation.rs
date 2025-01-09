use std::fs;
use std::path::PathBuf;
use std::process::Command;

const BEDSIZE: u32 = 1_000_000;

fn run_tabix_query(bed_file: &str, region: &str) -> Vec<String> {
    let output = Command::new("tabix")
        .arg(bed_file)
        .arg(region)
        .output()
        .expect("Failed to execute tabix");
    String::from_utf8_lossy(&output.stdout)
        .lines()
        .filter(|line| !line.trim().is_empty()) // Filter out empty lines
        .map(String::from)
        .collect()
}

fn run_hgindex_query(index_file: &str, region: &str) -> Vec<String> {
    let output = Command::new("cargo")
        .arg("run")
        .arg("--release")
        .arg("--features=cli")
        .arg("--")
        .arg("query")
        .arg(region)
        .arg(index_file)
        .output()
        .expect("Failed to execute hgindex");
    String::from_utf8_lossy(&output.stdout)
        .lines()
        .map(String::from)
        .collect()
}

#[test]
fn test_tabix_compatibility() {
    // Create a fixed test directory under target/
    let test_dir = PathBuf::from("target/test_files");
    fs::create_dir_all(&test_dir).expect("Failed to create test directory");

    // Generate paths for all files
    let test_bed = test_dir.join("test.bed");
    let bgzipped = test_dir.join("test.bed.bgz");
    let hgindex = test_dir.join("test.hgidx");

    // Print paths for inspection
    eprintln!("Test files location:");
    eprintln!("  BED file: {}", test_bed.display());
    eprintln!("  BGZip file: {}", bgzipped.display());
    eprintln!("  HGIndex file: {}", hgindex.display());

    // Generate test data
    let generate_status = Command::new("cargo")
        .arg("run")
        .arg("--release")
        .arg("--features=cli,dev")
        .arg("--")
        .arg("random-bed")
        .arg("-o")
        .arg(&test_bed)
        .arg("-n")
        .arg(BEDSIZE.to_string())
        .arg("-s")
        .arg("42")
        .status()
        .expect("Failed to generate test data");

    assert!(
        generate_status.success(),
        "Failed to generate test bed file"
    );
    assert!(test_bed.exists(), "Test bed file was not created");

    // Compress with bgzip
    let bgzip_status = Command::new("sh")
        .arg("-c")
        .arg(format!(
            "bgzip -c {} > {}",
            test_bed.display(),
            bgzipped.display()
        ))
        .status()
        .expect("Failed to execute bgzip");
    assert!(bgzip_status.success(), "bgzip command failed");
    assert!(bgzipped.exists(), "BGZipped file was not created");

    // Index with both tools
    let tabix_status = Command::new("tabix")
        .arg("-p")
        .arg("bed")
        .arg("-c")
        .arg("#")
        .arg(&bgzipped)
        .status()
        .expect("Failed to execute tabix");
    assert!(tabix_status.success(), "Tabix indexing failed");

    let hgindex_status = Command::new("cargo")
        .arg("run")
        .arg("--release")
        .arg("--features=cli")
        .arg("--")
        .arg("pack")
        .arg("-f")
        .arg(&test_bed)
        .arg("-o")
        .arg(&hgindex)
        .status()
        .expect("Failed to execute hgindex");
    assert!(hgindex_status.success(), "HGIndex packing failed");
    assert!(hgindex.exists(), "HGIndex file was not created");

    let test_regions = vec![
        // Original cases
        "chr1:15000-17000", // Spans 16kb boundary
        "chr1:15900-16100", // Small region crossing boundary
        "chr1:1-1000000",   // Large region
        "chr1:15999-16001", // Minimal region at boundary
        // Non-overlapping gaps
        "chr1:3066-3107",     // Clean gap
        "chr1:4897-4958",     // Gap between entries
        "chr1:5396-5431",     // Small gap
        "chr1:7464-7650",     // Medium gap
        "chr2:121534-121571", // Gap in chr2
        // Overlapping multiple entries
        "chr1:4000-5000",     // Dense overlapping region
        "chr1:9000-10000",    // Multiple overlaps
        "chr1:1000-2000",     // Start overlaps
        "chr2:123000-124000", // Chr2 overlaps
        // Edge cases
        "chr1:952-953",       // Minimal first entry overlap
        "chr2:134901-134902", // End of last entry
        "chr1:1-500",         // Before first entry
        "chr2:120000-120100", // Before first chr2 entry
        // Additional overlapping cases
        "chr1:1500-3000",     // Multiple entry overlap
        "chr1:5000-6000",     // Middle region overlap
        "chr2:122000-123000", // Chr2 middle overlap
        // Boundary cases
        "chr1:989-989",       // Single base overlap
        "chr2:120475-120475", // Single base in chr2
        "chr1:16395-16395",   // Single base at end
        // Large overlapping regions
        "chr1:1000-10000",    // Large overlapping region
        "chr2:120000-130000", // Large chr2 region
        // Cross-chromosome queries
        "chr1:1-2000000", // Whole chr1
        "chr2:1-2000000", // Whole chr2
        // More precise overlaps
        "chr1:2066-2147",     // Exact entry overlap
        "chr2:123105-123166", // Exact chr2 entry overlap
    ];

    for region in test_regions {
        let tabix_results = run_tabix_query(&bgzipped.to_string_lossy(), region);
        let hgindex_results = run_hgindex_query(&hgindex.to_string_lossy(), region);

        assert_eq!(
            tabix_results.len(),
            hgindex_results.len(),
            "Result count mismatch for region {}: tabix={}, hgindex={}",
            region,
            tabix_results.len(),
            hgindex_results.len()
        );

        // Compare sorted results
        let mut tabix_sorted = tabix_results.clone();
        let mut hgindex_sorted = hgindex_results.clone();
        tabix_sorted.sort();
        hgindex_sorted.sort();

        assert_eq!(
            tabix_sorted, hgindex_sorted,
            "Results differ for region {}",
            region
        );
    }
}

