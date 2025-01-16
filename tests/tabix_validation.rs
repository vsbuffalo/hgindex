use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

// const BEDSIZE: u32 = 50_000_000;
const BEDSIZE: u32 = 1_000_000;

fn check_tool_exists(tool: &str) -> Result<(), Box<dyn std::error::Error>> {
    match Command::new(tool).arg("--version").output() {
        Ok(output) => {
            println!(
                "{} version: {}",
                tool,
                String::from_utf8_lossy(&output.stdout)
            );
            Ok(())
        }
        Err(e) => Err(format!("{} not found: {}", tool, e).into()),
    }
}

#[test]
fn test_tabix_compatibility() -> Result<(), Box<dyn std::error::Error>> {
    check_tool_exists("bgzip")?;
    check_tool_exists("tabix")?;

    // Set up paths
    let test_dir = PathBuf::from("target/test_files");
    fs::create_dir_all(&test_dir)?;
    let test_dir = test_dir.canonicalize()?;
    let test_bed = test_dir.join("test.bed");
    let bgzipped = test_dir.join("test.bed.bgz");
    let hgindex = test_dir.join("test.hgidx");

    println!("Test files location:");
    println!("BED file: {}", test_bed.display());
    println!("BGZip file: {}", bgzipped.display());
    println!("HGIndex file: {}", hgindex.display());

    // Print tool versions
    let bgzip_version = Command::new("bgzip").arg("--version").output()?;
    println!(
        "BGZip version: {}",
        String::from_utf8_lossy(&bgzip_version.stdout)
    );

    let tabix_version = Command::new("tabix").arg("--version").output()?;
    println!(
        "Tabix version: {}",
        String::from_utf8_lossy(&tabix_version.stdout)
    );

    // Generate test data
    println!("Generating test BED file...");
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
        .status()?;

    if !generate_status.success() {
        return Err("Failed to generate test BED file".into());
    }

    // Verify BED file
    let bed_size = fs::metadata(&test_bed)?.len();
    println!("BED file size: {} bytes", bed_size);
    assert!(bed_size > 0, "BED file is empty");

    // Remove existing BGZip file if it exists
    if Path::new(&bgzipped).exists() {
        fs::remove_file(&bgzipped)?;
        println!("Removed existing BGZip file");
    }

    // Run BGZip compression
    println!("Running BGZip compression...");
    let bgzip_output = Command::new("sh")
        .arg("-c")
        .arg(format!(
            "bgzip -c {} > {}",
            test_bed.display(),
            bgzipped.display()
        ))
        .output()?;

    if !bgzip_output.status.success() {
        let stderr = String::from_utf8_lossy(&bgzip_output.stderr);
        println!("BGZip error: {}", stderr);
        return Err("BGZip compression failed".into());
    }

    // Verify BGZip file
    let bgzip_size = fs::metadata(&bgzipped)?.len();
    println!("BGZip file size: {} bytes", bgzip_size);
    assert!(bgzip_size > 0, "BGZip file is empty");

    // Verify file format
    let contents = fs::read(&bgzipped)?;
    if contents.len() < 2 || contents[0] != 0x1f || contents[1] != 0x8b {
        return Err("Invalid BGZip format - missing magic number".into());
    }

    // Run Tabix indexing
    println!("Running Tabix indexing...");
    let tabix_output = Command::new("tabix")
        .arg("-p")
        .arg("bed")
        .arg("-c")
        .arg("#")
        .arg(&bgzipped)
        .output()?;

    if !tabix_output.status.success() {
        let stderr = String::from_utf8_lossy(&tabix_output.stderr);
        println!("Tabix error: {}", stderr);
        return Err("Tabix indexing failed".into());
    }

    // Verify Tabix index file exists
    let index_file = format!("{}.tbi", bgzipped.display());
    assert!(
        Path::new(&index_file).exists(),
        "Tabix index file was not created"
    );

    // Run HGIndex pack
    println!("Running HGIndex pack...");
    let mut pack_command = Command::new("cargo");
    pack_command
        .arg("run")
        .arg("--release")
        .arg("--features=cli,dev")
        .arg("--")
        .arg("pack")
        .arg(&test_bed)
        .arg("--output")
        .arg(&hgindex)
        .arg("--force"); // Add force flag to overwrite if exists

    println!("Running pack command: {:?}", &pack_command);
    let pack_output = pack_command.output()?;

    if !pack_output.status.success() {
        let stderr = String::from_utf8_lossy(&pack_output.stderr);
        println!("HGIndex pack error: {}", stderr);
        return Err("HGIndex pack failed".into());
    }

    // Verify HGIndex file exists
    assert!(Path::new(&hgindex).exists(), "HGIndex file was not created");

    // Run queries with both tools and compare outputs
    let test_regions = vec![
        "chr1:15000-17000",
        "chr1:1-1000000",
        "chr2:500000-600000",
        "chrX:1000000-2000000",
    ];

    for region in test_regions {
        println!("Testing region: {}", region);
        // Get tabix results
        let tabix_output = Command::new("tabix").arg(&bgzipped).arg(region).output()?;
        if !tabix_output.status.success() {
            let stderr = String::from_utf8_lossy(&tabix_output.stderr);
            println!("Tabix query error: {}", stderr);
            return Err("Tabix query failed".into());
        }

        // Get HGIndex results
        let mut hgindex_command = Command::new("cargo");
        hgindex_command
            .arg("run")
            .arg("--release")
            .arg("--features=cli,dev")
            .arg("--")
            .arg("query")
            .arg("-i")
            .arg(&hgindex)
            .arg(region);

        println!("Running command: {:?}", &hgindex_command);
        let hgindex_output = hgindex_command.output()?;

        if !hgindex_output.status.success() {
            let stderr = String::from_utf8_lossy(&hgindex_output.stderr);
            println!("HGIndex query error: {}", stderr);
            return Err("HGIndex query failed".into());
        }

        // Compare outputs line by line to avoid any newline issues
        let tabix_data = String::from_utf8(tabix_output.stdout)?;
        let mut tabix_lines: Vec<&str> = tabix_data.lines().collect();
        let hgindex_data = String::from_utf8(hgindex_output.stdout)?;
        let mut hgindex_lines: Vec<&str> = hgindex_data.lines().collect();

        tabix_lines.sort();
        hgindex_lines.sort();

        assert_eq!(tabix_lines.len(), hgindex_lines.len());

        assert_eq!(
            tabix_lines, hgindex_lines,
            "Results don't match for region {}\nTabix output:\n{:#?}\nHGIndex output:\n{:#?}",
            region, tabix_lines, hgindex_lines
        );
    }

    Ok(())
}
