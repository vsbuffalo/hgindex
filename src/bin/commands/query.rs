// bin/commands/query.rs

use crate::commands::BedRecord;
use clap::Args;
use csv::ReaderBuilder;
use flate2::Compression;
use hgindex::error::HgIndexError;
use hgindex::io::OutputStream;
use hgindex::store::GenomicDataStore;
use std::fmt::Write;
use std::fs;
use std::fs::File;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Args)]
pub struct QueryArgs {
    /// Output file.
    #[arg(short, long, value_name = "overlaps.bed")]
    pub output: Option<String>,

    /// The query region, in the format seqname:start-end where start and end are
    /// 1-based inclusive coordinates (like tabix's region argument).
    #[arg(
        value_name = "chr17:7661779-7687538",
        required_unless_present = "regions"
    )]
    pub region: Option<String>,

    /// Input BED file for batch queries
    #[arg(long, value_name = "regions.bed", required_unless_present = "region")]
    pub regions: Option<PathBuf>,

    /// Input .hgidx directory. If not specified, a file with the suffix .hgidx
    /// will be looked for in the current directory. If a single match is found,
    /// it will be used.
    #[arg(short, long, value_name = "scores.hgidx")]
    pub input: Option<PathBuf>,
}

pub fn run(args: QueryArgs) -> Result<(), HgIndexError> {
    let duration_start = Instant::now();

    // Builder output file, possibly compressed
    let output_stream = OutputStream::builder()
        .filepath(args.output)
        .buffer_size(256 * 1024)
        .compression_level(None::<Compression>)
        .build();
    let mut output_writer = output_stream.writer()?;

    // Determine input path
    let input_path = match args.input {
        Some(path) => path,
        None => find_default_hgidx_file()?,
    };

    // Verify the input path exists
    if !input_path.exists() {
        return Err(format!("Input file {} does not exist.", input_path.display()).into());
    }

    // Open store once for all queries
    let mut store = GenomicDataStore::<BedRecord, ()>::open(&input_path, None)?;

    if let Some(region) = args.region {
        // Single region query
        eprintln!("Query region {} in {}", region, input_path.display());
        query_single_region(&mut store, &region, &mut output_writer)?;
    } else if let Some(regions_file) = args.regions {
        // Batch query from BED file
        eprintln!(
            "Querying regions from {} in {}",
            regions_file.display(),
            input_path.display()
        );
        query_bed_regions(&mut store, &regions_file, &mut output_writer)?;
    }

    let duration = duration_start.elapsed();
    eprintln!("Query completed in {:?}", duration);
    Ok(())
}

fn query_single_region<W: std::io::Write>(
    store: &mut GenomicDataStore<BedRecord, ()>,
    region: &str,
    output_writer: &mut W,
) -> Result<(), HgIndexError> {
    let (seqname, start, end) = parse_region(region)?;
    // let records = store.get_overlapping(seqname, start, end)?;
    let records = store.get_overlapping_batch(seqname, start, end)?;

    // Re-usable line buffer
    let mut line_buffer = String::new();

    eprintln!("{} records found.", records.len());

    // Write all matching records to output at once
    for record in records {
        line_buffer.clear(); // Clear the buffer for the next record
        write!(
            line_buffer,
            "{}\t{}\t{}\t{}",
            record.chrom, record.start, record.end, record.rest
        )
        .unwrap();
        writeln!(output_writer, "{}", line_buffer)?; // Write the record
    }

    Ok(())
}

fn query_bed_regions<W: std::io::Write>(
    store: &mut GenomicDataStore<BedRecord, ()>,
    regions_file: &PathBuf,
    output_writer: &mut W,
) -> Result<(), HgIndexError> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .comment(Some(b'#'))
        .from_reader(File::open(regions_file)?);

    let mut total_records = 0;
    // Re-usable line buffer
    let mut line_buffer = String::new();

    for record in reader.records() {
        let record = record?;
        let chrom = record.get(0).ok_or("Missing chrom")?;
        // Convert to 0-based coordinates
        let start: u32 = record
            .get(1)
            .ok_or("Missing start")?
            .parse()
            .map_err(|_| "Invalid start coordinate")?;
        let end: u32 = record
            .get(2)
            .ok_or("Missing end")?
            .parse()
            .map_err(|_| "Invalid end coordinate")?;

        let records = store.get_overlapping_batch(chrom, start, end)?;
        total_records += records.len();

        // Write all matching records to output at once
        for record in records {
            line_buffer.clear(); // Clear the buffer for the next record
            write!(
                line_buffer,
                "{}\t{}\t{}\t{}",
                record.chrom, record.start, record.end, record.rest
            )
            .unwrap();
            writeln!(output_writer, "{}", line_buffer)?; // Write the record
        }
    }

    eprintln!("{} total records found.", total_records);
    Ok(())
}

fn parse_region(region: &str) -> Result<(&str, u32, u32), HgIndexError> {
    let region_parts: Vec<&str> = region.split(':').collect();
    if region_parts.len() != 2 {
        return Err("Invalid region format. Expected seqname:start-end.".into());
    }

    let seqname = region_parts[0];
    let coords: Vec<&str> = region_parts[1].split('-').collect();
    if coords.len() != 2 {
        return Err("Invalid region format. Expected start-end.".into());
    }

    let tabix_start: u32 = coords[0].parse().map_err(|_| "Invalid start coordinate.")?;
    let tabix_end: u32 = coords[1].parse().map_err(|_| "Invalid end coordinate.")?;

    // Convert to 0-based exclusive coordinates
    let start = tabix_start
        .checked_sub(1)
        .ok_or("Start coordinate must be greater than 0")?;
    let end = tabix_end; // End remains the same as it's exclusive in 0-based

    Ok((seqname, start, end))
}

/// Utility function to find a .hgidx file in the current directory
fn find_default_hgidx_file() -> Result<PathBuf, Box<dyn std::error::Error>> {
    let current_dir = std::env::current_dir()?;
    let entries = fs::read_dir(current_dir)?;

    let mut hgidx_files: Vec<PathBuf> = Vec::new();
    for entry in entries {
        let entry = entry?;
        let path = entry.path();
        if path.extension().map(|e| e == "hgidx").unwrap_or(false) {
            hgidx_files.push(path);
        }
    }

    // If exactly one .hgidx file is found, return it
    if hgidx_files.len() == 1 {
        Ok(hgidx_files[0].clone())
    } else if hgidx_files.is_empty() {
        Err("No .hgidx file found in the current directory.".into())
    } else {
        Err("Multiple .hgidx files found, please specify one.".into())
    }
}
