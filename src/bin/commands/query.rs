// bin/commands/query.rs

use clap::Args;
use flate2::Compression;
use hgindex::error::HgIndexError;
use hgindex::io::OutputStream;
use hgindex::store::GenomicDataStore;
use hgindex::{BedRecord, BedRecordSlice};
use itoa;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use crate::commands::pack::build_tsv_reader;

#[derive(Args)]
pub struct QueryArgs {
    /// Output file.
    #[arg(short, long, value_name = "overlaps.bed")]
    pub output: Option<String>,

    /// Comment character to skip lines starting with this
    #[arg(long, default_value = "#")]
    pub comment: char,

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
        .buffer_size(1024 * 1024)
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
    let mut store = GenomicDataStore::<BedRecord>::open(&input_path, None)?;

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
        query_bed_regions(&mut store, &regions_file, &mut output_writer, &args.comment)?;
    }

    let duration = duration_start.elapsed();
    eprintln!("Query completed in {:?}", duration);
    Ok(())
}

fn query_single_region<W: std::io::Write>(
    store: &mut GenomicDataStore<BedRecord>,
    region: &str,
    output_writer: &mut W,
) -> Result<(), HgIndexError> {
    let (seqname, start, end) = parse_region(region)?;

    // Use `map_overlapping` for efficient ZCD
    let record_count = store.map_overlapping(seqname, start, end, |record_slice| {
        write_tsv_bytes(seqname, &record_slice, output_writer)?;
        Ok(())
    })?;

    eprintln!("{} records processed.", record_count);
    Ok(())
}

fn query_bed_regions<W: std::io::Write>(
    store: &mut GenomicDataStore<BedRecord>,
    regions_file: &PathBuf,
    output_writer: &mut W,
    comment_char: &char,
) -> Result<(), HgIndexError> {
    let mut reader = build_tsv_reader(
        regions_file,
        Some(*comment_char as u8),
        true,  // flexible
        false, // has_headers
    )?;

    let mut total_records = 0;
    // Initialize batch with reasonable starting capacity
    let mut batch = RecordBatch::with_capacity(64 * 1024);

    for record in reader.records() {
        let record = record?;
        let chrom = record.get(0).ok_or("Missing chrom")?.to_string();
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

        let records = store.get_overlapping_batch(&chrom, start, end)?;
        for record in records {
            batch.push_record(&chrom, &record);
            if batch.should_flush() {
                batch.write_batch(output_writer)?;
            }
            total_records += 1;
        }
    }

    // Flush any remaining records
    if batch.records_seen > 0 {
        batch.write_batch(output_writer)?;
    }

    eprintln!("Found {} total records.", total_records);
    Ok(())
}

#[inline(always)]
fn write_tsv_bytes<W: std::io::Write>(
    chrom: &str,
    record: &BedRecordSlice<'_>,
    writer: &mut W,
) -> Result<(), HgIndexError> {
    // Directly write to the writer without intermediate buffer
    write!(writer, "{}\t{}\t{}\t", chrom, record.start, record.end)?;
    writer.write_all(record.rest)?; // Raw bytes, no conversion
    writer.write_all(b"\n")?;
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

// Struct to batch records and minimize allocations
pub struct RecordBatch {
    // Pre-allocated buffer for collecting records
    buffer: Vec<u8>,
    // Capacity tracking to avoid too many resizes
    records_seen: usize,
    // Add dedicated number buffers to avoid allocations
    start_buffer: itoa::Buffer,
    end_buffer: itoa::Buffer,
}

impl RecordBatch {
    pub fn with_capacity(bytes: usize) -> Self {
        Self {
            buffer: Vec::with_capacity(bytes),
            records_seen: 0,
            start_buffer: itoa::Buffer::new(),
            end_buffer: itoa::Buffer::new(),
        }
    }

    #[inline(always)]
    pub fn push_record(&mut self, chrom: &str, record: &BedRecordSlice<'_>) {
        // Extend chrom bytes
        self.buffer.extend_from_slice(chrom.as_bytes());
        self.buffer.push(b'\t');

        // Format integers using dedicated buffers
        let start_str = self.start_buffer.format(record.start);
        let end_str = self.end_buffer.format(record.end);

        self.buffer.extend_from_slice(start_str.as_bytes());
        self.buffer.push(b'\t');
        self.buffer.extend_from_slice(end_str.as_bytes());
        self.buffer.push(b'\t');

        // Rest of record and newline
        self.buffer.extend_from_slice(record.rest);
        self.buffer.push(b'\n');

        self.records_seen += 1;
    }

    // Flush when batch is large enough
    #[inline(always)]
    pub fn should_flush(&self) -> bool {
        self.records_seen >= 1000 || self.buffer.len() >= 64 * 1024
    }

    #[inline(always)]
    pub fn write_batch<W: std::io::Write>(&mut self, writer: &mut W) -> Result<(), HgIndexError> {
        writer.write_all(&self.buffer)?;
        self.clear();
        Ok(())
    }

    pub fn clear(&mut self) {
        // Clears the underlying Vec<u8> without deallocating memory
        self.buffer.clear();
        // Reset the record counter
        self.records_seen = 0;
        // Note: we don't need to clear itoa::Buffer as it's reused in-place
    }
}
