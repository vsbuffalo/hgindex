// bin/commands/pack.rs

use clap::Args;
use csv::ReaderBuilder;
use hgindex::error::HgIndexError;
use hgindex::store::GenomicDataStore;
use hgindex::{BedRecord, InputStream};
use indicatif::{ProgressBar, ProgressStyle};
use std::fs::File;
use std::io::{BufRead, Read};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Args)]
pub struct PackArgs {
    /// Input TSV/BED file to pack and index (suffix should be data.hgidx)
    #[arg(value_name = "FILE")]
    pub input: PathBuf,

    /// Output path. If not specified, will append .hgidx to input path
    #[arg(short = 'o', long)]
    pub output: Option<PathBuf>,

    /// Comment character to skip lines starting with this
    #[arg(long, default_value = "#")]
    pub comment: char,

    /// Use 1-based coordinates instead of 0-based
    #[arg(long)]
    pub one_based: bool,

    /// Force overwrite of output file if it exists
    #[arg(short = 'f', long)]
    pub force: bool,

    /// Hierarchical binning schema to use
    #[arg(long, value_enum, default_value_t = hgindex::BinningSchema::Dense)]
    pub schema: hgindex::BinningSchema,
}

pub fn run(args: PackArgs) -> Result<(), HgIndexError> {
    // For timing the pack operation
    let start = Instant::now();

    // Create the output path by stemming the path.
    let output_path = args.output.unwrap_or_else(|| {
        let name = args.input.file_stem().unwrap_or_default().to_string_lossy();
        let parent = args.input.parent().unwrap_or_else(|| Path::new("."));
        parent.join(name.to_string()).with_extension("hgidx")
    });

    // Check if output exists and handle --force
    if output_path.exists() && !args.force {
        return Err("Output file exists. Use --force to overwrite.".into());
    }

    eprintln!(
        "Packing {} to {}",
        args.input.display(),
        output_path.display()
    );

    // Create store
    eprintln!("Index binning schema: {:?}", args.schema);
    let mut store =
        GenomicDataStore::<BedRecord>::create_with_schema(&output_path, None, &args.schema)?;

    let mut csv_reader = build_tsv_reader(
        &args.input,
        Some(args.comment as u8),
        true,  // flexible
        false, // has_headers
    )?;

    // Estimate total records
    let estimated_records =
        estimate_total_records(&args.input, Some(args.comment as u8), b'\t', false, true)?;

    // Set up the progress bar.
    let pb = ProgressBar::new(estimated_records).with_style(
            ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue}⟩ {pos}/{len} ({percent}%) [{eta_precise}]")?
            .progress_chars("=> ")
        );

    // At the start, after creating the progress bar:
    #[cfg(feature = "dev")]
    let start_time = std::time::Instant::now();
    #[cfg(feature = "dev")]
    let initial_estimated_records = estimated_records; // Capture for comparison

    // Duration estimation sampling stuff
    let update_frequency = 1000;
    let mut counter = 0;

    // Process records
    for result in csv_reader.byte_records() {
        let record = result?;

        // Safe conversion of chromosome name
        let chrom = String::from_utf8_lossy(&record[0]).into_owned();

        // Parse start and end positions
        let start: u32 = String::from_utf8_lossy(&record[1]).parse()?;
        let end: u32 = String::from_utf8_lossy(&record[2]).parse()?;

        // Handle coordinate system
        let (adj_start, adj_end) = if args.one_based {
            (start - 1, end)
        } else {
            (start, end)
        };

        // Join remaining fields using lossy UTF-8 conversion
        let rest = if record.len() > 3 {
            record
                .iter()
                .skip(3)
                .map(|bytes| String::from_utf8_lossy(bytes))
                .collect::<Vec<_>>()
                .join("\t")
        } else {
            String::new()
        };

        // Create BedRecord
        let bed_record = BedRecord {
            start: adj_start,
            end: adj_end,
            rest,
        };

        // Add to store
        store.add_record(&chrom, &bed_record)?;

        // Update progress bar less frequently
        counter += 1;
        if counter % update_frequency == 0 {
            pb.set_position(counter);
        }
    }
    // Finalize the store
    store.finalize()?;

    pb.finish_with_message("Packing complete!");

    // If --features=dev,report how off this is
    #[cfg(feature = "dev")]
    {
        let duration = start_time.elapsed();
        let estimate_diff = (counter as f64 - initial_estimated_records as f64)
            / initial_estimated_records as f64
            * 100.0;
        eprintln!("\n--- estimate_total_records() dev stats ---");
        eprintln!("  Estimated records: {}", initial_estimated_records);
        eprintln!("  Actual records:   {}", counter);
        eprintln!("  Estimation off by: {:.1}%", estimate_diff);
        eprintln!("  Processing time:  {:?}", duration);
        eprintln!(
            "  Records/second:   {:.0}",
            counter as f64 / duration.as_secs_f64()
        );
    }

    let duration = start.elapsed();
    eprintln!("Successfully packed and indexed the file in {:?}", duration);

    Ok(())
}

pub fn build_tsv_reader(
    filepath: impl Into<PathBuf>,
    comment_char: Option<u8>,
    flexible: bool,
    has_headers: bool,
) -> Result<csv::Reader<Box<dyn std::io::Read>>, Box<dyn std::error::Error>> {
    let filepath = filepath.into();
    let input_stream = InputStream::new(&filepath);

    let stream = input_stream.reader()?;

    // Box the BufReader directly as a Read trait object
    let boxed_reader: Box<dyn Read> = Box::new(stream);

    let csv_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(has_headers)
        .comment(comment_char)
        .flexible(flexible)
        .from_reader(boxed_reader);

    Ok(csv_reader)
}

pub fn estimate_total_records(
    path: &std::path::Path,
    comment_char: Option<u8>,
    delimiter: u8,
    has_headers: bool,
    flexible: bool,
) -> Result<u64, HgIndexError> {
    let input_stream = InputStream::new(path);
    let file = File::open(path)?;

    // For small files, just count exact records
    if file.metadata()?.len() < 1024 * 1024 {
        // 1MB threshold
        let mut csv_reader = build_tsv_reader(path, comment_char, flexible, has_headers)?;
        return Ok(csv_reader.byte_records().count() as u64);
    }

    // For larger files, sample decompressed content
    let mut reader = input_stream.buffered_reader()?;
    let mut valid_lines = 0;
    let mut total_bytes = 0;
    let mut line_buffer = Vec::new();
    let sample_size = 1000; // Increased sample size for better estimation

    // Sample more lines for better estimation
    while valid_lines < sample_size {
        line_buffer.clear();
        let bytes_read = reader.read_until(b'\n', &mut line_buffer)?;
        if bytes_read == 0 {
            break; // EOF
        }

        // Skip comment lines
        if let Some(comment) = comment_char {
            if !line_buffer.is_empty() && line_buffer[0] == comment {
                continue;
            }
        }

        // Verify minimum fields
        let field_count = line_buffer.iter().filter(|&&b| b == delimiter).count() + 1;
        if !flexible && field_count < 3 {
            continue;
        }

        total_bytes += bytes_read;
        valid_lines += 1;
    }

    if valid_lines == 0 {
        return Ok(0);
    }

    // Calculate average bytes per valid line
    let avg_bytes_per_line = total_bytes as f64 / valid_lines as f64;

    // Reset reader and count total decompressed bytes
    let mut total_decompressed_size = 0;
    let mut count_buffer = [0u8; 8192];
    let mut reader = input_stream.buffered_reader()?;

    while let Ok(n) = reader.read(&mut count_buffer) {
        if n == 0 {
            break;
        }
        total_decompressed_size += n;
    }

    // Estimate based on decompressed size
    let estimated_records = (total_decompressed_size as f64 / avg_bytes_per_line) as u64;

    // Add a buffer
    let buffer = 1.05;
    Ok((estimated_records as f64 * buffer) as u64)
}
