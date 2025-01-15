// bin/commands/pack.rs

use clap::Args;
use csv::ReaderBuilder;
use flate2::read::GzDecoder;
use glob::glob;
use hgindex::error::HgIndexError;
use hgindex::store::GenomicDataStore;
use hgindex::BedRecord;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::time::Instant;

use csv::Reader as CsvReader;

#[derive(Args)]
pub struct PackArgs {
    /// Input TSV/BED file to pack and index (conflicts with --glob)
    #[arg(value_name = "FILE", conflicts_with = "glob_pattern")]
    pub input: Option<PathBuf>,

    /// Glob pattern for chromosome-specific files (e.g. "data/{}.bed")
    #[arg(long, value_name = "PATTERN", conflicts_with = "input")]
    pub glob_pattern: Option<String>,

    /// Optional list of chromosomes to process. If not specified, will use all files matching pattern
    #[arg(long, value_name = "CHROM")]
    pub chromosomes: Option<Vec<String>>,

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
    #[arg(long, value_enum, default_value_t = hgindex::BinningSchema::Tabix)]
    pub schema: hgindex::BinningSchema,

    /// Number of threads to use for parallel processing (default: 1)
    #[arg(long, default_value_t = 1)]
    pub threads: usize,
}

impl PackArgs {
    pub fn validate(&self) -> Result<(), HgIndexError> {
        match (&self.input, &self.glob_pattern) {
            (None, None) => Err("Either --input or --glob-pattern must be specified".into()),
            (Some(_), Some(_)) => Err("Cannot specify both --input and --glob-pattern".into()),
            _ => Ok(()),
        }
    }

    pub fn get_input_files(&self) -> Result<Vec<PathBuf>, HgIndexError> {
        match (&self.input, &self.glob_pattern) {
            (Some(input), None) => Ok(vec![input.clone()]),
            (None, Some(pattern)) => {
                // Collect all paths or propagate errors
                let files: Vec<PathBuf> = glob(&pattern)
                    .map_err(|e| HgIndexError::StringError(format!("Invalid glob pattern: {}", e)))?
                    .collect::<Result<Vec<PathBuf>, glob::GlobError>>()?; // Convert to a Vec<PathBuf>

                if files.is_empty() {
                    return Err("No files found matching pattern".into());
                }

                Ok(files) // Now `files` is guaranteed to be a Vec<PathBuf>
            }
            _ => unreachable!("Should have been caught by validate()"),
        }
    }
}

pub fn run(args: PackArgs) -> Result<(), HgIndexError> {
    // For timing the pack operation
    let start = Instant::now();

    // Validate arguments
    args.validate()?;

    // Configure the number of threads for parallel processing
    ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .map_err(|e| HgIndexError::StringError(format!("Failed to set thread pool: {}", e)))?;

    // Get input files
    let input_files = args.get_input_files()?;

    // Determine output path
    let output_path = match args.output {
        Some(path) => path,
        None => {
            if args.glob_pattern.is_some() {
                return Err("Output path must be specified when using --glob".into());
            }
            input_files[0].with_extension("hgidx")
        }
    };

    // Create the genomic data store
    let mut store =
        GenomicDataStore::<BedRecord>::create_with_schema(&output_path, None, &args.schema)?;

    // Conditional progress bar handling based on thread count
    if args.threads == 1 {
        // Single-threaded: Use a simple progress bar
        let overall_pb = ProgressBar::new(input_files.len() as u64);
        overall_pb.set_style(
            ProgressStyle::default_bar()
                .template(
                    "[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({percent}%) Processing",
                )
                .expect("Failed to create progress bar style")
                .progress_chars("=>-"),
        );

        for input_file in input_files {
            overall_pb.set_message(format!("Processing {}", input_file.display()));

            // Estimate total records
            let estimated_records =
                estimate_total_records(&input_file, Some(args.comment as u8), b'\t', false, true)?;

            // Process the file
            let mut csv_reader =
                build_tsv_reader(&input_file, Some(args.comment as u8), true, false)?;
            for result in csv_reader.records() {
                let record = result?;

                if record.len() < 3 {
                    eprintln!(
                        "Warning: skipping insufficient fields in file {}",
                        input_file.display()
                    );
                    continue;
                }

                // Parse and process the record
                let chrom = record.get(0).unwrap().to_string();
                let start: u32 = record.get(1).unwrap().parse()?;
                let end: u32 = record.get(2).unwrap().parse()?;

                let (adj_start, adj_end) = if args.one_based {
                    (start - 1, end)
                } else {
                    (start, end)
                };

                let rest = if record.len() > 3 {
                    record.iter().skip(3).collect::<Vec<_>>().join("\t")
                } else {
                    String::new()
                };

                let bed_record = BedRecord {
                    chrom,
                    start: adj_start,
                    end: adj_end,
                    rest,
                };

                store.add_record(&bed_record.chrom, &bed_record)?;
            }

            overall_pb.inc(1);
        }

        overall_pb.finish_with_message("All files processed successfully!");
    } else {
        // Multi-threaded: Use MultiProgress
        let multi_pb = MultiProgress::new();
        let overall_pb = multi_pb.add(ProgressBar::new(input_files.len() as u64));
        overall_pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({percent}%) Overall Progress")
                .expect("Failed to create overall progress bar style")
                .progress_chars("=>-"),
        );

        // Process files in parallel
        input_files
            .par_iter()
            .try_for_each(|input_file| {
                // Estimate total records
                let estimated_records =
                    estimate_total_records(input_file, Some(args.comment as u8), b'\t', false, true)?;

                // Add a progress bar for the current file
                let file_pb = multi_pb.add(ProgressBar::new(estimated_records));
                file_pb.set_style(
                    ProgressStyle::default_bar()
                        .template("{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue} ⟩ {pos}/{len} ({percent}%) {msg}")
                        .expect("Failed to create file progress bar style")
                        .progress_chars("=>-"),
                );
                file_pb.set_message(format!("Processing {}", input_file.display()));

                let mut csv_reader = build_tsv_reader(input_file, Some(args.comment as u8), true, false)?;
                for result in csv_reader.records() {
                    let record = result?;
                    // Process the record...
                    file_pb.inc(1);
                }

                file_pb.finish_with_message(format!("Finished {}", input_file.display()));
                overall_pb.inc(1);

                Ok::<(), HgIndexError>(())
            })?;

        overall_pb.finish_with_message("All files processed successfully!");
    }

    // Finalize the genomic data store
    store.finalize()?;

    let duration = start.elapsed();
    eprintln!("Successfully packed and indexed the file in {:?}", duration);

    Ok(())
}

pub fn build_tsv_reader(
    filepath: impl Into<PathBuf>,
    comment_char: Option<u8>,
    flexible: bool,
    has_headers: bool,
) -> Result<CsvReader<Box<dyn Read>>, HgIndexError> {
    let filepath = filepath.into();
    let file = File::open(&filepath)?;
    let is_gzipped = is_gzipped_file(&filepath)?;
    let stream: Box<dyn Read> = if is_gzipped {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };

    let reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(has_headers)
        .comment(comment_char)
        .flexible(flexible)
        .from_reader(stream);

    Ok(reader)
}

// Helper function to check if file is gzipped
fn is_gzipped_file(filepath: &Path) -> Result<bool, HgIndexError> {
    let mut file = File::open(filepath)?;
    let mut magic = [0u8; 2];
    if file.read_exact(&mut magic).is_ok() {
        return Ok(magic[0] == 0x1f && magic[1] == 0x8b);
    }

    Ok(false)
}

pub fn estimate_total_records(
    path: &std::path::Path,
    comment_char: Option<u8>,
    delimiter: u8,
    has_headers: bool,
    flexible: bool,
) -> Result<u64, HgIndexError> {
    let file = File::open(path)?;
    let file_size = file.metadata()?.len();

    // If file is small, just count exact lines
    if file_size < 1024 * 1024 {
        // 1MB threshold
        let mut csv_reader = build_tsv_reader(path, comment_char, flexible, has_headers)?;

        return Ok(csv_reader.records().count() as u64);
    }

    // Sample the first 100 valid records to get average line size
    let mut reader = BufReader::new(file);
    let mut valid_lines = 0;
    let mut total_bytes = 0;
    let mut line = String::new();
    let sample_size = 100;

    while valid_lines < sample_size {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            break; // EOF
        }

        // Skip comment lines
        if let Some(comment) = comment_char {
            if line.starts_with(comment as char) {
                continue;
            }
        }

        // Check if line has minimum required fields
        let fields: Vec<&str> = line.trim().split(delimiter as char).collect();
        if !flexible && fields.len() < 3 {
            continue;
        }

        total_bytes += bytes_read;
        valid_lines += 1;
    }

    if valid_lines == 0 {
        return Ok(0);
    }

    // Calculate average bytes per valid line
    let avg_line_size = total_bytes as f64 / valid_lines as f64;

    // Estimate total records
    let estimated_records = (file_size as f64 / avg_line_size) as u64;

    // Add a small buffer for estimation errors
    Ok((estimated_records as f64 * 1.1) as u64)
}
