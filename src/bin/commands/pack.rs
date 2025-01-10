// bin/commands/pack.rs

#[cfg(feature = "cli")]
pub mod pack {
    use crate::commands::BedRecord;
    use clap::Args;
    use csv::ReaderBuilder;
    use flate2::read::GzDecoder;
    use hgindex::error::HgIndexError;
    use hgindex::store::GenomicDataStore;
    use indicatif::{ProgressBar, ProgressStyle};
    use std::fs::File;
    use std::io::{BufRead, BufReader, Read};
    use std::path::{Path, PathBuf};
    use std::time::Instant;

    use csv::Reader as CsvReader;

    #[derive(Args)]
    pub struct PackArgs {
        /// Input TSV/BED file to pack and index
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
        #[arg(long, value_enum, default_value_t = hgindex::BinningSchema::Tabix)]
        pub schema: hgindex::BinningSchema,
    }

    pub fn run(args: PackArgs) -> Result<(), HgIndexError> {
        // For timing the pack operation
        let start = Instant::now();

        // Get output path
        let output_path = args
            .output
            .unwrap_or_else(|| args.input.with_extension("hgidx"));

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
        let mut store = GenomicDataStore::<BedRecord, ()>::create_with_schema(
            &output_path,
            None,
            &args.schema,
        )?;

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

        // Process records
        for (line_no, result) in csv_reader.records().enumerate() {
            let record = result?;

            // Need at least 3 fields
            if record.len() < 3 {
                eprintln!(
                    "Warning: skipping line {} - insufficient fields",
                    line_no + 1
                );
                continue;
            }

            // Parse coordinates
            let chrom = record.get(0).unwrap().to_string();
            let start: u32 = record.get(1).unwrap().parse()?;
            let end: u32 = record.get(2).unwrap().parse()?;

            // Handle coordinate system
            let (adj_start, adj_end) = if args.one_based {
                (start - 1, end) // Convert 1-based to 0-based
            } else {
                (start, end)
            };

            // Join remaining fields
            let rest = if record.len() > 3 {
                record.iter().skip(3).collect::<Vec<_>>().join("\t")
            } else {
                String::new()
            };

            // Create BedRecord
            let bed_record = BedRecord {
                chrom: chrom.clone(),
                start: adj_start,
                end: adj_end,
                rest,
            };

            // Add to store
            store.add_record(&chrom, &bed_record)?;
            pb.inc(1);
        }

        // Finalize the store
        store.finalize()?;

        pb.finish_with_message("Packing complete!");

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
            let mut csv_reader = build_tsv_reader(&path, comment_char, flexible, has_headers)?;

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
}
