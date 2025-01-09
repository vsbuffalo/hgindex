// bin/commands/query.rs

#[cfg(feature = "cli")]
pub mod query {
    use crate::commands::BedRecord;
    use clap::Args;
    use flate2::Compression;
    use hgindex::error::HgIndexError;
    use hgindex::io::io::OutputStream;
    use hgindex::store::GenomicDataStore;
    use std::fs;
    use std::path::PathBuf;
    use std::time::Instant;

    #[derive(Args)]
    pub struct QueryArgs {
        /// Output file.
        #[arg(short, long, value_name = "overlaps.bed")]
        pub output: Option<String>,

        /// The query region, in the format seqname:start-end where start and end are
        /// 1-based inclusive coordinates (like tabix's region argument).
        #[arg(value_name = "chr17:7661779-7687538")]
        pub region: String,

        /// Input .hgidx directory. If not specified, a file with the suffix .hgidx
        /// will be looked for in the current directory. If a single match is found,
        /// it will be used.
        #[arg(value_name = "scores.hgidx")]
        pub input: Option<PathBuf>,
    }

    pub fn run(args: QueryArgs) -> Result<(), HgIndexError> {
        // For timing the query operation
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

        eprintln!("Query region {} in {}", args.region, input_path.display());

        // Now, use the input_path to create a GenomicDataStore
        let mut store = GenomicDataStore::<BedRecord, ()>::open(&input_path, None)?;

        // Parse the region into the appropriate format
        let region_parts: Vec<&str> = args.region.split(':').collect();
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

        // Query the store for the region
        let records = store.get_overlapping(seqname, start, end)?;
        eprintln!("{} records found.", records.len());

        // Output the records
        // TODO - add more features?
        for record in records {
            writeln!(output_writer, "{}", record)?;
        }

        let duration = duration_start.elapsed();
        eprintln!("Query completed in {:?}", duration);

        Ok(())
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
}
