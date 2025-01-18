use clap::Args;
use hgindex::error::HgIndexError;
use hgindex::index::BinningIndex;
use hgindex::stats::BinningStats;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Args)]
pub struct StatsArgs {
    /// Input index file to analyze (should be a .hgidx file)
    #[arg(value_name = "FILE")]
    pub input: PathBuf,

    /// Print bin indices for debugging purposes
    #[arg(long)]
    pub show_bins: bool,
}

pub fn run(args: StatsArgs) -> Result<(), HgIndexError> {
    let start = Instant::now();

    // Load the BinningIndex from the input file
    eprintln!("Loading index from {}...", args.input.display());
    let index = BinningIndex::open(&args.input)?;
    eprintln!("Index loaded successfully.");

    // Compute statistics
    eprintln!("Analyzing index structure and performance...");
    let stats = BinningStats::analyze(&index);

    // Print statistics summary
    eprintln!("\nIndex Analysis Summary:");
    stats.print_summary();

    // Print detailed performance report
    let report = stats.generate_performance_report();
    println!("{}", report);

    // Optionally print bin indices
    if args.show_bins {
        eprintln!("\nBin Indices:");
        for (chrom, sequence_index) in &index.sequences {
            println!("Chromosome: {}", chrom);
            let bins: Vec<_> = sequence_index.bins.keys().cloned().collect(); // Collect bin IDs
            println!("  Bins: {:?}", bins);
        }
    }

    let duration = start.elapsed();
    eprintln!("Analysis completed in {:?}", duration);

    Ok(())
}
