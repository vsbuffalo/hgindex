// src/bin/index-stats.rs

use clap::Parser;
use hgindex::index::BinningIndex;
use std::path::PathBuf;

#[derive(Parser)]
#[command(author, version, about = "Analyze genomic index statistics")]
struct Args {
    /// Path to the index file
    #[arg(short, long)]
    index: PathBuf,

    /// Show detailed bin statistics per chromosome
    #[arg(short, long)]
    detailed: bool,
}

#[derive(Clone, Default)]
struct LevelStats {
    bin_count: usize,
    feature_count: usize,
    max_features: usize,
}

fn analyze_index(index: &BinningIndex<()>, detailed: bool) {
    let binning = hgindex::index::HierarchicalBins::default();
    let mut stats_by_level: Vec<LevelStats> = vec![LevelStats::default(); binning.num_levels];
    let mut total_features = 0;

    // Calculate statistics for each chromosome
    for (chrom, seq_index) in &index.sequences {
        let mut chrom_stats_by_level: Vec<LevelStats> =
            vec![LevelStats::default(); binning.num_levels];

        for (&bin_id, features) in &seq_index.bins {
            // Find which level this bin belongs to
            let level = binning
                .bin_offsets
                .iter()
                .position(|&offset| bin_id >= offset)
                .unwrap_or(0);

            // Update statistics
            let stats = &mut stats_by_level[level];
            stats.bin_count += 1;
            stats.feature_count += features.len();
            stats.max_features = stats.max_features.max(features.len());

            if detailed {
                let chrom_stats = &mut chrom_stats_by_level[level];
                chrom_stats.bin_count += 1;
                chrom_stats.feature_count += features.len();
                chrom_stats.max_features = chrom_stats.max_features.max(features.len());
            }
        }

        total_features += seq_index.bins.values().map(|v| v.len()).sum::<usize>();

        // Print chromosome-specific statistics if detailed mode is enabled
        if detailed {
            println!("\nChromosome: {}", chrom);
            println!("Level |  Bins  | Features | Avg Features/Bin | Max Features");
            println!("------|-------|----------|-----------------|-------------");
            for (level, stats) in chrom_stats_by_level.iter().enumerate() {
                if stats.bin_count > 0 {
                    println!(
                        "{:5} | {:5} | {:8} | {:15.2} | {:12}",
                        level,
                        stats.bin_count,
                        stats.feature_count,
                        stats.feature_count as f64 / stats.bin_count as f64,
                        stats.max_features
                    );
                }
            }
            println!();
        }
    }

    // Print overall statistics
    println!("\nOverall Statistics:");
    println!("Total features: {}", total_features);
    println!("\nBin Level Statistics:");
    println!("Level | Bins | Features | Avg Features/Bin | Max Features");
    println!("------|-------|----------|-----------------|-------------");
    for (level, stats) in stats_by_level.iter().enumerate() {
        if stats.bin_count > 0 {
            println!(
                "{:5} | {:5} | {:8} | {:15.2} | {:12}",
                level,
                stats.bin_count,
                stats.feature_count,
                stats.feature_count as f64 / stats.bin_count as f64,
                stats.max_features
            );
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // Read the index file
    let index = BinningIndex::open(&args.index)?;

    // Analyze and print statistics
    analyze_index(&index, args.detailed);

    Ok(())
}
