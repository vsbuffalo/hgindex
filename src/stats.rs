// src/stats.rs

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::BinningIndex;

/// Detailed statistics about bin utilization and access patterns
#[derive(Debug, Serialize, Deserialize)]
pub struct BinningStats {
    // Overall stats
    pub total_features: u64,
    pub total_bins_used: u32,
    pub total_possible_bins: u32,
    pub bin_utilization: f64, // percentage of bins actually used

    // Per-level analysis
    pub level_stats: Vec<LevelStats>,

    // Distribution analysis
    pub bin_occupancy: HashMap<u32, usize>, // bin_id -> feature count
    pub feature_size_dist: SizeDistribution,

    // Query performance predictors
    pub bin_density: f64,     // avg features per used bin
    pub feature_overlap: f64, // avg number of bins a feature maps to
    pub level_overhead: f64,  // avg number of higher level bins checked

    // Schema details
    pub schema_type: String,
    pub base_shift: u32,
    pub level_shift: u32,
    pub num_levels: usize,
    pub linear_index_present: bool,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct LevelStats {
    pub level: usize,
    pub bins_used: u32,
    pub total_bins: u32,
    pub features_count: u64,
    pub utilization: f64,
    pub avg_features_per_bin: f64,
    pub max_features_in_bin: usize,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct SizeDistribution {
    pub min_size: u32,
    pub max_size: u32,
    pub mean_size: f64,
    pub median_size: f64,
    pub size_histogram: HashMap<u32, u32>, // size bucket -> count
}

impl BinningStats {
    /// Analyze a BinningIndex to generate comprehensive statistics
    pub fn analyze(index: &BinningIndex) -> Self {
        let mut stats = BinningStats {
            total_features: 0,
            total_bins_used: 0,
            total_possible_bins: 0,
            bin_utilization: 0.0,
            level_stats: Vec::new(),
            bin_occupancy: HashMap::new(),
            feature_size_dist: SizeDistribution::default(),
            bin_density: 0.0,
            feature_overlap: 0.0,
            level_overhead: 0.0,
            schema_type: index.bins.schema.to_string(),
            base_shift: index.bins.base_shift,
            level_shift: index.bins.level_shift,
            num_levels: index.bins.num_levels,
            linear_index_present: index.bins.uses_linear_index(),
        };

        // Calculate total possible bins
        stats.total_possible_bins = index.bins.levels.iter().sum();

        // Collect feature statistics across all sequences
        let mut all_sizes = Vec::new();
        let mut total_bins_hit = 0u64;

        for seq_index in index.sequences.values() {
            // Count features and bin usage
            for (bin_id, features) in &seq_index.bins {
                stats.bin_occupancy.insert(*bin_id, features.len());
                stats.total_features += features.len() as u64;

                // Collect size info
                for feature in features {
                    let size = feature.end - feature.start;
                    all_sizes.push(size);

                    // Count how many bins this feature hits
                    let bins_hit = index.bins.region_to_bins(feature.start, feature.end).len();
                    total_bins_hit += bins_hit as u64;
                }
            }
        }

        // Calculate bin utilization
        stats.total_bins_used = stats.bin_occupancy.len() as u32;
        stats.bin_utilization =
            (stats.total_bins_used as f64 / stats.total_possible_bins as f64) * 100.0;

        // Calculate per-level statistics
        stats.level_stats = Self::calculate_level_stats(index);

        // Calculate feature size distribution
        if !all_sizes.is_empty() {
            all_sizes.sort_unstable();
            stats.feature_size_dist = SizeDistribution {
                min_size: *all_sizes.first().unwrap(),
                max_size: *all_sizes.last().unwrap(),
                mean_size: all_sizes.iter().sum::<u32>() as f64 / all_sizes.len() as f64,
                median_size: all_sizes[all_sizes.len() / 2] as f64,
                size_histogram: Self::calculate_size_histogram(&all_sizes),
            };
        }

        // Calculate query performance predictors
        if stats.total_features > 0 {
            stats.bin_density = stats.total_features as f64 / stats.total_bins_used as f64;
            stats.feature_overlap = total_bins_hit as f64 / stats.total_features as f64;

            // Calculate average overhead from higher level bins
            let mut total_overhead = 0.0;
            for level_stat in &stats.level_stats {
                if level_stat.level > 0 {
                    total_overhead += level_stat.avg_features_per_bin * level_stat.bins_used as f64;
                }
            }
            stats.level_overhead = total_overhead / stats.total_features as f64;
        }

        stats
    }

    fn calculate_level_stats(index: &BinningIndex) -> Vec<LevelStats> {
        let mut level_stats = Vec::new();
        let mut current_offset = 0;

        for level in 0..index.bins.num_levels {
            let level_size = index.bins.levels[level];
            let mut level_bins = HashMap::new();

            // Collect bins for this level
            for (bin_id, features) in index
                .sequences
                .values()
                .flat_map(|seq| &seq.bins)
                .filter(|(id, _)| **id >= current_offset && **id < current_offset + level_size)
            {
                level_bins.insert(*bin_id, features.len());
            }

            if !level_bins.is_empty() {
                let max_features = level_bins.values().max().copied().unwrap_or(0);
                let total_features = level_bins.values().sum::<usize>() as u64;
                let bins_used = level_bins.len() as u32;

                level_stats.push(LevelStats {
                    level,
                    bins_used,
                    total_bins: level_size,
                    features_count: total_features,
                    utilization: (bins_used as f64 / level_size as f64) * 100.0,
                    avg_features_per_bin: total_features as f64 / bins_used as f64,
                    max_features_in_bin: max_features,
                });
            }

            current_offset += level_size;
        }

        level_stats
    }

    fn calculate_size_histogram(sizes: &[u32]) -> HashMap<u32, u32> {
        let mut histogram = HashMap::new();

        // Create log-scale buckets
        for &size in sizes {
            let bucket = (size as f64).log2().floor() as u32;
            *histogram.entry(bucket).or_default() += 1;
        }

        histogram
    }

    /// Generate a detailed report analyzing why this binning schema performs as it does
    pub fn generate_performance_report(&self) -> String {
        let mut report = String::new();

        // Overview section
        report.push_str("\nBinning Schema Performance Analysis\n");
        report.push_str("================================\n\n");
        report.push_str(&format!("Schema: {}\n", self.schema_type));
        report.push_str(&format!("Base shift: {}\n", self.base_shift));
        report.push_str(&format!("Level shift: {}\n", self.level_shift));
        report.push_str(&format!("Number of levels: {}\n\n", self.num_levels));

        // Key metrics
        report.push_str("Key Performance Metrics:\n");
        report.push_str(&format!(
            "- Total features indexed: {}\n",
            self.total_features
        ));
        report.push_str(&format!(
            "- Bin utilization: {:.2}%\n",
            self.bin_utilization
        ));
        report.push_str(&format!(
            "- Average features per used bin: {:.2}\n",
            self.bin_density
        ));
        report.push_str(&format!(
            "- Average bins per feature: {:.2}\n",
            self.feature_overlap
        ));
        report.push_str(&format!(
            "- Level traversal overhead: {:.2}\n\n",
            self.level_overhead
        ));

        // Level analysis

        report.push_str("Level-by-Level Analysis:\n");
        for level in &self.level_stats {
            // TODO
            //let bin_size_kb =
            //    1 << (self.base_shift + (level.level as u32 * self.level_shift)) >> 10;
            //report.push_str(&format!(
            //    "Level {} ({}kb bins):\n",
            //    level.level, bin_size_kb,
            //));
            report.push_str(&format!("  - Utilization: {:.2}%\n", level.utilization));
            report.push_str(&format!("  - Features: {}\n", level.features_count));
            report.push_str(&format!(
                "  - Avg features/bin: {:.2}\n",
                level.avg_features_per_bin
            ));
            report.push_str(&format!(
                "  - Max features in any bin: {}\n",
                level.max_features_in_bin
            ));
        }

        // Feature size analysis
        report.push_str("\nFeature Size Distribution:\n");
        report.push_str(&format!(
            "- Min size: {}\n",
            self.feature_size_dist.min_size
        ));
        report.push_str(&format!(
            "- Max size: {}\n",
            self.feature_size_dist.max_size
        ));
        report.push_str(&format!(
            "- Mean size: {:.2}\n",
            self.feature_size_dist.mean_size
        ));
        report.push_str(&format!(
            "- Median size: {}\n",
            self.feature_size_dist.median_size
        ));

        // Performance implications
        report.push_str("\nPerformance Analysis:\n");

        // Analyze bin density impact
        if self.bin_density > 100.0 {
            report.push_str(
                "- High bin density suggests potential bottleneck in heavily populated bins\n",
            );
        } else if self.bin_density < 10.0 {
            report.push_str(
                "- Low bin density indicates efficient space usage but possible memory overhead\n",
            );
        }

        // Analyze level overhead
        if self.level_overhead > 2.0 {
            report.push_str(
                "- High level overhead suggests many features span multiple bin levels\n",
            );
        } else {
            report.push_str("- Low level overhead indicates efficient bin level utilization\n");
        }

        // Schema-specific analysis
        match self.schema_type.as_str() {
            "Dense" => {
                report.push_str("\nDense Schema Analysis:\n");
                report.push_str("- Smaller bins enable more precise feature location\n");
                report.push_str("- Higher bin count provides better distribution\n");
                if self.bin_utilization < 30.0 {
                    report
                        .push_str("- Low bin utilization suggests schema might be too granular\n");
                }
            }
            "Sparse" => {
                report.push_str("\nSparse Schema Analysis:\n");
                report.push_str("- Larger bins reduce memory overhead\n");
                report.push_str("- Fewer levels minimize traversal time\n");
                if self.bin_density > 50.0 {
                    report.push_str("- High bin density suggests potential for bin hotspots\n");
                }
            }
            schema => {
                report.push_str(&format!("\nOther Schema ({}) Analysis:\n", schema).to_string());
                report.push_str("- Balanced approach between memory and lookup speed\n");
                if self.feature_overlap > 3.0 {
                    report.push_str("- High feature overlap suggests many large features\n");
                }
            }
        }

        report
    }

    /// Print a condensed summary of the most important stats
    pub fn print_summary(&self) {
        println!("\nBinning Stats Summary");
        println!("====================");
        println!("Schema: {} ({} levels)", self.schema_type, self.num_levels);
        println!("Total features: {}", self.total_features);
        println!("Bin utilization: {:.2}%", self.bin_utilization);
        println!("Avg features/bin: {:.2}", self.bin_density);
        println!("Avg bins/feature: {:.2}", self.feature_overlap);
        println!("Level overhead: {:.2}", self.level_overhead);
    }
}

impl Default for SizeDistribution {
    fn default() -> Self {
        Self {
            min_size: 0,
            max_size: 0,
            mean_size: 0.0,
            median_size: 0.0,
            size_histogram: HashMap::new(),
        }
    }
}

