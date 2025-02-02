#[cfg(feature = "cli")]
use clap::ValueEnum;
use serde::{Deserialize, Serialize};
use std::fmt;

/// index/binning.rs
///
/// # Hierarchical Binning
///
/// We can think about indexing a set of *reference ranges* (e.g. genomic features) such that one
/// can more efficiently look for reference range(s) that overlap a particular *query range* (or
/// set of query ranges).
///
/// We index the reference ranges using hierarchical bins, inspired by the schemed used by the UCSC
/// Genome Browser and tabix.
///
/// ## Hierarchical Binning and Indexing Strategy
///
/// A hierarchical binning strategy takes a reference range and places it into a unique *bin
/// level*, based on its width. The $b$ bin levels width exponents l_0, l_1, ..., l_b are a
/// geometrically decreasing sequence that create the bin widths for each level. For the UCSC
/// system, these levels' widths start at 512Mb (2^l_0 = 2^29) for level 0, and decrease by a
/// factor of 8 (2^3) at each level.
///
/// Each reference range [start, end) is placed in (or assigned) *the smallest bin it fully fits
/// into*. For looking for all reference ranges that overlap a query range, the highest level's
/// (the smallest bin) bins are searched first. These bins can be found by taking the query range's
/// start position and integer dividing it by that level's bin width (equivalent to a right bit
/// shift with the bin level's width exponent$). So, for example,
///
///   >>> 100_000_000 >> 17   # start range right shifted with level 0 (128kb bins)
/// > > > 762                     # this is the first bin index
///
///   >>> 100_191_121 >> 17   # end range right shifted with level 0
/// > > > 764                     # this is the last bin index.
///
/// The query range needs to be checked against *all levels*. If a query range is larger, it will
/// need to be checked against many of the smaller windows. This increases *total bin access* time,
/// but the savings are all the other many small bins that don't need to be checked for overlaps.
/// Note that feature density should impact things considerably, as that means more ranges per
/// small windows that need to be checked.
///
/// ## The Offset Scheme
///
/// Geometrically scaling bin widths have the advantage that indices nest neatly in each other.
/// Each level's offset is calculated to start *after* all bins in each previous level. The first
/// level has one bin, it's index is 0; the second level has eight bins, their indices are [0, 8],
/// etc
///
///  Level 0: 1 bin
///  Level 1: needs to start after bin 0     -> starts at 1
///  Level 2: needs to start after bin 8     -> starts at 9
///  Level 3: needs to start after bin 72    -> starts at 73
///  Level 4: needs to start after bin 584   -> starts at 585
///
///
///
/// ## Tuning
/// Tuning Strategy:
///
/// 1. Adjust `base_shift`:
///    - Decrease for dense small-range queries (e.g., `base_shift = 13`).
///    - Increase for large-range queries (e.g., `base_shift = 15`).
/// 2. Modify `level_shift`:
///    - Decrease for finer level granularity (e.g., `level_shift = 2`).
///    - Increase for faster bin traversal (e.g., `level_shift = 4`).
/// 3. Experiment with `num_levels`:
///    - Fewer levels (e.g., 6): Faster traversal, better for dense data.
///    - More levels (e.g., 10): Improved flexibility for mixed range sizes.
/// 4. Optimize `linear_shift`:
///    - Decrease for clustered dense queries (e.g., `linear_shift = 6`).
///    - Increase for sparse datasets (e.g., `linear_shift = 10`).

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct HierarchicalBins {
    pub schema: BinningSchema,
    /// For hierarchical bins
    pub base_shift: u32,
    /// For level scaling
    pub level_shift: u32,
    /// Number of bin levels
    pub num_levels: usize,
    /// Number of bins at each level
    pub levels: Vec<u32>,
    /// Offsets for each level
    pub bin_offsets: Vec<u32>,
    /// For optional linear index
    pub linear_shift: Option<u32>,
}

impl Default for HierarchicalBins {
    fn default() -> Self {
        Self::from_schema(&BinningSchema::default())
    }
}

/// Calculate number of bins at each level.
pub fn calc_level_sizes(next_shift: u32, nlevels: usize) -> Vec<u32> {
    let mut level_bins = Vec::with_capacity(nlevels);
    for i in 0..nlevels {
        let level = 1u32 << (next_shift * i as u32);
        level_bins.push(level);
    }
    level_bins
}

/// Calculate the bin offsets for these levels.
pub fn calc_offsets_from_levels(levels: &[u32]) -> Vec<u32> {
    // Create cumulative sum
    let mut offsets: Vec<u32> = levels
        .iter()
        .scan(0u32, |sum, &x| {
            let current = *sum;
            *sum += x;
            Some(current)
        })
        .collect();
    offsets.reverse();
    offsets
}

/// Calculate the bin offsets, from level 0 to level nlevels, shifting next_shift
/// each level.
#[allow(dead_code)]
pub fn calc_offsets(next_shift: u32, nlevels: usize) -> Vec<u32> {
    let level_bins = calc_level_sizes(next_shift, nlevels);
    calc_offsets_from_levels(&level_bins)
}

#[derive(Debug, Default, Deserialize, Serialize, PartialEq, Clone)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum BinningSchema {
    #[default]
    Tabix,
    TabixNoLinear,
    Ucsc,
    UcscNoLinear,
    Dense,
    Sparse,
}

impl fmt::Display for BinningSchema {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BinningSchema::Tabix => write!(f, "Tabix"),
            BinningSchema::TabixNoLinear => write!(f, "Tabix (No Linear Index)"),
            BinningSchema::Ucsc => write!(f, "UCSC"),
            BinningSchema::UcscNoLinear => write!(f, "UCSC (No Linear Index)"),
            BinningSchema::Dense => write!(f, "Dense"),
            BinningSchema::Sparse => write!(f, "Sparse"),
        }
    }
}

impl HierarchicalBins {
    pub fn from_schema(schema: &BinningSchema) -> Self {
        match schema {
            BinningSchema::Tabix => Self::tabix(),
            BinningSchema::TabixNoLinear => Self::tabix_no_linear(),
            BinningSchema::Ucsc => Self::ucsc(),
            BinningSchema::UcscNoLinear => Self::ucsc_no_linear(),
            BinningSchema::Dense => Self::dense(),
            BinningSchema::Sparse => Self::sparse(),
        }
    }

    pub fn new(
        schema_type: BinningSchema,
        base_shift: u32,
        level_shift: u32,
        num_levels: usize,
        linear_shift: Option<u32>,
    ) -> Self {
        let levels = calc_level_sizes(level_shift, num_levels);
        let bin_offsets = calc_offsets_from_levels(&levels);
        // TODO: check
        if base_shift + (num_levels as u32 - 1) * level_shift > 63 {
            panic!("Shift exceeds maximum allowable bit width!");
        }

        Self {
            schema: schema_type,
            base_shift,
            level_shift,
            num_levels,
            levels,
            bin_offsets,
            linear_shift,
        }
    }

    pub fn tabix() -> Self {
        Self::new(BinningSchema::Tabix, 14, 3, 6, Some(14))
    }

    pub fn tabix_no_linear() -> Self {
        Self::new(BinningSchema::TabixNoLinear, 14, 3, 6, None)
    }

    pub fn ucsc() -> Self {
        Self::new(BinningSchema::Ucsc, 17, 3, 5, Some(14))
    }

    pub fn ucsc_no_linear() -> Self {
        Self::new(BinningSchema::UcscNoLinear, 17, 3, 5, None)
    }

    pub fn dense() -> Self {
        Self::new(BinningSchema::Dense, 14, 3, 10, Some(8))
    }

    pub fn sparse() -> Self {
        Self::new(BinningSchema::Sparse, 20, 4, 4, Some(16))
    }

    pub fn uses_linear_index(&self) -> bool {
        self.linear_shift.is_some()
    }

    /// Compute the smallest bin fully containing the range `[start, end)`.
    pub fn region_to_bin(&self, start: u32, end: u32) -> u32 {
        let mut start_bin = start >> self.base_shift;
        let mut end_bin = (end - 1) >> self.base_shift;

        for &offset in &self.bin_offsets {
            if start_bin == end_bin {
                return offset + start_bin;
            }
            start_bin >>= self.level_shift;
            end_bin >>= self.level_shift;
        }

        panic!(
            "start {}, end {} out of range for region_to_bin",
            start, end
        );
    }

    /// Compute all bins potentially overlapping the range `[start, end)`.
    pub fn region_to_bins(&self, start: u32, end: u32) -> Vec<u32> {
        let mut bins = Vec::new();
        let mut start_bin = start >> self.base_shift;
        let mut end_bin = (end - 1) >> self.base_shift;

        for &offset in &self.bin_offsets {
            bins.extend(offset + start_bin..=offset + end_bin);
            start_bin >>= self.level_shift;
            end_bin >>= self.level_shift;
        }

        bins
    }

    pub fn region_to_bins_iter(&self, start: u32, end: u32) -> RegionToBins {
        let start_bin = start >> self.base_shift;
        let end_bin = (end - 1) >> self.base_shift;

        RegionToBins {
            current_level: 0,
            start_bin,
            end_bin,
            bin_offsets: &self.bin_offsets,
            level_shift: self.level_shift,
        }
    }
}

pub struct RegionToBins<'a> {
    current_level: usize,
    start_bin: u32,
    end_bin: u32,
    bin_offsets: &'a [u32],
    level_shift: u32,
}

impl Iterator for RegionToBins<'_> {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_level >= self.bin_offsets.len() {
            return None; // No more levels to iterate
        }

        // Return the current bin
        let current_bin = self.bin_offsets[self.current_level] + self.start_bin;

        if self.start_bin < self.end_bin {
            // Move to the next bin within the current level
            self.start_bin += 1;
        } else {
            // Move to the next level
            self.current_level += 1;

            if self.current_level < self.bin_offsets.len() {
                // Reset bins for the next level
                self.start_bin >>= self.level_shift;
                self.end_bin >>= self.level_shift;
            }
        }

        Some(current_bin)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use std::collections::HashSet;

    #[test]
    fn test_ucsc_calc_levels() {
        // Calculate bin counts for each level
        let bin_counts = vec![1, 8, 64, 512, 4096];
        assert_eq!(bin_counts, calc_level_sizes(3, 5));
    }

    #[test]
    fn test_ucsc_offsets() {
        let offsets = calc_offsets(3, 5);
        assert_eq!(offsets, vec![585, 73, 9, 1, 0]);
    }

    #[test]
    fn test_ucsc_extended_offsets() {
        // For the extended case that handles up to 2Gb
        let offsets = calc_offsets(3, 6);
        // These values come from binOffsetsExtended in the UCSC code
        assert_eq!(offsets, vec![4681, 585, 73, 9, 1, 0]);
    }

    #[test]
    fn test_bin_counts_per_level() {
        let offsets = calc_offsets(3, 5);

        // Verify the differences between offsets match UCSC bin counts
        // Level 0: 1 bin    (512Mb)
        // Level 1: 8 bins   (64Mb bins)
        // Level 2: 64 bins  (8Mb bins)
        // Level 3: 512 bins (1Mb bins)
        // Level 4: 4096 bins (128kb bins)

        assert_eq!(offsets[0] - offsets[1], 512); // Level 0 offset difference
        assert_eq!(offsets[1] - offsets[2], 64); // Level 1 offset difference
        assert_eq!(offsets[2] - offsets[3], 8); // Level 2 offset difference
        assert_eq!(offsets[3] - offsets[4], 1); // Level 3 offset difference
        assert_eq!(offsets[4], 0); // Level 4 starts at 0
    }

    #[test]
    fn test_non_overlapping_ranges() {
        let offsets = calc_offsets(3, 5);
        assert_eq!(offsets, vec![585, 73, 9, 1, 0]);

        // Calculate bin counts for each level
        let bin_counts = vec![1, 8, 64, 512, 4096];
        assert_eq!(bin_counts, calc_level_sizes(3, 5));

        // For each level, verify that:
        // 1. The offset is the start of that level's range
        // 2. The offset plus bin_count-1 gives the end of the range
        // 3. This matches UCSC's scheme
        assert_eq!(offsets[0], 585); // Level 0 starts at 585
        assert_eq!(offsets[0] + bin_counts[0] - 1, 585); // 1 bin

        assert_eq!(offsets[1], 73); // Level 1 starts at 73
        assert_eq!(offsets[1] + bin_counts[1] - 1, 80); // 8 bins

        assert_eq!(offsets[2], 9); // Level 2 starts at 9
        assert_eq!(offsets[2] + bin_counts[2] - 1, 72); // 64 bins

        assert_eq!(offsets[3], 1); // Level 3 starts at 1
        assert_eq!(offsets[3] + bin_counts[3] - 1, 512); // 512 bins

        assert_eq!(offsets[4], 0); // Level 4 starts at 0
        assert_eq!(offsets[4] + bin_counts[4] - 1, 4095); // 4096 bins
    }

    #[test]
    fn test_region_to_bin() {
        let index = HierarchicalBins::ucsc();

        // Test cases from UCSC example in documentation:
        // "100_000_000 >> 17" gives bin 762
        assert_eq!(index.region_to_bin(100_000_000, 100_000_100), 762 + 585);

        // Test different size ranges that should go into different bin levels

        // Small range (fits in level 4 - 128kb bins)
        assert_eq!(index.region_to_bin(0, 1000), 585); // Should be first bin at finest level

        // 1MB range (should go to level 3)
        assert_eq!(index.region_to_bin(1_000_000, 2_000_000), 9); // Level 2 offset + bin 0

        // 10MB range (should go to level 1 - 64MB bins)
        assert_eq!(index.region_to_bin(10_000_000, 20_000_000), 1); // Level 1 offset + bin 0

        // 100MB range (goes to level 0 - 512MB bins)
        assert_eq!(index.region_to_bin(100_000_000, 200_000_000), 0); // Level 0 offset + bin 0

        // 500MB range (should go to level 0)
        assert_eq!(index.region_to_bin(0, 500_000_000), 0); // Level 0 offset + bin number

        // Test edge cases

        #[allow(non_upper_case_globals)]
        const KiB: u32 = 1024;

        // Test exact bin boundaries
        assert_eq!(index.region_to_bin(0, 128 * KiB), 585); // Exactly one 128kb bin
        assert_eq!(index.region_to_bin(128 * KiB, 256 * KiB), 586); // Second 128kb bin

        // Test adjacent regions get different bins
        let bin1 = index.region_to_bin(0, 128_000);
        let bin2 = index.region_to_bin(128_000, 256_000);
        assert_ne!(bin1, bin2);
    }

    #[test]
    fn test_region_to_bins() {
        // Test with each schema type
        for schema in [
            BinningSchema::Tabix,
            BinningSchema::Dense,
            BinningSchema::Sparse,
            BinningSchema::Ucsc,
        ] {
            let index = HierarchicalBins::from_schema(&schema);

            // Test small range within a single smallest bin
            let bins = index.region_to_bins(1000, 2000);
            let smallest_level_offset = index.bin_offsets[0];
            assert!(!bins.is_empty());
            assert!(bins.contains(&smallest_level_offset)); // Should contain a bin from smallest level

            // Test range spanning multiple smallest bins
            let base_shift = index.base_shift;
            let bins = index.region_to_bins(1 << base_shift, 2 << base_shift); // Span exactly 2 bins
            assert!(bins.len() >= 2); // Should contain at least 2 bins

            // Test larger range that requires checking multiple levels
            let bins = index.region_to_bins(0, 1_000_000);
            // Should have bins from multiple levels
            let level_3_start = index.bin_offsets[index.bin_offsets.len() - 2];
            let level_4_start = index.bin_offsets[index.bin_offsets.len() - 1];
            assert!(bins.iter().any(|&x| x >= level_3_start)); // Should have higher level bins
            assert!(bins.iter().any(|&x| x >= level_4_start)); // Should have lower level bins

            // Test exact bin boundary
            let bins = index.region_to_bins(0, 1 << base_shift); // Exactly one smallest bin
            assert!(bins.contains(&smallest_level_offset)); // Should contain first bin

            // Test that bins are unique
            let bins = index.region_to_bins(0, 10_000_000);
            let unique_bins: HashSet<_> = bins.iter().cloned().collect();
            assert_eq!(bins.len(), unique_bins.len());
        }
    }

    fn test_with_all_configs<F>(test_fn: F)
    where
        F: Fn(&HierarchicalBins),
    {
        for config in [
            HierarchicalBins::default(),
            HierarchicalBins::ucsc(),
            HierarchicalBins::dense(),
            HierarchicalBins::sparse(),
        ] {
            test_fn(&config);
        }
    }

    #[test]
    fn test_bin_boundaries_all_configs() {
        test_with_all_configs(|index| {
            let bin_size = 1 << index.base_shift;
            let bin1 = index.region_to_bin(0, bin_size);
            let bin2 = index.region_to_bin(bin_size, 2 * bin_size);
            assert_ne!(bin1, bin2);
        });
    }

    proptest! {
        #[test]
        fn test_region_to_bins_properties(start in 0u32..1_000_000, len in 1u32..1_000_000) {
            test_with_all_configs(|index| {
                let end = start.saturating_add(len);
                let bins = index.region_to_bins(start, end);

                // Properties that should hold for all configs:
                assert!(!bins.is_empty()); // Should always return some bins

                // Bins should be unique
                let unique_bins: HashSet<_> = bins.iter().collect();
                assert_eq!(bins.len(), unique_bins.len());

                // Every bin should be within valid ranges for the scheme
                let max_valid_bin = index.bin_offsets[0] + index.levels[index.num_levels - 1];
                assert!(bins.iter().all(|&bin| bin < max_valid_bin),
                    "Found bin larger than max valid bin {} in {:?}", max_valid_bin, bins);
            });
        }
    }
}
