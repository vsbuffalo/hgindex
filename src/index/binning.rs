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
///   762                     # this is the first bin index
///
///   >>> 100_191_121 >> 17   # end range right shifted with level 0
///   764                     # this is the last bin index.
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

#[derive(Clone, Debug)]
pub struct HierarchicalBins {
    pub base_shift: u32,
    pub level_shift: u32,
    pub num_levels: usize,
    pub levels: Vec<u32>,
    pub bin_offsets: Vec<u32>,
}

impl Default for HierarchicalBins {
    fn default() -> Self {
        let levels = calc_level_sizes(3, 5);
        let bin_offsets = calc_offsets_from_levels(&levels);
        HierarchicalBins {
            base_shift: 17,
            level_shift: 3,
            num_levels: 5,
            levels,
            bin_offsets,
        }
    }
}

/// Calculate number of bins at each level.
pub fn calc_level_sizes(next_shift: u32, nlevels: usize) -> Vec<u32> {
    let mut level_bins = Vec::with_capacity(nlevels);
    for i in 0..nlevels {
        let level = 1u32 << (next_shift * i as u32);
        level_bins.push(level);
    }
    // dbg!(&level_bins);
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
pub fn calc_offsets(next_shift: u32, nlevels: usize) -> Vec<u32> {
    let level_bins = calc_level_sizes(next_shift, nlevels);
    calc_offsets_from_levels(&level_bins)
}

impl HierarchicalBins {
    // Create a new binning scheme with custom parameters
    pub fn new(base_shift: u32, level_shift: u32, num_levels: usize) -> Self {
        let levels = calc_level_sizes(level_shift, num_levels);
        let bin_offsets = calc_offsets_from_levels(&levels);
        Self {
            base_shift,
            level_shift,
            num_levels,
            levels,
            bin_offsets,
        }
    }

    // Create the classic UCSC configuration
    pub fn ucsc() -> Self {
        Self::new(17, 3, 5) // 128kb base bins, 8x scaling, 5 levels
    }

    // Create a denser binning scheme for higher resolution
    pub fn dense() -> Self {
        Self::new(14, 2, 6) // 16kb base bins, 4x scaling, 6 levels
    }

    // Create a sparser scheme for large regions
    pub fn sparse() -> Self {
        Self::new(20, 4, 4) // 1Mb base bins, 16x scaling, 4 levels
    }

    /// Find the smallest bin that fully contains this range,
    /// [start, end).
    pub fn region_to_bin(&self, start: u32, end: u32) -> u32 {
        let mut start_bin = start;
        let mut end_bin = end - 1; // Exclusive end, so subtract 1

        // First shift
        start_bin >>= self.base_shift;
        end_bin >>= self.base_shift;

        // Check each level
        for (i, &offset) in self.bin_offsets.iter().enumerate() {
            if start_bin == end_bin {
                dbg!(start, end, offset + start_bin);
                return offset + start_bin;
            }
            start_bin >>= self.level_shift;
            end_bin >>= self.level_shift;
        }

        panic!("start {}, end {} out of range in region_to_bin", start, end);
    }
    /// Find all bins that could contain features overlapping this range [start, end).
    /// This is used when searching for overlaps, as we need to check multiple bins
    /// across different levels.
    pub fn region_to_bins(&self, start: u32, end: u32) -> Vec<u32> {
        let mut bins = Vec::new();
        let mut start_bin = start;
        let mut end_bin = end - 1; // Exclusive end, so subtract 1

        // First shift with the initial shift amount (17 by default)
        start_bin >>= self.base_shift;
        end_bin >>= self.base_shift;

        // For each level
        for (level, &offset) in self.bin_offsets.iter().enumerate() {
            // Calculate all bins for this level between start and end
            let level_start = start_bin;
            let level_end = end_bin;

            // Add all bins in this level's range
            for bin in level_start..=level_end {
                bins.push(offset + bin);
            }

            // Early exit if we've reached a level where start and end are in the same bin
            if start_bin == end_bin {
                break;
            }

            // Shift to next level (dividing by 8, equivalent to right shift by 3)
            start_bin >>= self.level_shift;
            end_bin >>= self.level_shift;
        }

        bins
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ucsc_calc_levels() {
        // Calculate bin counts for each level
        let bin_counts = vec![1, 8, 64, 512, 4096];
        assert_eq!(bin_counts, calc_level_sizes(3, 5));
    }

    #[test]
    fn test_ucsc_offsets() {
        let offsets = calc_offsets(3, 5);
        dbg!(&offsets);
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
        let index = HierarchicalBins::default();

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
        let index = HierarchicalBins::default();

        // Test small range within a single smallest bin
        let bins = index.region_to_bins(1000, 2000);
        assert!(bins.contains(&(585 + 0))); // Should contain the level 4 bin

        // Test range spanning multiple smallest bins
        //let bins = index.region_to_bins(128_000, 256_000);
        //assert!(bins.contains(&(585 + 1))); // Should contain level 4 bins
        //assert!(bins.contains(&(585 + 2)));

        //// Test larger range that requires checking higher level bins
        //let bins = index.region_to_bins(0, 1_000_000);
        //// Should contain bins from multiple levels
        //assert!(bins.contains(&0)); // Level 0 bin
        //assert!(bins.contains(&1)); // Level 1 bin
        //assert!(bins.iter().any(|&x| x >= 585)); // Should have some level 4 bins

        //// Test exact bin boundary
        //let bins = index.region_to_bins(0, 128 * 1024); // Exactly one 128kb bin
        //assert!(bins.contains(&585)); // Should contain first level 4 bin

        //// Test that bins are unique
        //let bins = index.region_to_bins(0, 10_000_000);
        //let unique_bins: Vec<_> = bins.iter().cloned().collect();
        //assert_eq!(bins.len(), unique_bins.len());
    }
}
