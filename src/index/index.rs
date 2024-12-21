// index.rs

use super::binning::HierarchicalBins;
use indexmap::IndexMap;

#[derive(Debug)]
pub struct BinningIndex {
    // Map from chromosome name/id to its features
    sequences: IndexMap<String, SequenceIndex>,
}

#[derive(Debug)]
pub struct SequenceIndex {
    // Map from bin ID to list of features in that bin
    bins: IndexMap<u32, Vec<Range>>,
    // Linear index for quick region queries
    linear_index: Vec<u64>, // File offsets every 16kb
}

#[derive(Debug)]
pub struct Range {
    start: u32,
    end: u32,
    file_offset: u64,
}

impl BinningIndex {
    pub fn new() -> Self {
        BinningIndex {
            sequences: IndexMap::new(),
        }
    }

    pub fn add_feature(&mut self, chrom: &str, start: u32, end: u32, file_offset: u64) {
        let binning = HierarchicalBins::default();
        let bin_id = binning.region_to_bin(start, end);

        let chrom_index = self
            .sequences
            .entry(chrom.to_string())
            .or_insert_with(SequenceIndex::new);

        let feature = Range {
            start,
            end,
            file_offset,
        };

        chrom_index
            .bins
            .entry(bin_id)
            .or_insert_with(Vec::new)
            .push(feature);

        // Update linear index
        let linear_idx = start >> 14; // 16kb chunks (2^14)
        if linear_idx >= chrom_index.linear_index.len() as u32 {
            chrom_index
                .linear_index
                .resize(linear_idx as usize + 1, u64::MAX);
        }
        chrom_index.linear_index[linear_idx as usize] =
            chrom_index.linear_index[linear_idx as usize].min(file_offset);
    }

    pub fn find_overlapping(&self, chrom: &str, start: u32, end: u32) -> Vec<u64> {
        let Some(chrom_index) = self.sequences.get(chrom) else {
            return vec![];
        };

        let binning = HierarchicalBins::default();
        let bins_to_check = binning.region_to_bins(start, end);

        let mut overlapping = Vec::new();

        // Use linear index to find minimum file offset to start checking
        let min_offset = {
            let linear_idx = start >> 14;
            if linear_idx < chrom_index.linear_index.len() as u32 {
                chrom_index.linear_index[linear_idx as usize]
            } else {
                return vec![];
            }
        };

        for bin_id in bins_to_check {
            if let Some(features) = chrom_index.bins.get(&bin_id) {
                for feature in features {
                    if feature.file_offset >= min_offset
                        && feature.start < end
                        && feature.end > start
                    {
                        overlapping.push(feature.file_offset);
                    }
                }
            }
        }

        overlapping.sort_unstable();
        overlapping.dedup();
        overlapping
    }
}

impl SequenceIndex {
    pub fn new() -> Self {
        SequenceIndex {
            bins: IndexMap::new(),
            linear_index: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_index_storage() {
        let mut index = BinningIndex::new();

        // Add some features
        index.add_feature("chr1", 1000, 2000, 100);
        index.add_feature("chr1", 1500, 2500, 200);
        index.add_feature("chr1", 5000, 6000, 300);

        // Test finding overlaps
        let overlaps = index.find_overlapping("chr1", 1750, 2250);
        assert!(overlaps.contains(&100)); // First feature overlaps
        assert!(overlaps.contains(&200)); // Second feature overlaps
        assert!(!overlaps.contains(&300)); // Third feature doesn't overlap

        // Test linear index
        let chrom_index = index.sequences.get("chr1").unwrap();
        assert!(chrom_index.linear_index[0] == 100); // First offset in the 16kb chunk
    }
}
