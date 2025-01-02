// binning_index.rs

use std::{
    fs::File,
    io::{BufReader, BufWriter},
    path::Path,
};

use crate::SerdeType;

use super::binning::HierarchicalBins;
use indexmap::IndexMap;
use serde::{Deserialize, Serialize};

/// BinningIndex is the sequence-level (e.g. chromosome) container
/// for SequenceIndex objects that index the features.
#[derive(Debug, Serialize, PartialEq)]
pub struct BinningIndex<M>
where
    M: SerdeType,
{
    pub sequences: IndexMap<String, SequenceIndex>,
    pub metadata: Option<M>,
}

// Implement Deserialize manually to handle the lifetime bounds correctly
impl<'de, M> Deserialize<'de> for BinningIndex<M>
where
    M: SerdeType,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        // Create a temporary struct for deriving
        #[derive(Deserialize)]
        struct BinningIndexHelper<M> {
            sequences: IndexMap<String, SequenceIndex>,
            metadata: Option<M>,
        }

        let helper = BinningIndexHelper::deserialize(deserializer)?;
        Ok(BinningIndex {
            sequences: helper.sequences,
            metadata: helper.metadata,
        })
    }
}

/// SequenceIndex stores the bin indices to the features they
/// contain fully.
///
#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct SequenceIndex {
    // Map from bin ID to list of features in that bin
    pub bins: IndexMap<u32, Vec<Feature>>,
    // Linear index for quick region queries
    pub linear_index: Vec<u64>,
}

#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct Feature {
    /// Start position.
    start: u32,
    /// End position.
    end: u32,
    /// The feature index (e.g. a file offset).
    index: u64,
}

impl<M: SerdeType> Default for BinningIndex<M> {
    fn default() -> Self {
        Self::new()
    }
}

impl<M: SerdeType> BinningIndex<M> {
    pub fn new() -> Self {
        BinningIndex {
            sequences: IndexMap::new(),
            metadata: None,
        }
    }

    /// Create a new index object by reading a binary serialized version of disk.
    pub fn open(path: &Path) -> std::result::Result<Self, Box<dyn std::error::Error>> {
        let mut file = BufReader::new(File::open(path)?);
        let obj = bincode::deserialize_from(&mut file).unwrap();
        Ok(obj)
    }

    /// Set the metadata
    pub fn set_metadata(&mut self, metadata: M) {
        self.metadata = Some(metadata);
    }

    /// Add a feature, a range with a file
    pub fn add_feature(&mut self, chrom: &str, start: u32, end: u32, index: u64) {
        let binning = HierarchicalBins::default();
        let bin_id = binning.region_to_bin(start, end);

        // the chromosome-level index
        let bin_index = self.sequences.entry(chrom.to_string()).or_default();

        // make the feature
        let feature = Feature { start, end, index };

        // insert the feature
        bin_index.bins.entry(bin_id).or_default().push(feature);

        // Update linear index
        let linear_idx = start >> binning.linear_shift;
        if linear_idx >= bin_index.linear_index.len() as u32 {
            bin_index
                .linear_index
                .resize(linear_idx as usize + 1, u64::MAX);
        }
        bin_index.linear_index[linear_idx as usize] =
            bin_index.linear_index[linear_idx as usize].min(index);
    }

    /// Return the indices (e.g. file offsets) of all ranges that overlap with the supplied range.
    pub fn find_overlapping(&self, chrom: &str, start: u32, end: u32) -> Vec<u64> {
        let Some(chrom_index) = self.sequences.get(chrom) else {
            return vec![];
        };

        let binning = HierarchicalBins::default();
        let bins_to_check = binning.region_to_bins(start, end);

        let mut overlapping = Vec::new();

        // Use linear index to find minimum file offset to start checking
        let min_offset = {
            let linear_idx = start >> binning.linear_shift;
            let max_linear_idx = chrom_index.linear_index.len();
            if linear_idx < max_linear_idx as u32 {
                chrom_index.linear_index[linear_idx as usize]
            } else {
                // If we're beyond the end of the linear index, there can't be any features
                // that overlap this region, so return an empty vector
                return vec![];
            }
        };

        for bin_id in bins_to_check {
            if let Some(features) = chrom_index.bins.get(&bin_id) {
                for feature in features {
                    if feature.index >= min_offset && feature.start < end && feature.end > start {
                        overlapping.push(feature.index);
                    }
                }
            }
        }

        overlapping.sort_unstable();
        overlapping.dedup();
        overlapping
    }

    /// Write the BinningIndex to a path by binary serialization.
    pub fn write(&self, path: &Path) -> std::result::Result<(), Box<dyn std::error::Error>> {
        let mut file = BufWriter::new(File::create(path)?);
        bincode::serialize_into(&mut file, &self)?;
        Ok(())
    }
}

impl Default for SequenceIndex {
    fn default() -> Self {
        Self::new()
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

    fn make_test_case_01() -> BinningIndex<()> {
        let mut index = BinningIndex::new();

        // Add some features with the indices
        index.add_feature("chr1", 1000, 2000, 100);
        index.add_feature("chr1", 1500, 2500, 200);
        index.add_feature("chr1", 5000, 6000, 300);
        index
    }

    #[test]
    fn test_index_storage() {
        let index = make_test_case_01();

        // Test finding overlaps
        let overlaps = index.find_overlapping("chr1", 1750, 2250);
        assert!(overlaps.contains(&100)); // First feature overlaps
        assert!(overlaps.contains(&200)); // Second feature overlaps
        assert!(!overlaps.contains(&300)); // Third feature doesn't overlap

        // Test linear index
        let chrom_index = index.sequences.get("chr1").unwrap();
        assert!(chrom_index.linear_index[0] == 100); // First offset in the 16kb chunk
    }

    #[test]
    fn test_persistent_storage() {
        let index = make_test_case_01();
        let tmp_dir = tempfile::tempdir().expect("Error creating tempdir");
        let index_file = tmp_dir.path().join("index.hgi");

        // write the index
        index.write(&index_file).expect("Error writing index");

        // read the index into new object
        let obj = BinningIndex::open(&index_file).expect("Error reading index");
        assert_eq!(index, obj);
    }

    #[test]
    fn test_bin_boundaries() {
        let index = HierarchicalBins::default();

        // Test ranges exactly matching bin sizes
        let bin1 = index.region_to_bin(0, 128 * 1024);
        let bin2 = index.region_to_bin(128 * 1024, 256 * 1024);

        // Test off-by-one ranges
        let bin3 = index.region_to_bin(0, (128 * 1024) - 1);
        let bin4 = index.region_to_bin(0, (128 * 1024) + 1);

        assert_ne!(bin1, bin2);
        assert_eq!(bin3, bin1);
        assert_ne!(bin4, bin1);
    }

    #[test]
    fn test_bed_interval_semantics() {
        let mut index = BinningIndex::<()>::new();

        // Test single-base features
        index.add_feature("chr1", 100, 101, 1); // One base at position 100

        // Test adjacent features
        index.add_feature("chr1", 200, 300, 2); // 100 bases
        index.add_feature("chr1", 300, 400, 3); // Next 100 bases, should not overlap

        // Test overlapping features
        index.add_feature("chr1", 500, 600, 4); // 100 bases
        index.add_feature("chr1", 550, 650, 5); // Overlaps with previous

        // Test zero-based coordinate system
        index.add_feature("chr1", 0, 10, 6); // First 10 bases

        // Test queries

        // Single base queries
        let overlaps = index.find_overlapping("chr1", 100, 101);
        assert_eq!(overlaps, vec![1], "Single base feature");

        // Query empty space between adjacent features
        let overlaps = index.find_overlapping("chr1", 299, 300);
        assert_eq!(overlaps, vec![2], "Should only contain first feature");
        let overlaps = index.find_overlapping("chr1", 300, 301);
        assert_eq!(overlaps, vec![3], "Should only contain second feature");

        // Query exact feature boundaries
        let overlaps = index.find_overlapping("chr1", 200, 300);
        assert_eq!(overlaps, vec![2], "Exact feature bounds");

        // Query overlapping region
        let overlaps = index.find_overlapping("chr1", 540, 560);
        assert_eq!(overlaps, vec![4, 5], "Overlapping features");

        // Query start of coordinate system
        let overlaps = index.find_overlapping("chr1", 0, 5);
        assert_eq!(overlaps, vec![6], "Zero-based start");

        // Single-base query at feature boundaries
        let overlaps = index.find_overlapping("chr1", 500, 501);
        assert_eq!(overlaps, vec![4], "Start boundary");
        let overlaps = index.find_overlapping("chr1", 599, 600);
        assert_eq!(overlaps, vec![4, 5], "End boundary");
    }

    #[test]
    fn test_bed_adjacent_intervals() {
        let mut index = BinningIndex::<()>::new();

        // Two adjacent features
        index.add_feature("chr1", 100, 200, 1); // [100,200)
        index.add_feature("chr1", 200, 300, 2); // [200,300)

        // Query the boundary
        let overlaps = index.find_overlapping("chr1", 200, 201);
        assert_eq!(overlaps, vec![2], "Should only match second feature");

        let overlaps = index.find_overlapping("chr1", 199, 200);
        assert_eq!(overlaps, vec![1], "Should only match first feature");
    }
}
