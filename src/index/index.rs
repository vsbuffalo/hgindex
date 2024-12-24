// index.rs

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
    // Linear index for quick region querieskjj
    pub linear_index: Vec<u64>, // File offsets every 16kb
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
        let bin_index = self
            .sequences
            .entry(chrom.to_string())
            .or_insert_with(SequenceIndex::new);

        // make the feature
        let feature = Feature { start, end, index };

        // insert the feature
        bin_index
            .bins
            .entry(bin_id)
            .or_insert_with(Vec::new)
            .push(feature);

        // Update linear index
        let linear_idx = start >> 14; // 16kb chunks (2^14)
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
    use tempfile;

    fn make_test_case_01() -> BinningIndex<()> {
        let mut index = BinningIndex::new();

        // Add some features with the indices
        index.add_feature("chr1", 1000, 2000, 100);
        index.add_feature("chr1", 1500, 2500, 200);
        index.add_feature("chr1", 5000, 6000, 300);
        return index;
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
}
