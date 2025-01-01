// binning_index.rs

use std::{
    fs::File,
    io::{BufReader, BufWriter},
    path::Path,
};

use crate::{error::HgIndexError, SerdeType};

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
    pub fn add_feature(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
        index: u64,
    ) -> Result<(), HgIndexError> {
        if start >= end {
            return Err(HgIndexError::InvalidInterval { start, end });
        }
        let binning = HierarchicalBins::default();
        // println!(
        //     "Adding feature to index: {}:{}-{} at offset {}",
        //     chrom, start, end, index
        // );
        let bin_id = binning.region_to_bin(start, end);
        // println!("Assigned to bin: {}", bin_id);

        // Get or create a new SequenceIndex for this chromosome.
        let bin_index = self.sequences.entry(chrom.to_string()).or_default();

        // Add feature to bin
        let features = bin_index.bins.entry(bin_id).or_default();
        // println!("Bin {} now has {} features", bin_id, features.len() + 1);
        features.push(Feature { start, end, index });

        // Update linear index with bounds check
        let linear_idx = start >> binning.linear_shift;
        if linear_idx >= bin_index.linear_index.len() as u32 {
            // Need to resize to accommodate new start position
            let new_size = linear_idx as usize + 1;
            bin_index.linear_index.resize(new_size, u64::MAX);
        }
        bin_index.linear_index[linear_idx as usize] =
            bin_index.linear_index[linear_idx as usize].min(index);

        Ok(())
    }

    /// Return the indices (e.g. file offsets) of all ranges that overlap with the supplied range.
    pub fn find_overlapping(&self, chrom: &str, start: u32, end: u32) -> Vec<u64> {
        let Some(chrom_index) = self.sequences.get(chrom) else {
            return vec![];
        };

        let binning = HierarchicalBins::default();
        let bins_to_check = binning.region_to_bins(start, end);

        let mut overlapping = Vec::new();
        let min_offset = {
            let linear_idx = start >> binning.linear_shift;
            if linear_idx < chrom_index.linear_index.len() as u32 {
                let offset = chrom_index.linear_index[linear_idx as usize];
                if offset == u64::MAX {
                    0 // Use 0 for uninitialized entries
                } else {
                    offset
                }
            } else {
                0 // Use 0 for out of range
            }
        };

        for bin_id in bins_to_check {
            if let Some(features) = chrom_index.bins.get(&bin_id) {
                // println!("Bin {} has {} features", bin_id, features.len());
                for feature in features {
                    if feature.index >= min_offset && feature.start < end && feature.end > start {
                        // println!(
                        //     "Found overlapping feature: {}-{}",
                        //     feature.start, feature.end
                        // );
                        overlapping.push(feature.index);
                    }
                }
            } else {
                // println!("No features in bin {}", bin_id);
            }
        }

        overlapping
    }

    /// Write the BinningIndex to a path by binary serialization.
    pub fn write(&self, path: &Path) -> Result<(), HgIndexError> {
        let mut file = BufWriter::new(File::create(path)?);
        bincode::serialize_into(&mut file, &self)
            .map_err(|e| HgIndexError::SerializationError(e.to_string()))?;
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
    use tempfile;

    fn make_test_case_01() -> BinningIndex<()> {
        let mut index = BinningIndex::new();

        // Add some features with the indices
        index.add_feature("chr1", 1000, 2000, 100).unwrap();
        index.add_feature("chr1", 1500, 2500, 200).unwrap();
        index.add_feature("chr1", 5000, 6000, 300).unwrap();
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
    fn test_edge_case_binning() {
        let mut index = BinningIndex::<()>::new();
        // Add features testing various edge cases
        index.add_feature("chr1", 1000, 1001, 100).unwrap(); // 1-base instead of zero-length
        index.add_feature("chr1", 2000, 2001, 200).unwrap(); // 1-base
        index.add_feature("chr1", 3000, 3100, 300).unwrap(); // regular feature
        index.add_feature("chr1", 4000, 4001, 400).unwrap(); // 1-base instead of zero-length

        // Test one-base length queries
        let overlaps = index.find_overlapping("chr1", 1000, 1001);
        assert_eq!(overlaps.len(), 1);
        assert!(overlaps.contains(&100));

        // Test exact boundary matches
        let overlaps = index.find_overlapping("chr1", 2000, 2001);
        assert_eq!(overlaps.len(), 1);
        assert!(overlaps.contains(&200));

        // No overlaps here, since 2001 of the range is right exclusive
        let overlaps = index.find_overlapping("chr1", 2001, 2002);
        assert_eq!(overlaps.len(), 0);

        // Test regular overlaps
        let overlaps = index.find_overlapping("chr1", 3050, 3075);
        assert_eq!(overlaps.len(), 1);
        assert!(overlaps.contains(&300));

        // Test boundary overlaps
        let overlaps = index.find_overlapping("chr1", 3000, 3001);
        assert_eq!(overlaps.len(), 1);
        assert!(overlaps.contains(&300));

        let overlaps = index.find_overlapping("chr1", 3099, 3100);
        assert_eq!(overlaps.len(), 1);
        assert!(overlaps.contains(&300));
    }

    #[test]
    fn test_invalid_intervals() {
        let mut index = BinningIndex::<()>::new();

        // Test zero-width intervals
        let err = index.add_feature("chr1", 1000, 1000, 100).unwrap_err();
        assert!(matches!(
            err,
            HgIndexError::InvalidInterval {
                start: 1000,
                end: 1000
            }
        ));

        // Test where end < start
        let err = index.add_feature("chr1", 2000, 1999, 200).unwrap_err();
        assert!(matches!(
            err,
            HgIndexError::InvalidInterval {
                start: 2000,
                end: 1999
            }
        ));

        // Test boundary case where start is max value
        let err = index
            .add_feature("chr1", u32::MAX, u32::MAX, 300)
            .unwrap_err();
        assert!(matches!(
            err,
            HgIndexError::InvalidInterval {
                start: u32::MAX,
                end: u32::MAX
            }
        ));

        // Valid intervals should work
        assert!(index.add_feature("chr1", 1000, 1001, 100).is_ok()); // Minimum valid interval (1bp)
        assert!(index.add_feature("chr1", 0, 1, 200).is_ok()); // Valid at zero
        assert!(index.add_feature("chr1", 1000, 2000, 300).is_ok()); // Larger interval
    }

    //#[test]
    //fn test_bin_boundary_overlaps() {
    //    let mut index = BinningIndex::<()>::new();
    //    let binning = HierarchicalBins::default();
    //    let bin_size = 1 << binning.base_shift; // 128kb bin size

    //    // Add feature that spans bin boundary
    //    // If bin_size is 128kb (131072), this creates a feature spanning:
    //    // [131071, 131073) - which crosses from bin 0 to bin 1
    //    index.add_feature("chr1", bin_size - 1, bin_size + 1, 100);

    //    // Test overlap on left side of boundary
    //    // [131071, 131072) - should be in bin 0
    //    let overlaps = index.find_overlapping("chr1", bin_size - 1, bin_size);
    //    assert_eq!(overlaps.len(), 1);
    //    assert!(overlaps.contains(&100));

    //    // Test overlap on right side of boundary
    //    // [131072, 131073) - should be in bin 1
    //    let overlaps = index.find_overlapping("chr1", bin_size, bin_size + 1);
    //    assert_eq!(overlaps.len(), 1);
    //    assert!(overlaps.contains(&100));

    //    // Test overlap spanning boundary
    //    // [131071, 131073) - should span bins 0 and 1
    //    let overlaps = index.find_overlapping("chr1", bin_size - 1, bin_size + 1);
    //    assert_eq!(overlaps.len(), 1);
    //    assert!(overlaps.contains(&100));
    //}

    #[test]
    fn test_bed_coordinate_system() {
        let mut index = BinningIndex::<()>::new();

        // Create several strategic features
        index.add_feature("chr1", 10, 20, 1).unwrap(); // [10,20)
        index.add_feature("chr1", 20, 30, 2).unwrap(); // [20,30) - adjacent to first feature
        index.add_feature("chr1", 40, 41, 3).unwrap(); // [40,41) - single base feature

        // Test 1: Adjacent features should each overlap a query that spans their boundary
        let overlaps = index.find_overlapping("chr1", 19, 21);
        assert_eq!(overlaps.len(), 2); // Should find both features

        // Test 2: Point query exactly at feature boundary should only match feature that starts there
        let overlaps = index.find_overlapping("chr1", 20, 21);
        assert_eq!(overlaps.len(), 1); // Should only find second feature
        assert!(overlaps.contains(&2)); // Should be second feature
        assert!(!overlaps.contains(&1)); // Should not find first feature

        // Test 3: Point queries around single-base feature
        let overlaps = index.find_overlapping("chr1", 40, 41);
        assert_eq!(overlaps.len(), 1); // Should find the feature

        let overlaps = index.find_overlapping("chr1", 41, 42);
        assert_eq!(overlaps.len(), 0); // Should not find the feature

        // Test 4: One-base query at the end of first feature
        let overlaps = index.find_overlapping("chr1", 19, 20);
        assert_eq!(overlaps.len(), 1); // Should find first feature
        assert!(overlaps.contains(&1));

        // Test 5: Test non-overlapping ranges
        let overlaps = index.find_overlapping("chr1", 30, 40);
        assert_eq!(overlaps.len(), 0); // Should find no features
    }
}
