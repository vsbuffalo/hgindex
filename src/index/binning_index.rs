// binning_index.rs

use std::{fs::File, io::BufWriter, path::Path};

use super::binning::{BinningSchema, HierarchicalBins};
use crate::error::HgIndexError;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct LinearIndex {
    entries: Vec<u64>,
    shift: u32,
}

impl LinearIndex {
    // Create a new LinearIndex using a schema. Returns
    // None if this schema doesn't use a linear index.
    pub fn from_schema(bins: &HierarchicalBins) -> Option<Self> {
        bins.linear_shift.map(|shift| LinearIndex {
            entries: Vec::new(),
            shift,
        })
    }

    pub fn resize(&mut self, required_size: usize) {
        if self.entries.len() < required_size {
            self.entries.resize(required_size, u64::MAX);
        }
    }

    pub fn update(&mut self, start: u32, end: u32, offset: u64) {
        if end <= start {
            panic!(
                "Invalid range: start ({}) must be less than end ({})",
                start, end
            );
        }

        let start_window = start >> self.shift;
        let end_window = (end - 1) >> self.shift;
        self.resize((end_window + 1) as usize);

        for window in start_window..=end_window {
            self.entries[window as usize] = self.entries[window as usize].min(offset);
        }
    }

    pub fn get_min_offset(&self, start: u32) -> Option<u64> {
        let window = (start >> self.shift) as usize;
        self.entries.get(window).copied()
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }
}

/// BinningIndex is the sequence-level (e.g. chromosome) container
/// for SequenceIndex objects that index the features.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct BinningIndex {
    /// Initialized binning
    pub bins: HierarchicalBins,
    pub sequences: FxHashMap<String, SequenceIndex>,
    last_chrom: Option<String>,
    last_start: Option<u32>,
}

/// SequenceIndex stores the bin indices to the features they
/// contain fully.
#[derive(Debug, Serialize)]
pub struct SequenceIndex {
    // Map from bin ID to u64, which can be used as a VirtualOffset.
    pub bins: FxHashMap<u32, Vec<Feature>>,
    // Optional linear index for quick region queries
    pub linear_index: Option<LinearIndex>,
    // Buffer of (offsets, lengths)
    // overlap_buffer: Vec<(u64, u64)>,
}

impl Clone for SequenceIndex {
    fn clone(&self) -> Self {
        Self {
            bins: self.bins.clone(),
            linear_index: self.linear_index.clone(),
            // overlap_buffer: Vec::with_capacity(130),
        }
    }
}

impl PartialEq for SequenceIndex {
    fn eq(&self, other: &Self) -> bool {
        // Compare only bins and linear_index
        self.bins == other.bins && self.linear_index == other.linear_index
    }
}

impl<'de> Deserialize<'de> for SequenceIndex {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        // Helper struct to deserialize the parts we care about
        #[derive(Deserialize)]
        struct Helper {
            bins: FxHashMap<u32, Vec<Feature>>,
            linear_index: Option<LinearIndex>,
        }

        // Deserialize into helper
        let helper = Helper::deserialize(deserializer)?;

        // Construct full SequenceIndex with new buffer
        Ok(SequenceIndex {
            bins: helper.bins,
            linear_index: helper.linear_index,
            // overlap_buffer: Vec::with_capacity(130),
        })
    }
}

impl SequenceIndex {
    /// Create a new SequenceIndex from the pre-created HierarchicalBins.
    pub fn new(bins: &HierarchicalBins) -> Self {
        let linear_index = LinearIndex::from_schema(bins);
        SequenceIndex {
            bins: FxHashMap::default(),
            linear_index,
            // overlap_buffer: Vec::with_capacity(130),
        }
    }

    pub fn find_overlapping(
        &self,
        bins: &HierarchicalBins,
        start: u32,
        end: u32,
    ) -> Vec<(u64, u64)> {
        let min_offset = self
            .linear_index
            .as_ref()
            .and_then(|index| index.get_min_offset(start))
            .unwrap_or(0);

        // Pre-allocate results with an estimate based on bin count
        let estimated_capacity = bins.region_to_bins(start, end).len() * 10; // Assume ~10 features per bin
        let mut results = Vec::with_capacity(estimated_capacity);

        for &bin_id in bins.region_to_bins(start, end).iter() {
            if let Some(features) = self.bins.get(&bin_id) {
                // SIMD?
                // Filter features within the bin
                results.extend(features.iter().filter_map(|feature| {
                    if feature.index >= min_offset && feature.start < end && feature.end > start {
                        Some((feature.index, feature.length))
                    } else {
                        None
                    }
                }));
            }
        }

        results
    }

    /// Add a feature to the sequence index, ensuring it is in sorted order and updating bins and linear index.
    pub fn add_feature(
        &mut self,
        start: u32,
        end: u32,
        index: u64,
        bins: &HierarchicalBins,
        length: u64,
    ) -> Result<(), HgIndexError> {
        // Validate feature ordering
        if let Some(last_feature) = self.bins.values().flat_map(|f| f.iter()).last() {
            if start < last_feature.start {
                return Err(HgIndexError::UnsortedFeatures {
                    chrom: String::new(), // Chromosome validation occurs in BinningIndex
                    bin_id: 0,            // We could also calculate the bin ID here if helpful
                    previous: last_feature.start,
                    current: start,
                });
            }
        }

        // Determine the bin for the feature
        let bin_id = bins.region_to_bin(start, end);

        // Add the feature to the appropriate bin
        self.bins.entry(bin_id).or_default().push(Feature {
            start,
            end,
            index,
            length,
        });

        // Update the linear index
        if let Some(linear_index) = &mut self.linear_index {
            linear_index.update(start, end, index);
        }

        Ok(())
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Feature {
    /// Start position.
    pub start: u32,
    /// End position.
    pub end: u32,
    /// The feature index (e.g. a file offset).
    pub index: u64,
    /// The length of data in bytes.
    pub length: u64,
}

impl Default for BinningIndex {
    fn default() -> Self {
        let schema = BinningSchema::default();
        Self::new(&schema)
    }
}

impl BinningIndex {
    pub fn new(schema: &BinningSchema) -> Self {
        let bins = HierarchicalBins::from_schema(schema);
        BinningIndex {
            bins,
            sequences: FxHashMap::default(),
            last_chrom: None,
            last_start: None,
        }
    }

    pub fn get_sequence_index(&self, chrom: &str) -> Option<&SequenceIndex> {
        self.sequences.get(chrom)
    }

    pub fn disable_linear_index(&mut self) {
        // Clear out old linear indices.
        self.sequences
            .values_mut()
            .for_each(|f| f.linear_index = None);
        // Set the linear shift to None.
        self.bins.linear_shift = None;
    }

    pub fn has_linear_index(&self) -> bool {
        self.bins.linear_shift.is_some()
    }

    /// Create a new index object by reading a binary serialized version of disk.
    pub fn open(path: &Path) -> std::result::Result<Self, Box<dyn std::error::Error>> {
        let file = File::open(path)?;
        let mmap = unsafe { memmap2::Mmap::map(&file)? };
        let index: BinningIndex = bincode::deserialize(&mmap[..])?;
        Ok(index)
    }

    /// Add a feature, a range with a file
    pub fn add_feature(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
        index: u64,
        length: u64,
    ) -> Result<(), HgIndexError> {
        // Get or create the sequence index for the chromosome
        let sequence_index = self
            .sequences
            .entry(chrom.to_string())
            .or_insert_with(|| SequenceIndex::new(&self.bins));

        // Delegate the feature addition to SequenceIndex
        sequence_index.add_feature(start, end, index, &self.bins, length)?;

        Ok(())
    }

    /// Return the indices (e.g. file offsets) of all ranges that overlap with the supplied range.
    pub fn find_overlapping(&mut self, chrom: &str, start: u32, end: u32) -> Vec<(u64, u64)> {
        if let Some(chrom_index) = self.sequences.get_mut(chrom) {
            chrom_index.find_overlapping(&self.bins, start, end)
        } else {
            vec![]
        }
    }

    /// Write the BinningIndex to a path by binary serialization.
    pub fn serialize_to_path(
        &mut self,
        path: &Path,
    ) -> std::result::Result<(), Box<dyn std::error::Error>> {
        let mut file = BufWriter::new(File::create(path)?);
        bincode::serialize_into(&mut file, &self)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_feature_ordering() {
        let mut index = BinningIndex::default();

        // Test basic ordering works
        assert!(index.add_feature("chr1", 1000, 2000, 100, 7).is_ok());
        assert!(index.add_feature("chr1", 2000, 3000, 200, 7).is_ok());

        // Test same start different end
        assert!(index.add_feature("chr2", 1000, 2000, 300, 7).is_ok());
        assert!(index.add_feature("chr2", 1000, 2500, 400, 7).is_ok());

        // Test chromosome transitions
        assert!(index.add_feature("chr2", 5000, 6000, 500, 7).is_ok());
        assert!(index.add_feature("chr3", 1000, 2000, 600, 7).is_ok());

        // Test errors
        let mut index = BinningIndex::default();

        // Setup initial feature
        assert!(index.add_feature("chr1", 2000, 3000, 100, 7).is_ok());

        // Test earlier start fails
        assert!(matches!(
            index.add_feature("chr1", 1000, 2000, 200, 7),
            Err(HgIndexError::UnsortedFeatures {
                chrom: _,
                bin_id: _,
                previous: 2000,
                current: 1000
            })
        ));
    }

    #[test]
    fn test_feature_ordering_with_ties() {
        let mut index = BinningIndex::default();

        // Test ties are okay
        assert!(index.add_feature("chr1", 1000, 2000, 100, 0).is_ok());
        assert!(index.add_feature("chr1", 1000, 2000, 200, 0).is_ok()); // Exact tie ok
        assert!(index.add_feature("chr1", 1000, 2000, 300, 0).is_ok()); // Multiple ties ok

        // Test we can continue after ties
        assert!(index.add_feature("chr1", 1000, 2500, 400, 0).is_ok()); // Same start, longer end
        assert!(index.add_feature("chr1", 2000, 3000, 500, 0).is_ok()); // Progress after ties

        // Test chromosome transitions with ties
        assert!(index.add_feature("chr2", 1000, 2000, 600, 0).is_ok());
        assert!(index.add_feature("chr2", 1000, 2000, 700, 0).is_ok()); // Tie on new chrom
    }

    #[test]
    fn test_out_of_range_queries() {
        let mut index = BinningIndex::default();
        index.add_feature("chr1", 1000, 2000, 100, 0).unwrap();

        // Query completely out of range
        assert!(index.find_overlapping("chr1", 3000, 4000).is_empty());
        assert!(index.find_overlapping("chr1", 0, 500).is_empty());
    }

    #[test]
    fn test_spanning_features() {
        let mut index = BinningIndex::default();
        index.add_feature("chr1", 1000, 50000, 100, 0).unwrap(); // Spans multiple 16kb windows

        // Query a range covered by the feature
        let results = index.find_overlapping("chr1", 20000, 30000);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0], (100, 0));

        // Query a range outside the feature
        assert!(index.find_overlapping("chr1", 60000, 70000).is_empty());
    }

    #[test]
    fn test_large_scale_features() {
        let mut index = BinningIndex::default();

        // Add features in increments of 1kb
        for i in (0..1_000_000).step_by(1_000) {
            index.add_feature("chr1", i, i + 500, i as u64, 0).unwrap();
        }

        // Query a range covering multiple features
        let results = index.find_overlapping("chr1", 25000, 35000);
        assert_eq!(results.len(), 10); // Should find 10 features
    }

    #[test]
    fn test_disable_linear_index_consistency() {
        let mut index = BinningIndex::default();

        // Add features
        for i in (0..1_000_000).step_by(100_000) {
            index
                .add_feature("chr1", i, i + 50_000, i as u64, 0)
                .unwrap();
        }

        // Query with linear index
        let results_with_linear = index.find_overlapping("chr1", 100_000, 200_000).clone();
        assert!(index.has_linear_index());

        // Disable linear index
        index.disable_linear_index();
        assert!(!index.has_linear_index());

        // Query without linear index
        let results_without_linear = index.find_overlapping("chr1", 100_000, 200_000);

        // Ensure results are identical
        assert_eq!(results_with_linear, *results_without_linear);
    }

    #[test]
    fn test_linear_index_disable() {
        let mut index = BinningIndex::default();

        // First collect all features into a Vec
        let mut features = Vec::new();
        let mut offset = 0u64;

        // Create ~15000 features in the target region with some overlap
        for i in (100000..700000).step_by(17) {
            features.push((i, i + 100, offset));
            offset += 1;
            features.push((i, i + 100000, offset));
            offset += 1;
        }

        // Sort by start position
        features.sort_by_key(|(start, _, _)| *start);

        // Now add the sorted features
        for (start, end, offset) in features {
            index.add_feature("chr2", start, end, offset, 0).unwrap();
        }

        assert!(index.has_linear_index());
        let results_with_linear: Vec<(u64, u64)> =
            index.find_overlapping("chr2", 500000, 600000).to_vec();
        index.disable_linear_index();
        assert!(!index.has_linear_index());
        let results_without_linear: Vec<(u64, u64)> =
            index.find_overlapping("chr2", 500000, 600000).to_vec();

        assert_eq!(
            results_with_linear.len(),
            results_without_linear.len(),
            "Linear index results count {} vs no linear index count {}",
            results_with_linear.len(),
            results_without_linear.len()
        );
    }

    #[test]
    fn test_schema_persistence() {
        let schema = BinningSchema::Dense;
        let mut index = BinningIndex::new(&schema);
        let path = Path::new("test_index.hgidx");

        // Serialize
        index.serialize_to_path(&path).unwrap();

        // Deserialize
        let deserialized_index = BinningIndex::open(&path).unwrap();
        assert_eq!(deserialized_index.bins.schema, schema);

        // Clean up
        std::fs::remove_file(path).unwrap();
    }
}
