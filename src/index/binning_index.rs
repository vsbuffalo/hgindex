// binning_index.rs

use std::{fs::File, io::BufWriter, path::Path};

use super::binning::{BinningSchema, HierarchicalBins};
use crate::SerdeType;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};

/// BinningIndex is the sequence-level (e.g. chromosome) container
/// for SequenceIndex objects that index the features.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct BinningIndex<M>
where
    M: SerdeType + std::fmt::Debug,
{
    pub schema: BinningSchema,
    pub sequences: FxHashMap<String, SequenceIndex>,
    pub metadata: Option<M>,
    pub use_linear_index: bool,
}

// Implement Deserialize manually to handle the lifetime bounds correctly
impl<'de, M> Deserialize<'de> for BinningIndex<M>
where
    M: SerdeType + std::fmt::Debug,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let (schema, sequences, metadata, use_linear_index) =
            serde::Deserialize::deserialize(deserializer)?;

        Ok(BinningIndex {
            schema,
            sequences,
            metadata,
            use_linear_index,
        })
    }
}

/// SequenceIndex stores the bin indices to the features they
/// contain fully.
///
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SequenceIndex {
    // Map from bin ID to u64, which can be used as a VirtualOffset.
    pub bins: FxHashMap<u32, (u64, u64)>,
    // Linear index for quick region queries
    pub linear_index: Vec<u64>,
}

impl<M> Default for BinningIndex<M>
where
    M: SerdeType + std::fmt::Debug,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<M> BinningIndex<M>
where
    M: SerdeType + std::fmt::Debug,
{
    pub fn new() -> Self {
        BinningIndex {
            schema: BinningSchema::default(),
            sequences: FxHashMap::default(),
            metadata: None,
            use_linear_index: true,
        }
    }

    pub fn from_schema(schema: &BinningSchema) -> Self {
        BinningIndex {
            schema: schema.clone(),
            sequences: FxHashMap::default(),
            metadata: None,
            use_linear_index: true,
        }
    }

    pub fn disable_linear_index(&mut self) {
        self.use_linear_index = false;
    }

    pub fn enable_linear_index(&mut self) {
        self.use_linear_index = true;
    }

    /// Create a new index object by reading a binary serialized version of disk.
    pub fn open(path: &Path) -> std::result::Result<Self, Box<dyn std::error::Error>> {
        let file = File::open(path)?;
        let mmap = unsafe { memmap2::Mmap::map(&file)? };
        Ok(bincode::deserialize(&mmap[..])?)
    }

    /// Set the metadata
    pub fn set_metadata(&mut self, metadata: M) {
        self.metadata = Some(metadata);
    }

    /// Add a feature, a range with a file
    pub fn add_feature(&mut self, chrom: &str, start: u32, end: u32, index: u64) {
        let binning = HierarchicalBins::from_schema(&self.schema);
        let bin_id = binning.region_to_bin(start, end);

        // Get or create the sequence index.
        let bin_index = self.sequences.entry(chrom.to_string()).or_default();

        let bin_range = bin_index.bins.entry(bin_id).or_insert((u64::MAX, 0));
        bin_range.0 = bin_range.0.min(index); // Update the start offset of the range
        bin_range.1 = bin_range.1.max(index); // Update the end offset of the range

        if self.use_linear_index {
            let start_window = start >> binning.linear_shift;
            let end_window = (end - 1) >> binning.linear_shift;

            if end_window >= bin_index.linear_index.len() as u32 {
                bin_index
                    .linear_index
                    .resize(end_window as usize + 1, u64::MAX);
            }

            for window in start_window..=end_window {
                bin_index.linear_index[window as usize] =
                    bin_index.linear_index[window as usize].min(index);
            }
        }
    }

    /// Return the range of offsets that could contain an overlapping range.
    pub fn get_candidate_offsets(&self, chrom: &str, start: u32, end: u32) -> Option<(u64, u64)> {
        let chrom_index = self.sequences.get(chrom)?;
        let binning = HierarchicalBins::from_schema(&self.schema);
        let bins_to_check = binning.region_to_bins(start, end);

        // Calculate minimum offset using linear index
        let min_linear_offset = if self.use_linear_index {
            let start_window = start >> binning.linear_shift;
            let end_window = (end - 1) >> binning.linear_shift;
            let max_window = chrom_index.linear_index.len() as u32;

            let mut min = u64::MAX;
            for window in start_window..=end_window.min(max_window - 1) {
                min = min.min(chrom_index.linear_index[window as usize]);
            }
            min
        } else {
            0 // No linear index
        };

        // Initialize result range
        let mut min_offset = u64::MAX;
        let mut max_offset = 0;

        // Look through all potentially overlapping bins
        for bin_id in bins_to_check {
            if let Some(&(bin_min, bin_max)) = chrom_index.bins.get(&bin_id) {
                if bin_min != u64::MAX {
                    // Check if bin contains any features
                    min_offset = min_offset.min(bin_min);
                    max_offset = max_offset.max(bin_max);
                }
            }
        }

        // Use linear index minimum if available
        if min_linear_offset != u64::MAX {
            min_offset = min_offset.max(min_linear_offset);
        }

        if min_offset != u64::MAX && max_offset >= min_offset {
            Some((min_offset, max_offset))
        } else {
            None
        }
    }

    /// Write the BinningIndex to a path by binary serialization.
    pub fn write(&mut self, path: &Path) -> std::result::Result<(), Box<dyn std::error::Error>> {
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
            bins: FxHashMap::default(),
            linear_index: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_index() -> BinningIndex<()> {
        let mut index = BinningIndex::new();

        // Add features with the indices
        index.add_feature("chr1", 1000, 2000, 100);
        index.add_feature("chr1", 1500, 2500, 200);
        index.add_feature("chr1", 5000, 6000, 300);

        index
    }

    #[test]
    fn test_candidate_ranges() {
        let index = make_test_index();

        // Test finding candidates that might contain overlaps
        let range = index.get_candidate_offsets("chr1", 1750, 2250).unwrap();
        assert!(range.0 <= 100); // Should include offset of first feature
        assert!(range.1 >= 200); // Should include offset of second feature

        // Test region far from any bins - should return None
        let range = index.get_candidate_offsets("chr1", 1_000_000, 1_100_000);
        assert!(range.is_none());

        // Test chromosome boundaries
        let range = index.get_candidate_offsets("chr2", 1000, 2000);
        assert!(range.is_none());
    }

    #[test]
    fn test_linear_index() {
        let mut index = BinningIndex::<()>::new();

        // Add features with increasing offsets in bin-sized chunks
        for i in 0..10u32 {
            let start = i * 16384; // Use bin size as interval
            let end = (i + 1) * 16384;
            let offset = u64::from(i * 100);
            index.add_feature("chr1", start, end, offset);
        }

        // Test that linear index provides valid bounds
        let range = index
            .get_candidate_offsets("chr1", 5 * 16384, 6 * 16384)
            .unwrap();
        assert!(range.0 <= 500); // Should include the target feature's offset
        assert!(range.1 >= 500); // Should include at least this feature
    }

    #[test]
    fn test_bin_boundaries() {
        let mut index = BinningIndex::<()>::new();

        // Add features at bin boundaries
        index.add_feature("chr1", 0, 16384, 100); // Exactly one level 1 bin
        index.add_feature("chr1", 16384, 32768, 200); // Next level 1 bin

        // Query across bin boundary - both bins should be included as candidates
        let range = index.get_candidate_offsets("chr1", 16000, 17000).unwrap();
        assert!(range.0 <= 100); // Should include first feature
        assert!(range.1 >= 200); // Should include second feature
    }
}
