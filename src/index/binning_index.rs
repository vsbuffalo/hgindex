// binning_index.rs

use std::collections::HashMap;
use std::{fs::File, io::BufWriter, path::Path};

use super::binning::{BinningSchema, HierarchicalBins};
use crate::Metadata;
use crate::{error::HgIndexError, SerdeType};
use dashmap::DashMap;
use parking_lot::Mutex;
use rustc_hash::FxHashMap;
use serde::ser::SerializeStruct;
use serde::{Deserialize, Deserializer, Serialize};

/// State of chromosome reading / buffer for use as
/// DashMap value.
#[derive(Debug, Default)]
pub struct ChromState {
    buffer: Mutex<Vec<u64>>,
    last_start: Option<u32>,
}

impl ChromState {
    pub fn check_sorted(&self, start: u32, chrom: &str) -> Result<(), HgIndexError> {
        if let Some(last_start) = self.last_start {
            if last_start > start {
                return Err(HgIndexError::UnsortedFeatures {
                    chrom: chrom.to_string(),
                    bin_id: 0,
                    previous: last_start,
                    current: start,
                });
            }
        }
        Ok(())
    }

    pub fn update_position(&mut self, start: u32) {
        self.last_start = Some(start);
    }

    pub fn clear_buffer(&self) {
        self.buffer.lock().clear();
    }

    pub fn get_buffer(&self) -> parking_lot::MutexGuard<'_, Vec<u64>> {
        self.buffer.lock()
    }
}

/// BinningIndex is the sequence-level (e.g. chromosome) container
/// for SequenceIndex objects that index the features.
#[derive(Debug)]
pub struct BinningIndex {
    pub schema: BinningSchema,
    pub sequences: DashMap<String, SequenceIndex>,
    pub use_linear_index: bool,
}

impl Serialize for BinningIndex {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut state = serializer.serialize_struct("BinningIndex", 6)?;

        let items = self.sequences.iter();
        let map: HashMap<String, SequenceIndex> =
            HashMap::from_iter(items.map(|r| (r.key().clone(), r.value().clone())));

        state.serialize_field("schema", &self.schema)?;
        state.serialize_field("sequences", &map)?;
        state.serialize_field("use_linear_index", &self.use_linear_index)?;

        state.end()
    }
}

// Implement Deserialize manually to handle the lifetime bounds correctly
impl<'de> Deserialize<'de> for BinningIndex {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let (schema, sequences, use_linear_index): (
            BinningSchema,
            FxHashMap<String, SequenceIndex>, // Use FxHashMap for deserialization
            bool,
        ) = serde::Deserialize::deserialize(deserializer)?;

        // Convert FxHashMap to DashMap
        let sequences: DashMap<String, SequenceIndex> = sequences.into_iter().collect();

        Ok(BinningIndex {
            schema,
            sequences,
            use_linear_index,
        })
    }
}

/// SequenceIndex stores the bin indices to the features they
/// contain fully.
///
#[derive(Debug, Serialize, Deserialize)]
pub struct SequenceIndex {
    // Map from bin ID to u64, which can be used as a VirtualOffset.
    pub bins: FxHashMap<u32, Vec<Feature>>,
    // Linear index for quick region queries.
    pub linear_index: Vec<u64>,
    // Chromosome-specific state.
    #[serde(skip)]
    pub chrom_state: ChromState,
}

impl Clone for SequenceIndex {
    fn clone(&self) -> Self {
        Self {
            bins: self.bins.clone(),
            linear_index: self.linear_index.clone(),
            chrom_state: ChromState::default(),
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Feature {
    /// Start position.
    start: u32,
    /// End position.
    end: u32,
    /// The feature index (e.g. a file offset).
    index: u64,
}

impl Default for BinningIndex {
    fn default() -> Self {
        Self::new()
    }
}

impl BinningIndex {
    pub fn new() -> Self {
        BinningIndex {
            schema: BinningSchema::default(),
            sequences: DashMap::default(),
            use_linear_index: true,
        }
    }

    pub fn from_schema(schema: &BinningSchema) -> Self {
        BinningIndex {
            schema: schema.clone(),
            sequences: DashMap::default(),
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

    pub fn get_chrom_state(
        &self,
        chrom: &str,
    ) -> Option<dashmap::mapref::one::Ref<'_, String, SequenceIndex>> {
        self.sequences.get(chrom)
    }

    pub fn check_sorted(&self, chrom: &str, start: u32) -> Result<(), HgIndexError> {
        if let Some(sequence_index) = self.sequences.get(chrom) {
            sequence_index.chrom_state.check_sorted(start, chrom)?;
        }
        Ok(())
    }

    pub fn update_position_state(&self, chrom: &str, start: u32) {
        if let Some(mut sequence_index) = self.sequences.get_mut(chrom) {
            sequence_index.chrom_state.update_position(start);
        }
    }

    pub fn clear_chrom_buffer(&self, chrom: &str) {
        if let Some(sequence_index) = self.sequences.get_mut(chrom) {
            sequence_index.chrom_state.clear_buffer();
        }
    }

    /// Add a feature to the index at the `BinningIndex` level.
    pub fn add_feature(
        &self,
        chrom: &str,
        start: u32,
        end: u32,
        index: u64,
    ) -> Result<(), HgIndexError> {
        // Get or insert the `SequenceIndex` for the chromosome.
        let mut sequence_index = self
            .sequences
            .entry(chrom.to_string())
            .or_insert_with(SequenceIndex::new);

        // Add the feature to the `SequenceIndex`.
        let schema = &self.schema;
        sequence_index.add_feature(start, end, index, schema)
    }

    // Helper method: Compute the minimum offset from the linear index
    fn get_min_offset<'a>(
        &self,
        sequence_index: &SequenceIndex,
        start: u32,
        binning: &HierarchicalBins,
    ) -> u64 {
        let linear_idx = start >> binning.linear_shift;
        if linear_idx >= sequence_index.linear_index.len() as u32 {
            return u64::MAX;
        }
        sequence_index.linear_index[linear_idx as usize]
    }

    /// Return an iterator over the indices (e.g. file offsets) of all ranges that overlap with the supplied range.
    pub fn find_overlapping<'a>(
        &'a self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Box<dyn Iterator<Item = u64> + 'a> {
        if let Some(sequence_index) = self.sequences.get(chrom) {
            let binning = HierarchicalBins::from_schema(&self.schema);
            let bins_to_check = binning.region_to_bins(start, end);

            let iter = bins_to_check
                .into_iter()
                .filter_map(move |bin_id| sequence_index.bins.get(&bin_id).cloned())
                .flat_map(move |features| {
                    features
                        .into_iter()
                        .filter(move |feature| feature.start < end && feature.end > start)
                        .map(|feature| feature.index)
                });

            Box::new(iter)
        } else {
            Box::new(std::iter::empty())
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
            if let Some(features) = chrom_index.bins.get(&bin_id) {
                for feature in features {
                    // Check if bin contains any features
                    min_offset = min_offset.min(feature.index);
                    max_offset = max_offset.max(feature.index);
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
    // TODO write metdata
    pub fn serialize_to_path<M>(
        &self,
        path: &Path,
        metadata: Option<&M>,
    ) -> Result<(), HgIndexError>
    where
        M: Metadata,
    {
        let mut file = BufWriter::new(File::create(path)?);

        let index_store = BinningIndexStore::new(&self, metadata);
        bincode::serialize_into(&mut file, &index_store)?;
        Ok(())
    }
}

/// A store-variant of SequenceIndex (with metadata)
#[derive(Debug, Serialize)]
pub struct BinningIndexStore<'a, M>
where
    M: Metadata,
{
    pub store: &'a BinningIndex,
    pub metadata: Option<&'a M>,
}

impl<'a, M> BinningIndexStore<'a, M>
where
    M: SerdeType + std::fmt::Debug,
{
    pub fn new(store: &'a BinningIndex, metadata: Option<&'a M>) -> Self {
        BinningIndexStore { store, metadata }
    }
}

impl<'de, M> Deserialize<'de> for BinningIndexStore<'de, M>
where
    M: Metadata,
    for<'a> &'a M: Metadata,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct BinningIndexStoreData<M> {
            schema: BinningSchema,
            sequences: FxHashMap<String, SequenceIndex>,
            use_linear_index: bool,
            metadata: M,
        }

        let data = BinningIndexStoreData::deserialize(deserializer)?;

        // Reconstruct DashMap from the deserialized FxHashMap
        let sequences: DashMap<String, SequenceIndex> = data.sequences.into_iter().collect();

        let store = BinningIndex {
            schema: data.schema,
            sequences,
            use_linear_index: data.use_linear_index,
        };

        Ok(BinningIndexStore {
            store: Box::leak(Box::new(store)), // Lifetime adjustment
            metadata: *Box::leak(Box::new(data.metadata)), // Lifetime adjustment
        })
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
            chrom_state: ChromState::default(),
        }
    }

    /// Add a feature to the index
    // DEVNOTE: this was initially at the BinningIndex level,
    // but this method created a headache at that level. Moving
    // it down to SequenceIndex solved this.
    pub fn add_feature(
        &mut self,
        start: u32,
        end: u32,
        index: u64,
        schema: &BinningSchema,
    ) -> Result<(), HgIndexError> {
        // Check sorting
        self.chrom_state.check_sorted(start, "chrom")?;

        // Update tracking
        self.chrom_state.update_position(start);

        let binning = HierarchicalBins::from_schema(schema);
        let bin_id = binning.region_to_bin(start, end);

        // Add the feature
        let feature = Feature { start, end, index };
        self.bins.entry(bin_id).or_default().push(feature);

        // Update the linear index
        let linear_idx = start >> binning.linear_shift;
        if linear_idx >= self.linear_index.len() as u32 {
            self.linear_index.resize(linear_idx as usize + 1, u64::MAX);
        }
        self.linear_index[linear_idx as usize] = self.linear_index[linear_idx as usize].min(index);

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_index() -> BinningIndex {
        let index = BinningIndex::default();

        // Add features with the indices
        index.add_feature("chr1", 1000, 2000, 100).unwrap();
        index.add_feature("chr1", 1500, 2500, 200).unwrap();
        index.add_feature("chr1", 5000, 6000, 300).unwrap();

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
        let index = BinningIndex::new();

        // Add features with increasing offsets in bin-sized chunks
        for i in 0..10u32 {
            let start = i * 16384; // Use bin size as interval
            let end = (i + 1) * 16384;
            let offset = u64::from(i * 100);
            index.add_feature("chr1", start, end, offset).unwrap();
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
        let index = BinningIndex::new();

        // Add features at bin boundaries
        index.add_feature("chr1", 0, 16384, 100).unwrap(); // Exactly one level 1 bin
        index.add_feature("chr1", 16384, 32768, 200).unwrap(); // Next level 1 bin

        // Query across bin boundary - both bins should be included as candidates
        let range = index.get_candidate_offsets("chr1", 16000, 17000).unwrap();
        assert!(range.0 <= 100); // Should include first feature
        assert!(range.1 >= 200); // Should include second feature
    }

    #[test]
    fn test_feature_ordering() {
        let index = BinningIndex::new();

        // Test basic ordering works
        assert!(index.add_feature("chr1", 1000, 2000, 100).is_ok());
        assert!(index.add_feature("chr1", 2000, 3000, 200).is_ok());

        // Test same start different end
        assert!(index.add_feature("chr2", 1000, 2000, 300).is_ok());
        assert!(index.add_feature("chr2", 1000, 2500, 400).is_ok());

        // Test chromosome transitions
        assert!(index.add_feature("chr2", 5000, 6000, 500).is_ok());
        assert!(index.add_feature("chr3", 1000, 2000, 600).is_ok());

        // Test errors
        let index = BinningIndex::new();

        // Setup initial feature
        assert!(index.add_feature("chr1", 2000, 3000, 100).is_ok());

        // Test earlier start fails
        assert!(matches!(
            index.add_feature("chr1", 1000, 2000, 200),
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
        let index = BinningIndex::new();

        // Test ties are okay
        assert!(index.add_feature("chr1", 1000, 2000, 100).is_ok());
        assert!(index.add_feature("chr1", 1000, 2000, 200).is_ok()); // Exact tie ok
        assert!(index.add_feature("chr1", 1000, 2000, 300).is_ok()); // Multiple ties ok

        // Test we can continue after ties
        assert!(index.add_feature("chr1", 1000, 2500, 400).is_ok()); // Same start, longer end
        assert!(index.add_feature("chr1", 2000, 3000, 500).is_ok()); // Progress after ties

        // Test chromosome transitions with ties
        assert!(index.add_feature("chr2", 1000, 2000, 600).is_ok());
        assert!(index.add_feature("chr2", 1000, 2000, 700).is_ok()); // Tie on new chrom
    }
}
