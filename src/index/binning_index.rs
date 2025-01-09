// binning_index.rs

use std::{
    fs::File,
    io::{BufReader, BufWriter},
    path::Path,
};

use crate::SerdeType;

use super::binning::{BinningSchema, HierarchicalBins};
use indexmap::IndexMap;
use serde::{Deserialize, Serialize};

/// BinningIndex is the sequence-level (e.g. chromosome) container
/// for SequenceIndex objects that index the features.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct BinningIndex<M>
where
    M: SerdeType,
{
    pub schema: BinningSchema,
    pub sequences: IndexMap<String, SequenceIndex>,
    pub metadata: Option<M>,
    pub use_linear_index: bool,
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
        // Helper struct has same fields but can derive Deserialize
        #[derive(Deserialize)]
        struct BinningIndexHelper<M> {
            binning: BinningSchema,
            sequences: IndexMap<String, SequenceIndex>,
            metadata: Option<M>,
            use_linear_index: bool,
        }

        // Deserialize into the helper struct
        let helper = BinningIndexHelper::deserialize(deserializer)?;

        // Convert helper into our actual type
        Ok(BinningIndex {
            schema: helper.binning,
            sequences: helper.sequences,
            metadata: helper.metadata,
            use_linear_index: helper.use_linear_index,
        })
    }
}

/// SequenceIndex stores the bin indices to the features they
/// contain fully.
///
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SequenceIndex {
    // Map from bin ID to list of features in that bin
    pub bins: IndexMap<u32, Vec<Chunk>>,
    // Linear index for quick region queries
    pub linear_index: Vec<u64>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Chunk {
    pub start_offset: u64, // Virtual offset where chunk starts
    pub end_offset: u64,   // Virtual offset where chunk ends
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Feature {
    /// Start position.
    pub start: u32,
    /// End position.
    pub end: u32,
    /// The feature index (e.g. a file offset).
    pub index: u64,
}

impl<M: SerdeType> Default for BinningIndex<M> {
    fn default() -> Self {
        Self::new()
    }
}

impl<M: SerdeType> BinningIndex<M> {
    pub fn new() -> Self {
        BinningIndex {
            schema: BinningSchema::default(),
            sequences: IndexMap::new(),
            metadata: None,
            use_linear_index: true,
        }
    }

    pub fn from_schema(schema: &BinningSchema) -> Self {
        BinningIndex {
            schema: schema.clone(),
            sequences: IndexMap::new(),
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
        let mut file = BufReader::new(File::open(path)?);
        let obj = bincode::deserialize_from(&mut file).unwrap();
        Ok(obj)
    }

    /// Set the metadata
    pub fn set_metadata(&mut self, metadata: M) {
        self.metadata = Some(metadata);
    }

    /// Add a feature, a range with a file
    pub fn add_feature(&mut self, chrom: &str, start: u32, end: u32, offset: u64) {
        let binning = HierarchicalBins::from_schema(&self.schema);
        let bin_id = binning.region_to_bin(start, end);

        let seq_index = self.sequences.entry(chrom.to_string()).or_default();

        // Get or create chunks for this bin
        let chunks = seq_index.bins.entry(bin_id).or_default();

        match chunks.last_mut() {
            // If we have a chunk and this feature is close enough, extend it
            Some(chunk) if offset - chunk.end_offset <= 16384 => {
                chunk.end_offset = offset;
            }
            // Otherwise start a new chunk
            _ => {
                chunks.push(Chunk {
                    start_offset: offset,
                    end_offset: offset,
                });
            }
        }

        // Update linear index
        let start_window = start >> binning.linear_shift;
        let end_window = (end - 1) >> binning.linear_shift; // -1 because end is exclusive

        if end_window >= seq_index.linear_index.len() as u32 {
            seq_index
                .linear_index
                .resize(end_window as usize + 1, u64::MAX);
        }

        for window in start_window..=end_window {
            seq_index.linear_index[window as usize] =
                seq_index.linear_index[window as usize].min(offset);
        }
    }

    pub fn find_overlapping(&self, chrom: &str, start: u32, end: u32) -> Vec<u64> {
        let Some(chrom_index) = self.sequences.get(chrom) else {
            return vec![];
        };

        let binning = HierarchicalBins::from_schema(&self.schema);
        let bins_to_check = binning.region_to_bins(start, end);

        let mut chunks = Vec::new();

        // Calculate min_offset from all relevant linear index windows
        let min_offset = if self.use_linear_index {
            let start_window = start >> binning.linear_shift;
            let end_window = (end - 1) >> binning.linear_shift;
            let max_window = chrom_index.linear_index.len() as u32;

            let mut min = u64::MAX;
            for window in start_window..=end_window.min(max_window - 1) {
                min = min.min(chrom_index.linear_index[window as usize]);
            }
            min
        } else {
            0 // When linear index is disabled, check all chunks
        };

        // Collect relevant chunks from each bin
        for bin_id in bins_to_check {
            if let Some(bin_chunks) = chrom_index.bins.get(&bin_id) {
                for chunk in bin_chunks {
                    if chunk.end_offset >= min_offset {
                        chunks.push(chunk.clone());
                    }
                }
            }
        }

        // Sort and merge overlapping chunks
        if !chunks.is_empty() {
            chunks.sort_by_key(|c| c.start_offset);
            let mut merged = vec![chunks[0].clone()];

            for chunk in chunks.into_iter().skip(1) {
                let last = merged.last_mut().unwrap();
                if chunk.start_offset <= last.end_offset {
                    // Extend existing chunk
                    last.end_offset = last.end_offset.max(chunk.end_offset);
                } else {
                    // Add new chunk
                    merged.push(chunk);
                }
            }

            // Return the start offsets of chunks
            merged.into_iter().map(|c| c.start_offset).collect()
        } else {
            vec![]
        }
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

    #[test]
    fn test_linear_index_edge_cases() {
        let mut index = BinningIndex::<()>::new();

        // Add features - mix of overlapping and non-overlapping
        index.add_feature("chr1", 16000, 17000, 1); // Should overlap
        index.add_feature("chr1", 15900, 16500, 2); // Should overlap
        index.add_feature("chr1", 16300, 16400, 3); // Should overlap
        index.add_feature("chr1", 15000, 16384, 4); // Should overlap
        index.add_feature("chr1", 16384, 17000, 5); // Should overlap
        index.add_feature("chr1", 15000, 16200, 6); // Should NOT overlap - ends before query
        index.add_feature("chr1", 16500, 17000, 7); // Should NOT overlap - starts after query

        // Query region
        let query_start = 16300;
        let query_end = 16400;

        let mut overlaps_with_linear = index.find_overlapping("chr1", query_start, query_end);
        overlaps_with_linear.sort();

        // Just features 1-5 should overlap
        let mut all_overlapping = vec![1, 2, 3, 4, 5];
        all_overlapping.sort();

        assert_eq!(
            overlaps_with_linear,
            all_overlapping,
            "Expected overlapping features {:?} but got {:?}. Features 6 and 7 should NOT be included.", 
            all_overlapping,
            overlaps_with_linear
        );

        // Also verify same results with linear indexing disabled
        index.disable_linear_index();
        let mut overlaps_without_linear = index.find_overlapping("chr1", query_start, query_end);
        overlaps_without_linear.sort();

        assert_eq!(
            overlaps_with_linear, overlaps_without_linear,
            "Results differ with linear indexing disabled: {:?} vs {:?}",
            overlaps_with_linear, overlaps_without_linear
        );
    }
}
