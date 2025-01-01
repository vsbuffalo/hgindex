use std::{
    collections::HashMap,
    fs::{self, File, OpenOptions},
    io::{self, Seek, SeekFrom},
    marker::PhantomData,
    path::{Path, PathBuf},
};

use crate::{
    block::{write_block, BlockConfig, BlockReader, BlockWriter},
    error::HgIndexError,
    index::BinningIndex,
    SerdeType,
};

#[derive(Debug)]
pub struct GenomicDataStore<T, M>
where
    T: SerdeType,
    M: SerdeType,
{
    /// The hierarchical binning index.
    index: BinningIndex<M>,
    /// The directory of per-sequence files.
    directory: PathBuf,
    /// An optional hierarchical "key" for storing multi-level data.
    key: Option<String>,
    /// The block compressed writers for each sequence (e.g. chromosome).
    block_writers: HashMap<String, BlockWriter<T>>,
    _phantom: PhantomData<T>,
}

impl<T: SerdeType, M: SerdeType> GenomicDataStore<T, M> {
    const INDEX_FILENAME: &'static str = "index.bin";

    fn get_data_path(&self, chrom: &str) -> PathBuf {
        let mut path = self.directory.clone();
        if let Some(key) = &self.key {
            path = path.join(key);
        }
        path.join(format!("{}.bin", chrom))
    }

    pub fn create(directory: &Path, key: Option<String>) -> std::io::Result<Self> {
        // Create main directory and key subdirectory if needed
        let _target_dir = if let Some(ref key) = key {
            let key_dir = directory.join(key);
            fs::create_dir_all(&key_dir)?;
            key_dir
        } else {
            fs::create_dir_all(directory)?;
            directory.to_path_buf()
        };

        Ok(Self {
            index: BinningIndex::new(),
            directory: directory.to_path_buf(),
            key,
            block_writers: HashMap::new(),
            _phantom: PhantomData,
        })
    }

    /// Get a reference to the index's metadata.
    pub fn metadata(&self) -> Option<&M> {
        self.index.metadata.as_ref()
    }

    /// Get a mutable reference to the store's metadata
    pub fn metadata_mut(&mut self) -> Option<&mut M> {
        self.index.metadata.as_mut()
    }

    /// Set the store's metadata
    pub fn set_metadata(&mut self, metadata: M) {
        self.index.metadata = Some(metadata);
    }

    /// Take ownership of the metadata, leaving None in its place
    pub fn take_metadata(&mut self) -> Option<M> {
        self.index.metadata.take()
    }

    pub fn add_record(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
        record: T,
    ) -> Result<(), HgIndexError> {
        // println!("Adding record: {}:{}-{}", chrom, start, end);

        let path = self.get_data_path(chrom);
        let mut file = OpenOptions::new().create(true).append(true).open(&path)?;

        // Get or create BlockWriter for this chromosome
        let writer = self
            .block_writers
            .entry(chrom.to_string())
            .or_insert_with(|| {
                BlockWriter::new(BlockConfig::default()).expect("Failed to create block writer")
            });

        // Add to block writer
        if let Some(block) = writer.add_record(start, end, record)? {
            let offset = file.stream_position()?;
            // println!( "Created block at offset {} with range {}-{}", offset, block.start, block.end);

            write_block(&mut file, &block)?;
            self.index.add_feature(chrom, start, end, offset)?;
        }

        Ok(())
    }

    pub fn finalize(&mut self) -> Result<(), HgIndexError> {
        // First collect all paths to avoid borrow issues
        let mut path_map: HashMap<String, PathBuf> = HashMap::new();
        for chrom in self.block_writers.keys() {
            path_map.insert(chrom.clone(), self.get_data_path(chrom));
        }

        // Now handle the blocks
        for (chrom, writer) in self.block_writers.iter_mut() {
            // Force flush any remaining records
            if let Some(block) = writer.finish()? {
                // Call finish() to get remaining records
                let path = &path_map[chrom];
                let mut file = OpenOptions::new().create(true).append(true).open(path)?;
                let offset = file.stream_position()?;
                // println!( "In finalize: Writing block for {} at offset {} with range {}-{}", chrom, offset, block.start, block.end);
                write_block(&mut file, &block)?;
                self.index
                    .add_feature(chrom, block.start, block.end, offset)?;
            }
        }

        // Write index to file
        let index_path = if let Some(ref key) = self.key {
            self.directory.join(key).join(Self::INDEX_FILENAME)
        } else {
            self.directory.join(Self::INDEX_FILENAME)
        };

        self.index.write(index_path.as_path())?;
        Ok(())
    }

    pub fn open(
        directory: &Path,
        key: Option<String>,
    ) -> std::result::Result<Self, Box<dyn std::error::Error>> {
        let target_dir = if let Some(ref key) = key {
            directory.join(key)
        } else {
            directory.to_path_buf()
        };

        // Read index file
        let index_path = target_dir.join(Self::INDEX_FILENAME);
        let index = BinningIndex::open(&index_path)?;

        Ok(Self {
            index,
            directory: directory.to_path_buf(),
            key,
            block_writers: HashMap::new(), // Empty HashMap for reading
            _phantom: PhantomData,
        })
    }

    pub fn get_overlapping(&self, chrom: &str, start: u32, end: u32) -> io::Result<Vec<T>> {
        let mut results = Vec::new();

        // println!("Querying: {}:{}-{}", chrom, start, end);

        // Early return if chromosome not in index
        if !self.index.sequences.contains_key(chrom) {
            // println!("Chromosome {} not found in index", chrom);
            return Ok(results);
        }

        // Open the file for reading
        let path = self.get_data_path(chrom);
        let file = File::open(path)?;
        let mut reader = BlockReader::new(file);

        let offsets = self.index.find_overlapping(chrom, start, end);
        // println!("Found {} potential blocks to check", offsets.len());

        for offset in offsets {
            reader.seek(SeekFrom::Start(offset))?;
            let block = reader.read_header()?;
            // println!("Checking block: {}-{}", block.start, block.end);

            let records = reader.read_records(&block)?;
            for (rec_start, rec_end, record) in records {
                // println!("Checking record: {}-{}", rec_start, rec_end);
                if rec_start < end && rec_end >= start {
                    results.push(record);
                }
            }
        }

        Ok(results)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::test_utils::TestDir;
    use serde::{Deserialize, Serialize};

    // A simple test record type
    #[derive(Debug, Serialize, Deserialize, PartialEq)]
    struct TestRecord {
        name: String,
        score: f64,
        tags: Vec<String>,
    }

    fn make_test_records() -> Vec<(String, u32, u32, TestRecord)> {
        vec![
            (
                "chr1".to_string(),
                1000,
                2000,
                TestRecord {
                    name: "feature1".to_string(),
                    score: 0.5,
                    tags: vec!["exon".to_string(), "coding".to_string()],
                },
            ),
            (
                "chr1".to_string(),
                1500,
                2500,
                TestRecord {
                    name: "feature2".to_string(),
                    score: 0.8,
                    tags: vec!["promoter".to_string()],
                },
            ),
            (
                "chr2".to_string(),
                50000,
                60000,
                TestRecord {
                    name: "feature3".to_string(),
                    score: 0.3,
                    tags: vec!["intron".to_string()],
                },
            ),
        ]
    }

    #[test]
    fn test_store_and_retrieve() {
        let test_dir = TestDir::new("store_and_retrieve").expect("Failed to create test dir");
        let store_path = test_dir.path().join("test.gidx");

        // An example key
        let key = "example-key".to_string();

        // Create store and add records
        let mut store = GenomicDataStore::<TestRecord, ()>::create(&store_path, Some(key.clone()))
            .expect("Failed to create store");

        for (chrom, start, end, record) in make_test_records() {
            store
                .add_record(&chrom, start, end, record)
                .expect("Failed to add record");
        }

        store.finalize().expect("Failed to finalize store");

        // Open store and query
        let store = GenomicDataStore::<TestRecord, ()>::open(&store_path, Some(key))
            .expect("Failed to open store");

        // Test overlapping query
        let results = store
            .get_overlapping("chr1", 1200, 1800)
            .expect("Failed to get overlapping records");
        assert_eq!(results.len(), 2); // Should get both chr1 features
        assert_eq!(results[0].name, "feature1");
        assert_eq!(results[1].name, "feature2");

        // Test non-overlapping region
        let results = store
            .get_overlapping("chr1", 3000, 4000)
            .expect("Failed to get overlapping records");
        assert_eq!(results.len(), 0);

        // Test different chromosome
        let results = store
            .get_overlapping("chr2", 55000, 58000)
            .expect("Failed to get overlapping records");
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].name, "feature3");
    }

    #[test]
    fn test_empty_regions() {
        let test_dir = TestDir::new("empty_regions").expect("Failed to create test dir");
        let store_path = test_dir.path().join("empty.gidx");

        let mut store = GenomicDataStore::<TestRecord, ()>::create(&store_path, None)
            .expect("Failed to create store");

        store.finalize().expect("Failed to finalize store");

        // Query empty store
        let store = GenomicDataStore::<TestRecord, ()>::open(&store_path, None)
            .expect("Failed to open store");

        let results = store
            .get_overlapping("chr1", 0, 1000)
            .expect("Failed to get overlapping records");
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_metadata_handling() {
        let test_dir = TestDir::new("metadata").expect("Failed to create test dir");
        let store_path = test_dir.path().join("meta.gidx");

        // Create a store with initial metadata
        let mut store = GenomicDataStore::<TestRecord, String>::create(&store_path, None)
            .expect("Failed to create store");

        store.set_metadata("initial metadata".to_string());

        // Test metadata access
        assert_eq!(store.metadata(), Some(&"initial metadata".to_string()));

        // Test metadata mutation
        if let Some(meta) = store.metadata_mut() {
            meta.push_str(" - updated");
        }
        assert_eq!(
            store.metadata(),
            Some(&"initial metadata - updated".to_string())
        );

        // Test metadata replacement
        store.set_metadata("new metadata".to_string());
        assert_eq!(store.metadata(), Some(&"new metadata".to_string()));

        // Test metadata removal
        let taken = store.take_metadata();
        assert_eq!(taken, Some("new metadata".to_string()));
        assert_eq!(store.metadata(), None);

        store.finalize().expect("Failed to finalize store");

        // Test metadata persistence
        let reopened = GenomicDataStore::<TestRecord, String>::open(&store_path, None)
            .expect("Failed to open store");
        assert_eq!(reopened.metadata(), None);
    }

    // A minimal test record
    #[derive(Debug, Serialize, Deserialize, PartialEq)]
    struct MinimalTestRecord {
        score: f64,
    }

    #[test]
    fn test_concurrent_reads() {
        use std::sync::Arc;
        use std::thread;

        let test_dir = TestDir::new("concurrent").expect("Failed to create test dir");
        let store_path = test_dir.path().join("test.gidx");

        // Create and populate store
        {
            let mut store = GenomicDataStore::<MinimalTestRecord, ()>::create(&store_path, None)
                .expect("Failed to create store");

            // Add some overlapping records
            for i in 0..10 {
                store
                    .add_record(
                        "chr1",
                        i * 1000,
                        (i + 2) * 1000, // Overlapping regions
                        MinimalTestRecord { score: i as f64 },
                    )
                    .expect("Failed to add record");
            }
            store.finalize().expect("Failed to finalize");
        }

        // Create path that can be shared between threads
        let path = Arc::new(store_path);

        // Spawn multiple reader threads
        let handles: Vec<_> = (0..4)
            .map(|i| {
                let path = Arc::clone(&path);
                thread::spawn(move || {
                    let store = GenomicDataStore::<MinimalTestRecord, ()>::open(&path, None)
                        .expect("Failed to open store");

                    // Each thread queries a different but overlapping region
                    let start = i * 500;
                    let end = start + 2000;
                    let results = store
                        .get_overlapping("chr1", start, end)
                        .expect("Failed to get overlapping records");

                    // Results should not be empty due to overlapping regions
                    assert!(!results.is_empty());
                    results.len()
                })
            })
            .collect();

        // Verify all threads completed successfully
        let result_counts: Vec<_> = handles
            .into_iter()
            .map(|h| h.join().expect("Thread panicked"))
            .collect();

        // Verify that at least some threads got different numbers of results
        // due to querying different regions
        assert!(result_counts.iter().any(|&x| x != result_counts[0]));
    }

    #[test]
    fn test_edge_case_overlaps() {
        let test_dir = TestDir::new("edge_overlaps").expect("Failed to create test dir");
        let store_path = test_dir.path().join("test.gidx");

        let mut store = GenomicDataStore::<MinimalTestRecord, ()>::create(&store_path, None)
            .expect("Failed to create store");

        // Test cases with non-round numbers
        let test_cases = vec![
            // Small intervals
            ("chr1", 5010061, 5010100, 1.0),
            ("chr1", 5010062, 5010099, 2.0),
            // Cross block boundary cases
            ("chr1", 65535, 65537, 3.0),
            ("chr1", 131071, 131073, 4.0),
            // Edge of bin cases
            ("chr1", 4095, 4097, 5.0),
            ("chr1", 32767, 32769, 6.0),
            // Odd sized intervals
            ("chr1", 1234567, 1234601, 7.0),
            ("chr1", 7654321, 7654399, 8.0),
        ];

        // Add all records
        for (chrom, start, end, score) in &test_cases {
            store
                .add_record(chrom, *start, *end, MinimalTestRecord { score: *score })
                .expect("Failed to add record");
        }

        store.finalize().expect("Failed to finalize store");

        // Open store for querying
        let store = GenomicDataStore::<MinimalTestRecord, ()>::open(&store_path, None)
            .expect("Failed to open store");

        // Test various overlap scenarios
        let overlap_tests = vec![
            // Exact matches with overlapping ranges - expect both matches
            ("chr1", 5010061, 5010100, 2), // Changed from 1 to 2
            // Partial overlaps
            ("chr1", 5010050, 5010080, 2), // Should find both overlapping records
            ("chr1", 5010080, 5010120, 2), // Should find both overlapping records
            // Multiple overlaps
            ("chr1", 5010060, 5010101, 2), // Should find both overlapping records
            // Block boundary overlaps
            ("chr1", 65530, 65540, 1),
            ("chr1", 65530, 65600, 1),
            // Bin boundary overlaps
            ("chr1", 4090, 4100, 1),
            ("chr1", 32760, 32780, 1),
        ];

        for (chrom, start, end, expected_count) in overlap_tests {
            let results = store
                .get_overlapping(chrom, start, end)
                .expect("Failed to get overlapping records");
            assert_eq!(
                results.len(),
                expected_count,
                "Failed for range {}:{}-{}, expected {} overlaps, got {}",
                chrom,
                start,
                end,
                expected_count,
                results.len()
            );
        }
    }

    #[test]
    fn test_block_boundary_handling() {
        let test_dir = TestDir::new("block_boundaries").expect("Failed to create test dir");
        let store_path = test_dir.path().join("test.gidx");
        let mut store = GenomicDataStore::<MinimalTestRecord, ()>::create(&store_path, None)
            .expect("Failed to create store");

        // Create features that cross block boundaries
        // Each feature will span 200 bases centered on a 64KB boundary
        let boundaries = [65536, 131072, 196608]; // 64KB, 128KB, 192KB boundaries
        for (i, boundary) in boundaries.iter().enumerate() {
            let start = boundary - 100; // 100 bases before boundary
            let end = boundary + 100; // 100 bases after boundary
            store
                .add_record("chr1", start, end, MinimalTestRecord { score: i as f64 })
                .expect("Failed to add record");
        }

        store.finalize().expect("Failed to finalize store");

        let store = GenomicDataStore::<MinimalTestRecord, ()>::open(&store_path, None)
            .expect("Failed to open store");

        // Test different types of queries across boundaries
        let test_cases = vec![
            // Query exactly at boundary
            (65536 - 1, 65536 + 1, 1), // Should find the record crossing 64KB
            // Query fully containing a boundary-crossing feature
            (65436, 65636, 1), // Should find the 64KB crossing record
            // Query overlapping multiple boundary regions
            (65436, 131636, 2), // Should find both 64KB and 128KB records
            // Query between boundaries
            (100000, 100100, 0), // Should find no records
        ];

        for (query_start, query_end, expected_count) in test_cases {
            let results = store
                .get_overlapping("chr1", query_start, query_end)
                .expect("Failed to get overlapping records");
            assert_eq!(
                results.len(),
                expected_count,
                "Query {}..{} found {} overlaps, expected {}",
                query_start,
                query_end,
                results.len(),
                expected_count,
            );
        }
    }

    #[test]
    fn test_tiny_overlaps() {
        let test_dir = TestDir::new("tiny_overlaps").expect("Failed to create test dir");
        let store_path = test_dir.path().join("test.gidx");
        let mut store = GenomicDataStore::<MinimalTestRecord, ()>::create(&store_path, None)
            .expect("Failed to create store");

        // Add records with 1-base overlaps at different positions
        let test_cases = vec![
            (1000, 2000, 1999, 2001), // 1-base overlap at end
            (3000, 4000, 2999, 3001), // 1-base overlap at start
            (5000, 5001, 5000, 5001), // Exact 1-base overlap
        ];

        for (start1, end1, start2, end2) in &test_cases {
            store
                .add_record("chr1", *start1, *end1, MinimalTestRecord { score: 1.0 })
                .expect("Failed to add first record");
            store
                .add_record("chr1", *start2, *end2, MinimalTestRecord { score: 2.0 })
                .expect("Failed to add second record");
        }

        store.finalize().expect("Failed to finalize store");

        let store = GenomicDataStore::<MinimalTestRecord, ()>::open(&store_path, None)
            .expect("Failed to open store");

        // Test each overlap case
        for (start1, end1, start2, end2) in test_cases {
            let results = store
                .get_overlapping("chr1", start1, end1)
                .expect("Failed to get overlapping records");
            assert_eq!(
                results.len(),
                2,
                "Failed to find tiny overlap at {}-{}",
                start1,
                end1
            );

            let results = store
                .get_overlapping("chr1", start2, end2)
                .expect("Failed to get overlapping records");
            assert_eq!(
                results.len(),
                2,
                "Failed to find tiny overlap at {}-{}",
                start2,
                end2
            );
        }
    }
    #[test]
    fn test_block_file_offsets() {
        let test_dir = TestDir::new("block_offsets").expect("Failed to create test dir");
        let store_path = test_dir.path().join("test.gidx");
        let mut store = GenomicDataStore::<MinimalTestRecord, ()>::create(&store_path, None)
            .expect("Failed to create store");

        // Add records with distinctive scores so we can track them
        store
            .add_record("chr1", 1000, 2000, MinimalTestRecord { score: 1.0 })
            .unwrap();

        store
            .add_record("chr1", 5000000, 5000100, MinimalTestRecord { score: 2.0 })
            .unwrap();

        store
            .add_record("chr1", 3000, 4000, MinimalTestRecord { score: 3.0 })
            .unwrap();

        store.finalize().expect("Failed to finalize store");

        // Open store for querying
        let store = GenomicDataStore::<MinimalTestRecord, ()>::open(&store_path, None)
            .expect("Failed to open store");

        // Add assertions
        assert!(store.get_overlapping("chr1", 1000, 2000).unwrap().len() > 0);
        assert!(
            store
                .get_overlapping("chr1", 5000000, 5000100)
                .unwrap()
                .len()
                > 0
        );
        assert!(store.get_overlapping("chr1", 3000, 4000).unwrap().len() > 0);
    }
}
