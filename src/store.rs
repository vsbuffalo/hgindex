// store.rs

use std::{
    collections::HashMap,
    fs::{self, File},
    io,
    marker::PhantomData,
    path::{Path, PathBuf},
};

use serde::{Deserialize, Serialize};

use crate::{
    block::{BlockReader, BlockWriter, VirtualOffset},
    error::HgIndexError,
    index::{BinningIndex, BinningSchema},
    GenomicCoordinates, SerdeType,
};

#[derive(Debug)]
pub struct ChromosomeReader<T>
where
    T: GenomicCoordinates + std::fmt::Debug,
{
    reader: BlockReader<T>,
}

impl<T> ChromosomeReader<T>
where
    T: GenomicCoordinates + for<'de> Deserialize<'de> + std::fmt::Debug,
{
    pub fn new(path: impl AsRef<Path>) -> Result<Self, HgIndexError> {
        let file = File::open(path)?;
        Ok(Self {
            reader: BlockReader::new(file)?,
        })
    }

    //pub fn read_at(&mut self, voffset: VirtualOffset) -> Result<T, HgIndexError> {
    //    self.reader.read_record(voffset)
    //}
}

#[derive(Debug)]
pub struct ChromosomeWriter {
    writer: BlockWriter,
}

impl ChromosomeWriter {
    pub fn new(path: impl AsRef<Path>) -> Result<Self, HgIndexError> {
        let file = File::create(path)?;
        Ok(Self {
            writer: BlockWriter::new(file)?,
        })
    }

    pub fn add_record<T: Serialize + std::fmt::Debug>(
        &mut self,
        record: &T,
    ) -> Result<VirtualOffset, HgIndexError> {
        self.writer.add_record(&record)
    }
}

pub struct GenomicDataStore<T, M>
where
    T: GenomicCoordinates + SerdeType + std::fmt::Debug,
    M: SerdeType + std::fmt::Debug,
{
    directory: PathBuf,
    index: Option<BinningIndex<M>>,
    writers: HashMap<String, ChromosomeWriter>,
    readers: HashMap<String, ChromosomeReader<T>>,
    schema: BinningSchema,
    key: Option<String>,
    // For sort validation
    last_chrom: Option<String>,
    last_pos: Option<u32>,
    // Store coordinates and offsets during writing
    coordinates: HashMap<String, Vec<(u32, u32, VirtualOffset)>>,
    _phantom: PhantomData<T>,
}

impl<T, M> GenomicDataStore<T, M>
where
    T: GenomicCoordinates + SerdeType + std::fmt::Debug,
    M: SerdeType + std::fmt::Debug,
{
    pub fn create(directory: &Path, key: Option<String>) -> io::Result<Self> {
        Self::create_with_schema(directory, key, &BinningSchema::default())
    }

    pub fn open(directory: &Path, key: Option<String>) -> Result<Self, HgIndexError> {
        // Load existing index
        let index_path = if let Some(ref key) = key {
            directory.join(key).join("index.bin")
        } else {
            directory.join("index.bin")
        };

        let index = BinningIndex::open(&index_path)?;
        let schema = index.schema.clone();

        Ok(Self {
            directory: directory.to_path_buf(),
            index: Some(index),
            writers: HashMap::new(),
            readers: HashMap::new(),
            schema,
            key,
            last_chrom: None,
            last_pos: None,
            coordinates: HashMap::new(),
            _phantom: PhantomData,
        })
    }

    // Metadata methods
    pub fn metadata(&self) -> Option<&M> {
        self.index.as_ref().and_then(|idx| idx.metadata.as_ref())
    }

    pub fn metadata_mut(&mut self) -> Option<&mut M> {
        self.index.as_mut().and_then(|idx| idx.metadata.as_mut())
    }

    fn set_metadata(&mut self, metadata: M) {
        if let Some(idx) = self.index.as_mut() {
            idx.metadata = Some(metadata);
        } else {
            panic!("Cannot set index-level metadata before index has been generated.");
        }
    }

    pub fn take_metadata(&mut self) -> Option<M> {
        self.index.as_mut().and_then(|idx| idx.metadata.take())
    }

    pub fn create_with_schema(
        directory: &Path,
        key: Option<String>,
        schema: &BinningSchema,
    ) -> io::Result<Self> {
        fs::create_dir_all(directory)?;
        if let Some(ref key) = key {
            fs::create_dir_all(directory.join(key))?;
        }

        Ok(Self {
            directory: directory.to_path_buf(),
            index: None,
            writers: HashMap::new(),
            readers: HashMap::new(),
            schema: schema.clone(),
            key,
            last_chrom: None,
            last_pos: None,
            coordinates: HashMap::new(),
            _phantom: PhantomData,
        })
    }

    fn get_data_path(&self, chrom: &str) -> PathBuf {
        let mut path = self.directory.clone();
        if let Some(key) = &self.key {
            path = path.join(key);
        }
        path.join(format!("{}.bin", chrom))
    }

    fn get_or_create_writer(&mut self, chrom: &str) -> Result<&mut ChromosomeWriter, HgIndexError> {
        if !self.writers.contains_key(chrom) {
            let path = self.get_data_path(chrom);
            std::fs::create_dir_all(path.parent().unwrap())?;
            let writer = ChromosomeWriter::new(path)?;
            self.writers.insert(chrom.to_string(), writer);
        }
        Ok(self.writers.get_mut(chrom).unwrap())
    }

    pub fn add_record(&mut self, chrom: &str, record: &T) -> Result<(), HgIndexError> {
        let start = record.start();
        let end = record.end();

        // Validate coordinates
        if end <= start {
            return Err(HgIndexError::ZeroLengthFeature(start, end));
        }

        // Validate sorting
        if let (Some(last_chrom), Some(last_pos)) = (self.last_chrom.as_ref(), self.last_pos) {
            match last_chrom.as_str().cmp(chrom) {
                std::cmp::Ordering::Greater => {
                    return Err(HgIndexError::from("Records must be sorted by chromosome"));
                }
                std::cmp::Ordering::Equal if start < last_pos => {
                    return Err(HgIndexError::from(
                        "Records must be sorted by position within chromosome",
                    ));
                }
                _ => {}
            }
        }

        // Update last position
        self.last_chrom = Some(chrom.to_string());
        self.last_pos = Some(start);

        // Write record and get its offset
        let writer = self.get_or_create_writer(chrom)?;
        let offset = writer.add_record(record)?;

        // Store coordinates and offset for indexing phase
        self.coordinates
            .entry(chrom.to_string())
            .or_default()
            .push((start, end, offset));

        Ok(())
    }

    fn index(&mut self) -> Result<(), HgIndexError> {
        let mut index = BinningIndex::from_schema(&self.schema);

        // Process each chromosome's coordinates
        for (chrom, coords) in &self.coordinates {
            for (start, end, offset) in coords {
                index.add_feature(chrom, *start, *end, (*offset).into());
            }
        }

        self.index = Some(index);
        Ok(())
    }

    pub fn finalize(&mut self) -> Result<(), HgIndexError> {
        // First flush all writers
        for writer in self.writers.values_mut() {
            writer.writer.flush_block()?;
        }

        // Clear writers
        self.writers.clear();

        // Build and save the index
        // let start = Instant::now();
        self.index()?;
        // let duration = start.elapsed();
        // eprintln!("  Building index took {:?}", duration);

        // Write index to disk
        if let Some(index) = &mut self.index {
            let index_path = if let Some(ref key) = self.key {
                self.directory.join(key).join("index.bin")
            } else {
                self.directory.join("index.bin")
            };
            index.write(&index_path)?;
        }

        Ok(())
    }

    pub fn finalize_with_metadata(&mut self, metadata: M) -> Result<(), HgIndexError> {
        // First flush all writers
        for writer in self.writers.values_mut() {
            writer.writer.flush_block()?;
        }

        // Clear writers
        self.writers.clear();

        // Build and save the index
        // let start = Instant::now();
        self.index()?;
        // let duration = start.elapsed();
        // eprintln!("  Building index took {:?}", duration);

        // Now we can set metadata.
        self.set_metadata(metadata);

        // Write index to disk
        if let Some(index) = &mut self.index {
            let index_path = if let Some(ref key) = self.key {
                self.directory.join(key).join("index.bin")
            } else {
                self.directory.join("index.bin")
            };
            index.write(&index_path)?;
        }

        Ok(())
    }

    fn get_or_create_reader(
        &mut self,
        chrom: &str,
    ) -> Result<&mut ChromosomeReader<T>, HgIndexError> {
        if !self.readers.contains_key(chrom) {
            let path = self.get_data_path(chrom);
            let reader = ChromosomeReader::new(path)?;
            self.readers.insert(chrom.to_string(), reader);
        }
        Ok(self.readers.get_mut(chrom).unwrap())
    }

    pub fn get_overlapping(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<Vec<T>, HgIndexError> {
        if end <= start {
            return Err(HgIndexError::InvalidInterval { start, end });
        }

        let index = self
            .index
            .as_ref()
            .ok_or_else(|| HgIndexError::from("Index not built yet"))?;

        // Get the range of offsets that might contain overlapping features
        let Some((min_offset, max_offset)) = index.get_candidate_offsets(chrom, start, end) else {
            return Ok(Vec::new());
        };

        let reader = self.get_or_create_reader(chrom)?;

        // Convert to VirtualOffset and read records
        reader.reader.read_records_between(
            VirtualOffset::from(min_offset),
            VirtualOffset::from(max_offset),
            start,
            end,
        )
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use crate::test_utils::test_utils::TestDir;

    use super::*;
    use serde::{Deserialize, Serialize};

    // A simple test record type
    #[derive(Debug, Serialize, Deserialize, PartialEq)]
    struct TestRecord {
        start: u32,
        end: u32,
        name: String,
        score: f64,
        tags: Vec<String>,
    }

    impl GenomicCoordinates for TestRecord {
        fn start(&self) -> u32 {
            self.start
        }

        fn end(&self) -> u32 {
            self.end
        }
    }

    fn make_test_records() -> Vec<(String, TestRecord)> {
        vec![
            (
                "chr1".to_string(),
                TestRecord {
                    start: 1000,
                    end: 2000,
                    name: "feature1".to_string(),
                    score: 0.5,
                    tags: vec!["exon".to_string(), "coding".to_string()],
                },
            ),
            (
                "chr1".to_string(),
                TestRecord {
                    start: 1500,
                    end: 2500,
                    name: "feature2".to_string(),
                    score: 0.8,
                    tags: vec!["promoter".to_string()],
                },
            ),
            (
                "chr2".to_string(),
                TestRecord {
                    start: 50000,
                    end: 60000,
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
        let base_dir = test_dir.path(); // Don't add test.gidx

        // An example key
        let key = "example-key".to_string();

        // Create store and add records
        let mut store = GenomicDataStore::<TestRecord, ()>::create(base_dir, Some(key.clone()))
            .expect("Failed to create store");
        for (chrom, record) in make_test_records() {
            store
                .add_record(&chrom, &record)
                .expect("Failed to add record");
        }

        store.finalize().expect("Failed to finalize store");

        let mut store = GenomicDataStore::<TestRecord, ()>::open(&base_dir, Some(key.clone()))
            .expect("Failed to open store");

        // Test overlapping query
        let results = store.get_overlapping("chr1", 1200, 1800).unwrap();
        assert_eq!(results.len(), 2); // Should get both chr1 features
        assert_eq!(results[0].name, "feature1");
        assert_eq!(results[1].name, "feature2");

        // Test non-overlapping region
        let results = store.get_overlapping("chr1", 3000, 4000).unwrap();
        assert_eq!(results.len(), 0);

        // Test different chromosome
        let results = store.get_overlapping("chr2", 55000, 58000).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].name, "feature3");
    }

    #[test]
    fn test_invalid_file() {
        let test_dir = TestDir::new("invalid_file").expect("Failed to create test dir");
        let bad_file = test_dir.path().join("bad.gidx");

        // Create file with invalid magic number
        let mut file = File::create(&bad_file).expect("Failed to create file");
        file.write_all(b"BAD!").expect("Failed to write");

        // Attempt to open should fail
        let result = GenomicDataStore::<TestRecord, ()>::open(&bad_file, None);
        assert!(result.is_err());
    }

    #[test]
    fn test_empty_regions() {
        let test_dir = TestDir::new("empty_regions").expect("Failed to create test dir");
        let store_path = test_dir.path().join("empty.gidx");

        let mut store = GenomicDataStore::<TestRecord, ()>::create(&store_path, None)
            .expect("Failed to create store");

        store.finalize().expect("Failed to finalize store");

        // Query empty store
        let mut store = GenomicDataStore::<TestRecord, ()>::open(&store_path, None)
            .expect("Failed to open store");

        let results = store.get_overlapping("chr1", 0, 1000).unwrap();
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_metadata_handling() {
        let test_dir = TestDir::new("metadata").expect("Failed to create test dir");
        let store_path = test_dir.path().join("meta.gidx");

        // Create a store with initial metadata
        let mut store = GenomicDataStore::<TestRecord, String>::create(&store_path, None)
            .expect("Failed to create store");

        // Index must be created before metadata.
        store.index().unwrap();
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
        start: u32,
        end: u32,
        score: f64,
    }

    impl GenomicCoordinates for MinimalTestRecord {
        fn start(&self) -> u32 {
            self.start
        }

        fn end(&self) -> u32 {
            self.end
        }
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
                let start = i * 1000;
                let end = (i + 2) * 1000; // Overlapping regions
                store
                    .add_record(
                        "chr1",
                        &MinimalTestRecord {
                            start,
                            end,
                            score: i as f64,
                        },
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
                    let mut store = GenomicDataStore::<MinimalTestRecord, ()>::open(&path, None)
                        .expect("Failed to open store");

                    // Each thread queries a different but overlapping region
                    let start = i * 500;
                    let end = start + 2000;
                    let results = store.get_overlapping("chr1", start, end).unwrap();

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
    fn test_zero_length_features() {
        let tmp_dir = tempfile::tempdir().expect("Failed to create temp dir");
        let store_path = tmp_dir.path();

        let mut store = GenomicDataStore::<TestRecord, ()>::create(store_path, None)
            .expect("Failed to create store");

        // Create a test record to use
        fn make_test_record(start: u32, end: u32) -> TestRecord {
            TestRecord {
                start,
                end,
                name: "test_feature".to_string(),
                score: 42.0,
                tags: vec!["test_tag".to_string(), "zero_length".to_string()],
            }
        }

        // Test case 1: start equals end
        let result = store.add_record("chr1", &make_test_record(100, 100));
        assert!(matches!(
            result,
            Err(HgIndexError::ZeroLengthFeature(100, 100))
        ));

        // Test case 2: start greater than end
        let result = store.add_record("chr1", &make_test_record(200, 100));
        assert!(matches!(
            result,
            Err(HgIndexError::ZeroLengthFeature(200, 100))
        ));

        // Test case 3: valid record should work
        let result = store.add_record("chr1", &make_test_record(300, 310));
        assert!(result.is_ok());

        // Test case 4: minimum valid length (1)
        let result = store.add_record("chr1", &make_test_record(300, 301));
        assert!(result.is_ok());

        // Verify we can retrieve the valid records
        if let Ok(records) = store.get_overlapping("chr1", 100, 200) {
            assert_eq!(records.len(), 1);
            assert_eq!(records[0], make_test_record(100, 107));
        }
    }
}
