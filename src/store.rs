use std::io;
use std::{
    collections::HashMap,
    fs::{self, File},
    io::{BufWriter, Read, Seek, SeekFrom, Write},
    marker::PhantomData,
    path::{Path, PathBuf},
};

use crate::GenomicCoordinates;
use crate::{error::HgIndexError, index::BinningIndex, BinningSchema, SerdeType};

#[derive(Debug)]
pub struct GenomicDataStore<T, M>
where
    T: GenomicCoordinates + SerdeType,
    M: SerdeType + std::fmt::Debug,
{
    index: BinningIndex<M>,
    data_files: HashMap<String, File>,
    directory: PathBuf,
    key: Option<String>,
    _phantom: PhantomData<T>,
}

impl<T: GenomicCoordinates + SerdeType, M: SerdeType + std::fmt::Debug> GenomicDataStore<T, M> {
    const MAGIC: [u8; 4] = *b"GIDX";
    const INDEX_FILENAME: &'static str = "index.bin";

    fn get_data_path(&self, chrom: &str) -> PathBuf {
        let mut path = self.directory.clone();
        if let Some(key) = &self.key {
            path = path.join(key);
        }
        path.join(format!("{}.bin", chrom))
    }

    pub fn create(directory: &Path, key: Option<String>) -> io::Result<Self> {
        Self::create_with_schema(directory, key, &BinningSchema::default())
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
            index: BinningIndex::from_schema(schema),
            data_files: HashMap::new(),
            directory: directory.to_path_buf(),
            key,
            _phantom: PhantomData,
        })
    }

    fn get_or_create_file(&mut self, chrom: &str) -> std::io::Result<&mut File> {
        if !self.data_files.contains_key(chrom) {
            let data_path = self.get_data_path(chrom);
            let file = File::create(&data_path)?;
            let mut writer = BufWriter::new(&file);
            writer.write_all(&Self::MAGIC)?;
            writer.flush()?;
            let _data_file = writer.into_inner()?;
            self.data_files.insert(chrom.to_string(), file);
        }
        Ok(self.data_files.get_mut(chrom).unwrap())
    }

    // Add these methods to your existing impl block
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

    pub fn add_record(&mut self, chrom: &str, record: &T) -> std::io::Result<()> {
        let file = self.get_or_create_file(chrom)?;

        // Create a scope to ensure writer is dropped before we use self.index
        let offset = {
            let mut writer = BufWriter::new(file);
            let offset = writer.stream_position()?;

            // Serialize record and get its length
            let record_data = bincode::serialize(record).unwrap();
            let length = record_data.len() as u64;

            // Write length followed by record
            writer.write_all(&length.to_le_bytes())?;
            writer.write_all(&record_data)?;
            writer.flush()?;

            offset
        }; // writer is dropped here, releasing the borrow

        // Add the feature to the index now.
        let start = record.start();
        let end = record.end();
        self.index.add_feature(chrom, start, end, offset);

        Ok(())
    }

    pub fn finalize(&mut self) -> std::result::Result<(), Box<dyn std::error::Error>> {
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

        // Initialize without opening any chromosome files yet
        Ok(Self {
            index,
            data_files: HashMap::new(),
            directory: directory.to_path_buf(),
            key,
            _phantom: PhantomData,
        })
    }

    pub fn open_chrom_file(&mut self, chrom: &str) -> std::io::Result<()> {
        if !self.data_files.contains_key(chrom) {
            let data_path = self.get_data_path(chrom);
            let mut file = File::open(&data_path)?;

            // Verify magic number
            let mut magic = [0u8; 4];
            file.read_exact(&mut magic)?;
            if magic != Self::MAGIC {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "Invalid file format",
                ));
            }

            self.data_files.insert(chrom.to_string(), file);
        }
        Ok(())
    }

    /// The features overlapping the range with this start and end position.
    pub fn get_overlapping(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<Vec<T>, HgIndexError> {
        if end <= start {
            return Err(HgIndexError::InvalidInterval { start, end });
        }
        // Early return empty vec if chromosome not in index
        if !self.index.sequences.contains_key(chrom) {
            return Ok(Vec::new());
        }
        let mut results = Vec::new();

        // Early return if chromosome not in index
        if !self.index.sequences.contains_key(chrom) {
            return Ok(results);
        }

        // Open chromosome file if not already open
        if self.open_chrom_file(chrom).is_err() {
            return Ok(results);
        }

        let file = self.data_files.get_mut(chrom).unwrap();

        for offset in self.index.find_overlapping(chrom, start, end) {
            file.seek(SeekFrom::Start(offset)).unwrap();

            // First read the length (8 bytes)
            let mut length_bytes = [0u8; 8];
            if file.read_exact(&mut length_bytes).is_ok() {
                let length = u64::from_le_bytes(length_bytes);

                // Now read the actual record
                let mut record_data = vec![0u8; length as usize];
                if file.read_exact(&mut record_data).is_ok() {
                    if let Ok(record) = bincode::deserialize(&record_data) {
                        results.push(record);
                    }
                }
            }
        }
        Ok(results)
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
        store.finalize().unwrap();
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
}
