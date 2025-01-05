use std::{
    collections::HashMap,
    fs::{self, File},
    marker::PhantomData,
    path::{Path, PathBuf},
};

use serde::{Deserialize, Serialize};

use crate::{
    block::{BlockReader, BlockWriter, VirtualOffset},
    error::HgIndexError,
    index::BinningIndex,
    SerdeType,
};

#[derive(Debug)]
pub struct ChromosomeReader<T> {
    reader: BlockReader<T>,
}

impl<T> ChromosomeReader<T>
where
    T: for<'de> Deserialize<'de>,
{
    pub fn new(path: impl AsRef<Path>) -> Result<Self, HgIndexError> {
        let file = File::open(path)?;
        Ok(Self {
            reader: BlockReader::new(file)?,
        })
    }

    pub fn read_at(&mut self, voffset: VirtualOffset) -> Result<T, HgIndexError> {
        self.reader.read_record(voffset)
    }
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

    pub fn write_record<T: Serialize>(
        &mut self,
        record: &T,
    ) -> Result<VirtualOffset, HgIndexError> {
        self.writer.add_record(record)
    }
}

#[derive(Debug)]
pub struct GenomicDataStore<T, M>
where
    T: SerdeType,
    M: SerdeType,
{
    index: BinningIndex<M>,
    writers: HashMap<String, ChromosomeWriter>,
    readers: HashMap<String, ChromosomeReader<T>>,
    directory: PathBuf,
    key: Option<String>,
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
        // Create base directory
        fs::create_dir_all(directory)?;

        // Create key subdirectory if needed
        if let Some(ref key) = key {
            fs::create_dir_all(directory.join(key))?;
        }

        Ok(Self {
            index: BinningIndex::new(),
            writers: HashMap::new(),
            readers: HashMap::new(),
            directory: directory.to_path_buf(),
            key,
            _phantom: PhantomData,
        })
    }

    fn get_or_create_writer(&mut self, chrom: &str) -> Result<&mut ChromosomeWriter, HgIndexError> {
        if !self.writers.contains_key(chrom) {
            let path = self.get_data_path(chrom);
            // Need to create parent directories here!
            std::fs::create_dir_all(path.parent().unwrap())?; // <-- Add this
            let writer = ChromosomeWriter::new(path.to_str().unwrap())?;
            self.writers.insert(chrom.to_string(), writer);
        }
        Ok(self.writers.get_mut(chrom).unwrap())
    }

    fn get_or_create_reader(
        &mut self,
        chrom: &str,
    ) -> Result<&mut ChromosomeReader<T>, HgIndexError> {
        if !self.readers.contains_key(chrom) {
            let path = self.get_data_path(chrom);
            let reader = ChromosomeReader::new(path.to_str().unwrap())?;
            self.readers.insert(chrom.to_string(), reader);
        }
        Ok(self.readers.get_mut(chrom).unwrap())
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

    pub fn add_record(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
        record: &T,
    ) -> Result<(), HgIndexError> {
        let store = self.get_or_create_writer(chrom)?;

        // Add the record to the block-compressed writer for this
        // chromosome.
        let virtual_offset = store.writer.add_record(&record)?;

        // Add the feature to the index now.
        self.index
            .add_feature(chrom, start, end, virtual_offset.into());

        Ok(())
    }

    pub fn finalize(&mut self) -> std::result::Result<(), Box<dyn std::error::Error>> {
        // Finish each writer explicitly before clearing
        for (_, writer) in self.writers.drain() {
            writer.writer.finish()?;
        }

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
        // Keep original path
        let index_path = if let Some(ref key) = key {
            directory.join(key).join(Self::INDEX_FILENAME)
        } else {
            directory.join(Self::INDEX_FILENAME)
        };

        let index = BinningIndex::open(&index_path)?;

        Ok(Self {
            index,
            directory: directory.to_path_buf(), // Use original directory
            writers: HashMap::new(),
            readers: HashMap::new(),
            key,
            _phantom: PhantomData,
        })
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

        let voffsets = self.index.find_overlapping(chrom, start, end);

        // Open chromosome file if not already open, returning empty results
        // if the chromosome is not found (assumed to mean no records for that
        // chromosome).
        match self.get_or_create_reader(chrom) {
            Ok(reader) => {
                for voffset in voffsets {
                    results.push(reader.read_at(voffset.into())?);
                }
            }
            Err(_) => return Ok(results),
        };

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
        let base_dir = test_dir.path(); // Don't add test.gidx

        // An example key
        let key = "example-key".to_string();

        // Create store and add records
        let mut store = GenomicDataStore::<TestRecord, ()>::create(base_dir, Some(key.clone()))
            .expect("Failed to create store");
        for (chrom, start, end, record) in make_test_records() {
            store
                .add_record(&chrom, start, end, &record)
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
                        &MinimalTestRecord { score: i as f64 },
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
