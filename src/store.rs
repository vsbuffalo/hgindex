// store.rs

use std::io;
use std::{
    collections::HashMap,
    fs::{self, File},
    io::{BufWriter, Seek, Write},
    marker::PhantomData,
    path::{Path, PathBuf},
};

use memmap2::Mmap;
use serde::{Deserialize, Serialize};

use crate::{error::HgIndexError, index::BinningIndex, BinningSchema};
use crate::{Record, RecordSlice};

#[derive(Debug)]
enum FileHandle {
    Write(File),
    Read(Mmap),
}

#[derive(Debug)]
pub struct GenomicDataStore<T>
where
    T: Record,
{
    index: BinningIndex,
    data_files: HashMap<String, FileHandle>,
    directory: PathBuf,
    key: Option<String>,
    results_buffer: Vec<T>,
    _phantom: PhantomData<T>,
}

impl<T: Record> GenomicDataStore<T> {
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
            index: BinningIndex::new(schema),
            data_files: HashMap::new(),
            directory: directory.to_path_buf(),
            key,
            results_buffer: Vec::with_capacity(1000),
            _phantom: PhantomData,
        })
    }

    fn get_or_create_file(&mut self, chrom: &str) -> std::io::Result<&mut File> {
        if !self.data_files.contains_key(chrom) {
            let data_path = self.get_data_path(chrom);
            let file = File::create(&data_path)?;
            let mut writer = BufWriter::new(file);
            writer.write_all(&Self::MAGIC)?;
            writer.flush()?;
            let file = writer.into_inner()?;
            self.data_files
                .insert(chrom.to_string(), FileHandle::Write(file));
        }

        match self.data_files.get_mut(chrom).unwrap() {
            FileHandle::Write(file) => Ok(file),
            FileHandle::Read(_) => Err(io::Error::new(
                io::ErrorKind::Other,
                "File is open for reading",
            )),
        }
    }

    pub fn add_record(&mut self, chrom: &str, record: &T) -> Result<(), HgIndexError> {
        if !self.data_files.contains_key(chrom) {
            self.data_files.retain(|k, _| k == chrom);
        }

        let file = self.get_or_create_file(chrom)?;

        let length;
        let offset = {
            let mut writer = BufWriter::new(file);
            let offset = writer.stream_position()?;

            // Use Record trait instead of bincode
            let record_data = record.to_bytes();
            length = record_data.len() as u64;

            writer.write_all(&length.to_le_bytes())?;
            writer.write_all(&record_data)?;
            writer.flush()?;

            offset
        };

        self.index
            .add_feature(chrom, record.start(), record.end(), offset, length)?;
        Ok(())
    }

    // Add a method to explicitly close files
    fn close_files(&mut self) -> io::Result<()> {
        self.data_files.clear();
        Ok(())
    }

    pub fn finalize(&mut self) -> std::result::Result<(), Box<dyn std::error::Error>> {
        self.close_files()?;

        // Write index to file
        let index_path = if let Some(ref key) = self.key {
            self.directory.join(key).join(Self::INDEX_FILENAME)
        } else {
            self.directory.join(Self::INDEX_FILENAME)
        };

        self.index.finalize(index_path.as_path())?;
        Ok(())
    }

    // Get metadata if it exists
    pub fn metadata<M: for<'de> Deserialize<'de>>(&self) -> Option<M> {
        self.index.metadata()
    }

    pub fn finalize_with_metadata<M>(
        &mut self,
        metadata: &M,
    ) -> std::result::Result<(), Box<dyn std::error::Error>>
    where
        M: Serialize + for<'de> Deserialize<'de>,
    {
        self.close_files()?;

        // Write index to file
        let index_path = if let Some(ref key) = self.key {
            self.directory.join(key).join(Self::INDEX_FILENAME)
        } else {
            self.directory.join(Self::INDEX_FILENAME)
        };

        self.index
            .finalize_with_metadata(index_path.as_path(), &metadata)?;
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
            results_buffer: Vec::with_capacity(1000),
            _phantom: PhantomData,
        })
    }

    // NOTE: currently this is not faster than the version below, but
    // it maybe in some cases — needs future benchmarking.
    // pub fn open_chrom_file(&mut self, chrom: &str) -> std::io::Result<()> {
    //     if !self.data_files.contains_key(chrom) {
    //         let data_path = self.get_data_path(chrom);
    //         let mmap = unsafe {
    //             // Add MAP_POPULATE to preload pages
    //             let file = File::open(&data_path)?;
    //             let mut options = memmap2::MmapOptions::new();
    //             let mmap_opts = options.populate();
    //             mmap_opts.map(&file)?
    //         };
    //
    //         if mmap[0..4] != Self::MAGIC {
    //             return Err(std::io::Error::new(
    //                 std::io::ErrorKind::InvalidData,
    //                 "Invalid file format",
    //             ));
    //         }
    //         self.data_files
    //             .insert(chrom.to_string(), FileHandle::Read(mmap));
    //     }
    //     Ok(())
    // }

    pub fn open_chrom_file(&mut self, chrom: &str) -> std::io::Result<()> {
        if !self.data_files.contains_key(chrom) {
            let data_path = self.get_data_path(chrom);
            let file = File::open(&data_path)?;
            let mmap = unsafe { Mmap::map(&file)? };

            if mmap[0..4] != Self::MAGIC {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "Invalid file format",
                ));
            }
            self.data_files
                .insert(chrom.to_string(), FileHandle::Read(mmap));
        }
        Ok(())
    }

    // Rename to just map_overlapping since there's no batching
    pub fn map_overlapping<F>(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
        mut fun: F,
    ) -> Result<usize, HgIndexError>
    where
        F: FnMut(T::Slice<'_>) -> Result<(), HgIndexError>,
    {
        if end <= start {
            return Err(HgIndexError::InvalidInterval { start, end });
        }

        if !self.index.sequences.contains_key(chrom) {
            return Ok(0);
        }

        if self.open_chrom_file(chrom).is_err() {
            return Ok(0);
        }

        let mmap = match self.data_files.get(chrom).unwrap() {
            FileHandle::Read(mmap) => mmap,
            FileHandle::Write(_) => {
                return Err(HgIndexError::StringError("File is open for writing".into()));
            }
        };

        let offsets = self.index.find_overlapping(chrom, start, end);
        if offsets.is_empty() {
            return Ok(0);
        }

        let mut count = 0;
        for (offset, length) in offsets {
            let offset = offset as usize;
            let length = length as usize;

            if offset + 8 > mmap.len() {
                continue;
            }

            if offset + 8 + length > mmap.len() {
                continue;
            }

            // Use RecordSlice for zero-copy parsing
            let record = T::Slice::from_bytes(&mmap[offset + 8..offset + 8 + length]);
            fun(record)?;
            count += 1;
        }

        Ok(count)
    }

    pub fn get_overlapping(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<&[T], HgIndexError> {
        self.results_buffer.clear();

        if end <= start {
            return Err(HgIndexError::InvalidInterval { start, end });
        }

        if !self.index.sequences.contains_key(chrom) {
            return Ok(&self.results_buffer);
        }

        if self.open_chrom_file(chrom).is_err() {
            return Ok(&self.results_buffer);
        }

        let mmap = match self.data_files.get(chrom).unwrap() {
            FileHandle::Read(mmap) => mmap,
            FileHandle::Write(_) => {
                return Err(HgIndexError::StringError("File is open for writing".into()));
            }
        };

        let offsets = self.index.find_overlapping(chrom, start, end);
        if offsets.is_empty() {
            return Ok(&self.results_buffer);
        }

        for (offset, length) in offsets {
            let offset = offset as usize;
            let length = length as usize;

            if offset + 8 > mmap.len() {
                continue;
            }

            if offset + 8 + length > mmap.len() {
                continue;
            }

            // Parse as slice then convert to owned
            let slice = T::Slice::from_bytes(&mmap[offset + 8..offset + 8 + length]);
            self.results_buffer.push(slice.into())
        }

        Ok(&self.results_buffer)
    }

    pub fn get_overlapping_batch<'a>(
        &'a mut self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<Vec<T::Slice<'a>>, HgIndexError> {
        let mut results = Vec::new();
        if end <= start {
            return Err(HgIndexError::InvalidInterval { start, end });
        }
        if !self.index.sequences.contains_key(chrom) {
            return Ok(results);
        }
        if self.open_chrom_file(chrom).is_err() {
            return Ok(results);
        }

        let mmap = match self.data_files.get(chrom).unwrap() {
            FileHandle::Read(mmap) => mmap,
            FileHandle::Write(_) => {
                return Err(HgIndexError::StringError("File is open for writing".into()))
            }
        };

        // Get all overlapping records at once
        let offsets = self.index.find_overlapping(chrom, start, end);

        // Pre-allocate to avoid resizing
        results.reserve(offsets.len());

        // Needs more extensive benchmarking:
        let chunk = false;
        if chunk {
            // Process in chunks to improve cache utilization
            const CHUNK_SIZE: usize = 32;
            for chunk in offsets.chunks(CHUNK_SIZE) {
                for &(offset, length) in chunk {
                    let offset = offset as usize;
                    let length = length as usize;
                    let record = T::Slice::from_bytes(&mmap[offset + 8..offset + 8 + length]);
                    results.push(record);
                }
            }
        } else {
            for (offset, length) in offsets {
                let offset = offset as usize;
                let length = length as usize;
                let record = T::Slice::from_bytes(&mmap[offset + 8..offset + 8 + length]);
                results.push(record);
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

    // --- Test types ---

    // A simple test record type
    #[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
    struct TestRecord {
        start: u32,
        end: u32,
        name: String,
        score: f64,
        tags: Vec<String>,
    }

    // The slice variant of TestRecord
    #[derive(Clone, Debug, Deserialize)]
    struct TestRecordSlice<'a> {
        start: u32,
        end: u32,
        name: &'a str,
        score: f64,
        tags: Vec<&'a str>,
    }

    // Then implement Record and RecordSlice
    impl Record for TestRecord {
        type Slice<'a> = TestRecordSlice<'a>;
        fn start(&self) -> u32 {
            self.start
        }
        fn end(&self) -> u32 {
            self.end
        }
        fn to_bytes(&self) -> Vec<u8> {
            // we can use bincode here for simplicity,
            // rather than manual serialization
            bincode::serialize(self).unwrap()
        }
    }

    impl<'a> RecordSlice<'a> for TestRecordSlice<'a> {
        type Owned = TestRecord;

        fn start(&self) -> u32 {
            self.start
        }

        fn end(&self) -> u32 {
            self.end
        }

        fn from_bytes(bytes: &'a [u8]) -> Self {
            bincode::deserialize(bytes)
                .map_err(|e| HgIndexError::StringError(e.to_string()))
                .unwrap()
        }

        fn to_owned(self) -> Self::Owned {
            Self::Owned {
                start: self.start,
                end: self.end,
                name: self.name.to_owned(),
                score: self.score,
                tags: self.tags.into_iter().map(|v| v.to_string()).collect(),
            }
        }
    }

    impl From<TestRecordSlice<'_>> for TestRecord {
        fn from(slice: TestRecordSlice<'_>) -> Self {
            Self {
                start: slice.start,
                end: slice.end,
                name: slice.name.to_owned(),
                score: slice.score,
                tags: slice.tags.iter().map(|&s| s.to_owned()).collect(),
            }
        }
    }

    #[derive(Debug, Serialize, Deserialize, PartialEq)]
    struct MinimalTestRecord {
        start: u32,
        end: u32,
        score: f64,
    }

    #[derive(Debug, Deserialize)]
    struct MinimalTestRecordSlice<'a> {
        start: u32,
        end: u32,
        score: f64,
        _lifetime: PhantomData<&'a ()>,
    }

    impl Record for MinimalTestRecord {
        type Slice<'a> = MinimalTestRecordSlice<'a>;
        fn start(&self) -> u32 {
            self.start
        }
        fn end(&self) -> u32 {
            self.end
        }
        fn to_bytes(&self) -> Vec<u8> {
            bincode::serialize(self).unwrap()
        }
    }

    impl<'a> RecordSlice<'a> for MinimalTestRecordSlice<'a> {
        type Owned = MinimalTestRecord;
        fn start(&self) -> u32 {
            self.start
        }
        fn end(&self) -> u32 {
            self.end
        }
        fn from_bytes(bytes: &'a [u8]) -> Self {
            bincode::deserialize(bytes)
                .map_err(|e| HgIndexError::StringError(e.to_string()))
                .unwrap()
        }

        fn to_owned(self) -> Self::Owned {
            Self::Owned {
                start: self.start,
                end: self.end,
                score: self.score,
            }
        }
    }

    impl From<MinimalTestRecordSlice<'_>> for MinimalTestRecord {
        fn from(slice: MinimalTestRecordSlice<'_>) -> Self {
            Self {
                start: slice.start,
                end: slice.end,
                score: slice.score,
            }
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

    // --- Test Functions ---

    #[test]
    fn test_store_and_retrieve() {
        let test_dir = TestDir::new("store_and_retrieve").expect("Failed to create test dir");
        let base_dir = test_dir.path(); // Don't add test.gidx

        // An example key
        let key = "example-key".to_string();

        // Create store and add records
        let mut store = GenomicDataStore::<TestRecord>::create(base_dir, Some(key.clone()))
            .expect("Failed to create store");
        for (chrom, record) in make_test_records() {
            store
                .add_record(&chrom, &record)
                .expect("Failed to add record");
        }

        store.finalize().expect("Failed to finalize store");

        let mut store = GenomicDataStore::<TestRecord>::open(&base_dir, Some(key.clone()))
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
        let result = GenomicDataStore::<TestRecord>::open(&bad_file, None);
        assert!(result.is_err());
    }

    #[test]
    fn test_empty_regions() {
        let test_dir = TestDir::new("empty_regions").expect("Failed to create test dir");
        let store_path = test_dir.path().join("empty.gidx");

        let mut store = GenomicDataStore::<TestRecord>::create(&store_path, None)
            .expect("Failed to create store");

        store.finalize().expect("Failed to finalize store");

        // Query empty store
        let mut store =
            GenomicDataStore::<TestRecord>::open(&store_path, None).expect("Failed to open store");

        let results = store.get_overlapping("chr1", 0, 1000).unwrap();
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_concurrent_reads() {
        use std::sync::Arc;
        use std::thread;

        let test_dir = TestDir::new("concurrent").expect("Failed to create test dir");
        let store_path = test_dir.path().join("test.gidx");

        // Create and populate store
        {
            let mut store = GenomicDataStore::<MinimalTestRecord>::create(&store_path, None)
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
                    let mut store = GenomicDataStore::<MinimalTestRecord>::open(&path, None)
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
    fn test_map_vs_get_consistency() {
        let test_dir = TestDir::new("map_vs_get_consistency").expect("Failed to create test dir");
        let base_dir = test_dir.path();

        // Create the store and add test records
        let mut store = GenomicDataStore::<TestRecord>::create(base_dir, None)
            .expect("Failed to create GenomicDataStore");
        for (chrom, record) in make_test_records() {
            store
                .add_record(&chrom, &record)
                .expect("Failed to add record");
        }
        store.finalize().expect("Failed to finalize store");

        // Reopen the finalized store
        let mut store = GenomicDataStore::<TestRecord>::open(base_dir, None)
            .expect("Failed to open GenomicDataStore");

        // Define test queries
        let queries = vec![
            ("chr1", 1200, 1800),
            ("chr1", 0, 3000),
            ("chr2", 50000, 60000),
            ("chr2", 55000, 58000),
            ("chr3", 0, 10000),
        ];

        for (chrom, start, end) in queries {
            // Get overlapping records
            let get_results = store.get_overlapping(chrom, start, end).unwrap().to_vec();

            // Map overlapping records
            let mut map_results = Vec::new();
            store
                .map_overlapping(chrom, start, end, |record| {
                    map_results.push(record.to_owned());
                    Ok(())
                })
                .unwrap();

            // Assert that both results are identical
            assert_eq!(
                get_results, map_results,
                "Mismatch for chrom: {}, start: {}, end: {}",
                chrom, start, end
            );
        }
    }

    #[test]
    fn test_metadata_storage_and_retrieval() {
        use std::collections::HashMap;
        let test_dir = TestDir::new("metadata_test").expect("Failed to create test dir");
        let base_dir = test_dir.path();

        // Create some test metadata (using a simple struct)
        #[derive(Debug, Serialize, Deserialize, PartialEq)]
        struct TestMetadata {
            name: String,
            values: HashMap<String, i32>,
        }

        let original_metadata = TestMetadata {
            name: "test".to_string(),
            values: {
                let mut m = HashMap::new();
                m.insert("key1".to_string(), 42);
                m.insert("key2".to_string(), 100);
                m
            },
        };

        // Create and populate store
        {
            let mut store = GenomicDataStore::<TestRecord>::create(base_dir, None)
                .expect("Failed to create store");

            // Add some test records
            let record = TestRecord {
                start: 1000,
                end: 2000,
                name: "feature1".to_string(),
                score: 0.5,
                tags: vec!["test".to_string()],
            };
            store
                .add_record("chr1", &record)
                .expect("Failed to add record");

            // Finalize with metadata
            store
                .finalize_with_metadata(&original_metadata)
                .expect("Failed to finalize with metadata");
        }

        // Reopen and check metadata
        {
            let store =
                GenomicDataStore::<TestRecord>::open(base_dir, None).expect("Failed to open store");

            let retrieved_metadata: Option<TestMetadata> = store.metadata();
            assert!(retrieved_metadata.is_some());

            let retrieved_metadata = retrieved_metadata.unwrap();
            assert_eq!(retrieved_metadata, original_metadata);
            assert_eq!(retrieved_metadata.name, "test");
            assert_eq!(retrieved_metadata.values.get("key1"), Some(&42));
            assert_eq!(retrieved_metadata.values.get("key2"), Some(&100));
        }
    }
}
