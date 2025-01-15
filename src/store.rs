use dashmap::DashMap;
use std::io;
use std::io::Write as IoWrite;
use std::sync::Arc;
use std::{
    fs::{self, File},
    io::{BufWriter, Seek},
    marker::PhantomData,
    path::{Path, PathBuf},
};

use rayon::prelude::*;

use memmap2::Mmap;

use crate::{error::HgIndexError, index::BinningIndex, BinningSchema, SerdeType};
use crate::{GenomicCoordinates, Metadata};

#[derive(Debug)]
enum FileHandle {
    Write(File),
    Read(Mmap),
}

#[derive(Debug)]
pub struct GenomicDataStore<T>
where
    T: GenomicCoordinates + SerdeType + Send + Sync,
{
    index: Arc<BinningIndex>,
    data_files: Arc<DashMap<String, FileHandle>>,
    directory: PathBuf,
    key: Option<String>,
    _phantom: PhantomData<T>,
}

impl<T> GenomicDataStore<T>
where
    T: GenomicCoordinates + SerdeType + Send + Sync,
{
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
            index: Arc::new(BinningIndex::from_schema(schema)),
            data_files: Arc::new(DashMap::new()),
            directory: directory.to_path_buf(),
            key,
            _phantom: PhantomData,
        })
    }

    fn get_or_create_file(
        &self,
        chrom: &str,
    ) -> io::Result<dashmap::mapref::one::RefMut<'_, String, FileHandle>> {
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

        let file_ref = self
            .data_files
            .get_mut(chrom)
            .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Failed to get file handle"))?;

        match &*file_ref {
            FileHandle::Write(_) => Ok(file_ref),
            FileHandle::Read(_) => Err(io::Error::new(
                io::ErrorKind::Other,
                "File is open for reading",
            )),
        }
    }

    pub fn add_record(&self, chrom: &str, record: &T) -> Result<(), HgIndexError> {
        let file_ref = self.get_or_create_file(chrom)?;

        let offset = match &*file_ref {
            FileHandle::Write(file) => {
                let mut writer = BufWriter::new(file);
                let offset = writer.stream_position()?;
                let record_data = bincode::serialize(record).unwrap();
                let length = record_data.len() as u64;
                writer.write_all(&length.to_le_bytes())?;
                writer.write_all(&record_data)?;
                writer.flush()?;
                offset
            }
            _ => {
                return Err(HgIndexError::StringError(
                    "File not opened for writing".into(),
                ))
            }
        };

        // Add feature directly to the chromosome-level index
        self.index
            .add_feature(chrom, record.start(), record.end(), offset)?;
        Ok(())
    }

    /// Get the current index path
    pub fn index_path(&self) -> PathBuf {
        if let Some(ref key) = self.key {
            self.directory.join(key).join(Self::INDEX_FILENAME)
        } else {
            self.directory.join(Self::INDEX_FILENAME)
        }
    }

    pub fn finalize(&self) -> Result<(), HgIndexError> {
        // Clear all file handles before finalizing
        self.data_files.clear();

        // Get the current index and serialize it with the provided metadata
        let index_path = self.index_path();
        self.index.serialize_to_path(&index_path, None::<&()>)?;
        Ok(())
    }

    pub fn finalize_with_metadata<M>(&self, metadata: &M) -> Result<(), Box<dyn std::error::Error>>
    where
        M: Metadata,
    {
        // Clear all file handles before finalizing
        self.data_files.clear();

        // Get the current index and serialize it with the provided metadata
        let index_path = self.index_path();
        self.index.serialize_to_path(&index_path, Some(metadata))?;

        Ok(())
    }

    pub fn open(directory: &Path, key: Option<String>) -> Result<Self, Box<dyn std::error::Error>> {
        let target_dir = if let Some(ref key) = key {
            directory.join(key)
        } else {
            directory.to_path_buf()
        };

        let index_path = target_dir.join(Self::INDEX_FILENAME);
        let index = BinningIndex::open(&index_path)?;

        Ok(Self {
            index: Arc::new(index),
            data_files: Arc::new(DashMap::new()),
            directory: directory.to_path_buf(),
            key,
            _phantom: PhantomData,
        })
    }

    fn open_chrom_file(&self, chrom: &str) -> io::Result<()> {
        if !self.data_files.contains_key(chrom) {
            let data_path = self.get_data_path(chrom);
            let file = File::open(&data_path)?;
            let mmap = unsafe { Mmap::map(&file)? };

            if mmap[0..4] != Self::MAGIC {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Invalid file format",
                ));
            }

            self.data_files
                .insert(chrom.to_string(), FileHandle::Read(mmap));
        }
        Ok(())
    }

    fn get_mmap_for_chrom(
        &self,
        chrom: &str,
    ) -> Result<dashmap::mapref::one::Ref<'_, String, FileHandle>, HgIndexError> {
        self.open_chrom_file(chrom)?;

        let guard = self
            .data_files
            .get(chrom)
            .ok_or_else(|| HgIndexError::StringError("Chromosome not found".into()))?;

        match &*guard {
            FileHandle::Read(_) => Ok(guard),
            _ => Err(HgIndexError::StringError("File is open for writing".into())),
        }
    }

    /// Alias
    pub fn get_overlapping(
        &self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<Vec<T>, HgIndexError> {
        self.get_overlapping_batch(chrom, start, end)
    }

    pub fn get_overlapping_batch(
        &self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<Vec<T>, HgIndexError> {
        // Validate the input interval
        if end <= start {
            return Err(HgIndexError::InvalidInterval { start, end });
        }

        // Clone the index for thread safety
        let index = self.index.clone();

        // Check if the chromosome exists in the index
        if !index.sequences.contains_key(chrom) {
            return Ok(Vec::new());
        }

        // Get the memory map for the specified chromosome
        let mmap_ref = self.get_mmap_for_chrom(chrom)?;
        let mmap = match &*mmap_ref {
            FileHandle::Read(mmap) => mmap,
            _ => return Err(HgIndexError::StringError("File is open for writing".into())),
        };

        // Find the offsets for the overlapping regions
        let offsets: Vec<u64> = index.find_overlapping(chrom, start, end).collect();
        if offsets.is_empty() {
            return Ok(Vec::new());
        }

        // Get the first and last offsets
        let first_offset = offsets[0] as usize;
        let last_offset = *offsets.last().unwrap() as usize;

        // Determine the length of the last record
        if last_offset + 8 > mmap.len() {
            return Err(HgIndexError::StringError("Invalid offset bounds".into()));
        }
        let last_length =
            u64::from_le_bytes(mmap[last_offset..last_offset + 8].try_into().map_err(|_| {
                HgIndexError::StringError("Failed to parse length of the last record".into())
            })?) as usize;

        // Define the end boundary of the region to process
        let region_end = last_offset + 8 + last_length;

        // Collect the results
        let mut results = Vec::with_capacity(offsets.len());
        let mut current_offset = first_offset;

        while current_offset < region_end {
            // Ensure we stay within the bounds of the memory map
            if current_offset + 8 > mmap.len() {
                break;
            }

            // Get the length of the current record
            let length = u64::from_le_bytes(
                mmap[current_offset..current_offset + 8]
                    .try_into()
                    .map_err(|_| {
                        HgIndexError::StringError("Failed to parse record length".into())
                    })?,
            ) as usize;

            // Ensure the full record fits within the memory map
            if current_offset + 8 + length > mmap.len() {
                break;
            }

            // Check if the current offset is in the list of offsets
            if offsets.contains(&(current_offset as u64)) {
                // Deserialize the record and add it to the results
                if let Ok(record) = bincode::deserialize::<T>(
                    &mmap[current_offset + 8..current_offset + 8 + length],
                ) {
                    results.push(record);
                }
            }

            // Move to the next record
            current_offset += 8 + length;
        }

        Ok(results)
    }

    pub fn map_overlapping_batch<F>(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
        mut fun: F,
    ) -> Result<usize, HgIndexError>
    where
        F: FnMut(&[u8]) -> Result<(), HgIndexError>,
    {
        // Validate interval
        if end <= start {
            return Err(HgIndexError::InvalidInterval { start, end });
        }

        // Early exit if chromosome not found
        if !self.index.sequences.contains_key(chrom) {
            return Ok(0);
        }

        // Open and validate the mmap
        self.open_chrom_file(chrom)?;
        let file_ref = self.data_files.get(chrom).unwrap();
        let mmap = match &*file_ref {
            // Deref the DashMap::Ref here
            FileHandle::Read(mmap) => mmap,
            FileHandle::Write(_) => {
                return Err(HgIndexError::StringError("File is open for writing".into()))
            }
        };

        // Get overlapping offsets and collect them since we can't check is_empty on iterator
        let offsets: Vec<u64> = self.index.find_overlapping(chrom, start, end).collect();
        if offsets.is_empty() {
            return Ok(0);
        }

        let mut count = 0;
        for offset in offsets {
            let offset = offset as usize;
            if offset + 8 > mmap.len() {
                continue;
            }

            // Get record length
            let length =
                u64::from_le_bytes(mmap[offset..offset + 8].try_into().map_err(|_| {
                    HgIndexError::StringError("Failed to parse record length".into())
                })?) as usize;

            if offset + 8 + length > mmap.len() {
                continue;
            }

            // Process record bytes with provided function
            fun(&mmap[offset + 8..offset + 8 + length])?;
            count += 1;
        }

        Ok(count)
    }

    pub fn map_overlapping_parallel<F>(
        &self,
        regions: Vec<(String, u32, u32)>,
        fun: F,
    ) -> Result<Vec<usize>, HgIndexError>
    where
        F: Fn(&[u8]) -> Result<(), HgIndexError> + Sync,
    {
        // Parallel processing of regions
        let results: Result<Vec<usize>, HgIndexError> = regions
            .into_par_iter()
            .map(|(chrom, start, end)| {
                // For each region, process overlaps with `fun`
                self.map_overlapping_batch(&chrom, start, end, fun.clone())
            })
            .collect();

        results
    }

    // New method for parallel querying of multiple regions
    pub fn get_overlapping_parallel(
        &self,
        regions: Vec<(String, u32, u32)>,
    ) -> Result<Vec<Vec<T>>, HgIndexError> {
        use rayon::prelude::*;

        regions
            .into_par_iter()
            .map(|(chrom, start, end)| self.get_overlapping_batch(&chrom, start, end))
            .collect()
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
        let store = GenomicDataStore::<TestRecord>::create(base_dir, Some(key.clone()))
            .expect("Failed to create store");
        for (chrom, record) in make_test_records() {
            store
                .add_record(&chrom, &record)
                .expect("Failed to add record");
        }

        store.finalize().expect("Failed to finalize store");

        let store = GenomicDataStore::<TestRecord>::open(&base_dir, Some(key.clone()))
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

        let store = GenomicDataStore::<TestRecord>::create(&store_path, None)
            .expect("Failed to create store");

        store.finalize().expect("Failed to finalize store");

        // Query empty store
        let store =
            GenomicDataStore::<TestRecord>::open(&store_path, None).expect("Failed to open store");

        let results = store.get_overlapping("chr1", 0, 1000).unwrap();
        assert_eq!(results.len(), 0);
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
            let store = GenomicDataStore::<MinimalTestRecord>::create(&store_path, None)
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
                    let store = GenomicDataStore::<MinimalTestRecord>::open(&path, None)
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
