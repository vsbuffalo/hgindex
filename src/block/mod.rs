// block/mod.rs

mod reader;
mod writer;

pub use reader::BlockReader;
pub use writer::BlockWriter;

pub const BLOCK_SIZE: usize = 65536; // 2^16 bytes like BGZF
pub const COMPRESSION_LEVEL: i32 = 3;
const LENGTH_PREFIX_SIZE: usize = 8; // size of u64 length prefix in bytes

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct VirtualOffset {
    value: u64,
}

impl VirtualOffset {
    pub fn new(file_offset: u64, block_offset: u16) -> Self {
        Self {
            value: (file_offset << 16) | (block_offset as u64),
        }
    }

    pub fn file_offset(&self) -> u64 {
        self.value >> 16
    }

    pub fn block_offset(&self) -> u16 {
        (self.value & 0xFFFF) as u16
    }
}

impl From<VirtualOffset> for u64 {
    fn from(value: VirtualOffset) -> Self {
        value.value
    }
}

impl From<u64> for VirtualOffset {
    fn from(value: u64) -> Self {
        VirtualOffset { value }
    }
}

#[cfg(test)]
mod tests {
    use std::fs::File;

    use crate::{error::HgIndexError, test_utils::test_utils::TestDir};

    use super::*;
    use rand::Rng;
    use serde::{Deserialize, Serialize};

    #[test]
    fn test_virtual_offset() {
        let voffset = VirtualOffset::new(1234, 56);
        assert_eq!(voffset.file_offset(), 1234);
        assert_eq!(voffset.block_offset(), 56);
    }

    #[derive(Debug, Serialize, Deserialize, PartialEq)]
    struct TestRecord {
        id: u32,
        data: String,
    }

    fn create_test_record(id: u32, data_size: usize) -> TestRecord {
        TestRecord {
            id,
            data: "x".repeat(data_size),
        }
    }

    #[test]
    fn test_single_small_record() -> Result<(), HgIndexError> {
        let test_dir = TestDir::new("single_small_record").unwrap();
        let path = test_dir.path().join("test.idx");

        // Write a single small record
        let record = create_test_record(1, 100);
        let file = File::create(&path)?;
        let mut writer = BlockWriter::new(file)?;
        let offset = writer.add_record(&record)?;
        writer.finish()?;

        // Read it back
        let file = File::open(&path)?;
        let mut reader = BlockReader::new(file)?;
        let read_record: TestRecord = reader.read_record(offset)?;

        assert_eq!(record, read_record);
        Ok(())
    }

    #[test]
    fn test_multiple_records_single_block() -> Result<(), HgIndexError> {
        let test_dir = TestDir::new("multiple_records").unwrap();
        let path = test_dir.path().join("test.idx");

        // Write multiple small records that should fit in one block
        let records: Vec<_> = (0..5).map(|i| create_test_record(i, 1000)).collect();
        let file = File::create(&path)?;
        let mut writer = BlockWriter::new(file)?;

        let offsets: Vec<_> = records
            .iter()
            .map(|rec| writer.add_record(rec))
            .collect::<Result<_, _>>()?;
        writer.finish()?;

        // Read them back
        let file = File::open(&path)?;
        let mut reader = BlockReader::new(file)?;

        for (record, offset) in records.iter().zip(offsets.iter()) {
            let read_record: TestRecord = reader.read_record(*offset)?;
            assert_eq!(record, &read_record);
        }

        Ok(())
    }

    #[test]
    fn test_block_spanning_record() -> Result<(), HgIndexError> {
        let test_dir = TestDir::new("block_spanning").unwrap();
        let path = test_dir.path().join("test.idx");

        // Create a record larger than BLOCK_SIZE to ensure it spans blocks
        let large_record = create_test_record(1, BLOCK_SIZE * 2);

        let file = File::create(&path)?;
        let mut writer = BlockWriter::new(file)?;
        let offset = writer.add_record(&large_record)?;
        writer.finish()?;

        // Read it back
        let file = File::open(&path)?;
        let mut reader = BlockReader::new(file)?;
        let read_record: TestRecord = reader.read_record(offset)?;

        assert_eq!(large_record, read_record);
        Ok(())
    }

    #[test]
    fn test_mixed_record_sizes() -> Result<(), HgIndexError> {
        let test_dir = TestDir::new("mixed_sizes").unwrap();
        let path = test_dir.path().join("test.idx");

        // Create mix of small and large records
        let records = vec![
            create_test_record(1, 100),            // small
            create_test_record(2, BLOCK_SIZE * 2), // spans blocks
            create_test_record(3, 200),            // small
            create_test_record(4, BLOCK_SIZE),     // exactly one block
        ];

        let file = File::create(&path)?;
        let mut writer = BlockWriter::new(file)?;
        let offsets: Vec<_> = records
            .iter()
            .map(|rec| writer.add_record(rec))
            .collect::<Result<_, _>>()?;
        writer.finish()?;

        // Read them back in random order to test seeking
        let file = File::open(&path)?;
        let mut reader = BlockReader::new(file)?;

        // Read in reverse order to test seeking
        for (record, offset) in records.iter().zip(offsets.iter()).rev() {
            let read_record: TestRecord = reader.read_record(*offset)?;
            assert_eq!(record, &read_record);
        }

        Ok(())
    }

    #[test]
    fn test_length_prefix_spanning() -> Result<(), HgIndexError> {
        let test_dir = TestDir::new("length_prefix_span").unwrap();
        let path = test_dir.path().join("test.idx");

        // Create a record that will cause the length prefix to span a block boundary
        let file = File::create(&path)?;
        let mut writer = BlockWriter::new(file)?;

        // First fill up most of a block
        let filler = create_test_record(1, BLOCK_SIZE - 4); // leave 4 bytes
        writer.add_record(&filler)?;

        // Now add another record - its length prefix should span blocks
        let record = create_test_record(2, 1000);
        let offset = writer.add_record(&record)?;
        writer.finish()?;

        // Read it back
        let file = File::open(&path)?;
        let mut reader = BlockReader::new(file)?;
        let read_record: TestRecord = reader.read_record(offset)?;

        assert_eq!(record, read_record);
        Ok(())
    }

    #[test]
    fn test_flush_partial_block() -> Result<(), HgIndexError> {
        let test_dir = TestDir::new("partial_block").unwrap();
        let path = test_dir.path().join("test.idx");

        let file = File::create(&path)?;
        let mut writer = BlockWriter::new(file)?;

        // Write a small record that won't fill block
        let small_record = create_test_record(1, 100);
        let offset = writer.add_record(&small_record)?; // Save the offset
                                                        // Finalize immediately - should handle partial block correctly
        writer.finish()?;

        // Read it back to verify
        let file = File::open(&path)?;
        let mut reader = BlockReader::new(file)?;
        let read_record: TestRecord = reader.read_record(offset)?;
        assert_eq!(small_record, read_record);
        Ok(())
    }

    #[test]
    fn test_many_small_records() -> Result<(), HgIndexError> {
        let test_dir = TestDir::new("many_records").unwrap();
        let path = test_dir.path().join("test.idx");

        let file = File::create(&path)?;
        let mut writer = BlockWriter::new(file)?;

        // Create many small records
        let num_records = 10_000;
        let records: Vec<_> = (0..num_records)
            .map(|i| create_test_record(i as u32, 10))
            .collect();

        // Write them all
        let offsets: Vec<_> = records
            .iter()
            .map(|rec| writer.add_record(rec))
            .collect::<Result<_, _>>()?;
        writer.finish()?;

        // Read random access to verify
        let file = File::open(&path)?;
        let mut reader = BlockReader::new(file)?;

        // Check random sampling
        let mut rng = rand::thread_rng();
        for _ in 0..100 {
            let idx = rng.gen_range(0..num_records);
            let read_record: TestRecord = reader.read_record(offsets[idx])?;
            assert_eq!(records[idx], read_record);
        }

        Ok(())
    }

    #[test]
    fn test_block_spanning_with_overflow() -> Result<(), HgIndexError> {
        let test_dir = TestDir::new("block_overflow_spanning").unwrap();
        let path = test_dir.path().join("test.idx");

        // Create a record that will be written near the end of a block
        let record = create_test_record(1, 100);

        let file = File::create(&path)?;
        let mut writer = BlockWriter::new(file)?;

        // First fill up most of a block
        let filler = create_test_record(0, BLOCK_SIZE - 10); // leave just a few bytes
        writer.add_record(&filler)?;

        // Add our test record which should span blocks
        let offset = writer.add_record(&record)?;
        writer.finish()?;

        // Read it back
        let file = File::open(&path)?;
        let mut reader = BlockReader::new(file)?;

        // This should now work without error
        let read_record: TestRecord = reader.read_record(offset)?;
        assert_eq!(record, read_record);

        Ok(())
    }
}
