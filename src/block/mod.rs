use std::{
    fs::File,
    io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write},
    marker::PhantomData,
};

use serde::{Deserialize, Serialize};

use crate::error::HgIndexError;

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

#[derive(Debug)]
pub struct BlockWriter {
    inner: BufWriter<File>,
    buffer: Vec<u8>,
    current_block_size: usize,
    current_file_offset: u64,
}

impl BlockWriter {
    pub fn new(file: File) -> io::Result<Self> {
        Ok(Self {
            inner: BufWriter::new(file),
            buffer: Vec::with_capacity(BLOCK_SIZE),
            current_block_size: 0,
            current_file_offset: 0,
        })
    }

    /// Add a new record.
    pub fn add_record<T>(&mut self, record: &T) -> Result<VirtualOffset, HgIndexError>
    where
        T: Serialize,
    {
        let record_data = bincode::serialize(record)
            .map_err(|e| HgIndexError::SerializationError(e.to_string()))?;
        let length = record_data.len() as u64;

        // Save starting position before we write anything
        let start_vo = VirtualOffset::new(self.current_file_offset, self.current_block_size as u16);

        // First write the length prefix
        let length_bytes = length.to_le_bytes();
        let mut bytes_written = 0;

        // Write length prefix, might itself span blocks
        while bytes_written < length_bytes.len() {
            if self.current_block_size >= BLOCK_SIZE {
                self.flush_block()?;
            }
            let can_write =
                (BLOCK_SIZE - self.current_block_size).min(length_bytes.len() - bytes_written);
            self.buffer
                .extend_from_slice(&length_bytes[bytes_written..bytes_written + can_write]);
            self.current_block_size += can_write;
            bytes_written += can_write;
        }

        // Then write actual record data
        bytes_written = 0;
        while bytes_written < record_data.len() {
            if self.current_block_size >= BLOCK_SIZE {
                self.flush_block()?;
            }
            let can_write =
                (BLOCK_SIZE - self.current_block_size).min(record_data.len() - bytes_written);
            self.buffer
                .extend_from_slice(&record_data[bytes_written..bytes_written + can_write]);
            self.current_block_size += can_write;
            bytes_written += can_write;
        }

        Ok(start_vo)
    }

    fn flush_block(&mut self) -> Result<(), HgIndexError> {
        if self.buffer.is_empty() {
            return Ok(());
        }

        let compressed = zstd::encode_all(&self.buffer[..], COMPRESSION_LEVEL)
            .map_err(|e| HgIndexError::IOError(io::Error::new(io::ErrorKind::Other, e)))?;

        // Write block header
        let uncomp_size = self.buffer.len() as u32;
        let comp_size = compressed.len() as u32;

        self.inner.write_all(&uncomp_size.to_le_bytes())?;
        self.inner.write_all(&comp_size.to_le_bytes())?;
        self.inner.write_all(&compressed)?;

        self.current_file_offset += (8 + compressed.len()) as u64;
        self.buffer.clear();
        self.current_block_size = 0;

        Ok(())
    }

    pub fn finish(mut self) -> Result<(), HgIndexError> {
        self.flush_block()?;
        Ok(self.inner.flush()?)
    }
}

#[derive(Debug)]
pub struct BlockReader<T> {
    inner: BufReader<File>,
    current_block: Vec<u8>,
    current_block_offset: u64,
    last_virtual_offset: Option<VirtualOffset>,
    _phantom: PhantomData<T>,
}

impl<T> BlockReader<T>
where
    T: for<'de> Deserialize<'de>,
{
    /// Create a new [`BlockReader`].
    pub fn new(file: File) -> io::Result<Self> {
        Ok(Self {
            inner: BufReader::new(file),
            current_block: Vec::with_capacity(BLOCK_SIZE),
            current_block_offset: 0,
            last_virtual_offset: None,
            _phantom: PhantomData,
        })
    }

    // Seek a block at the provided virtual offset.
    fn seek_virtual(&mut self, voffset: VirtualOffset) -> io::Result<()> {
        // Skip if we're already at the right block for performance
        if let Some(last_offset) = self.last_virtual_offset {
            if last_offset.file_offset() == voffset.file_offset() {
                return Ok(());
            }
        }

        let file_offset = voffset.file_offset();
        self.inner.seek(SeekFrom::Start(file_offset))?;

        // Read 8-byte block header:
        // - 4 bytes: uncompressed size (u32 little-endian)
        // - 4 bytes: compressed size (u32 little-endian)
        let mut header = [0u8; 8];
        self.inner.read_exact(&mut header)?;

        let uncomp_size = u32::from_le_bytes(header[0..4].try_into().unwrap());
        let comp_size = u32::from_le_bytes(header[4..8].try_into().unwrap());

        // Read the actual compressed data block
        let mut compressed = vec![0u8; comp_size as usize];
        self.inner.read_exact(&mut compressed)?;

        // Decompress the block
        self.current_block = zstd::decode_all(&compressed[..])
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

        if self.current_block.len() != uncomp_size as usize {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Block size mismatch after decompression",
            ));
        }

        self.current_block_offset = file_offset;
        self.last_virtual_offset = Some(voffset);

        Ok(())
    }

    // Read the provided length of bytes, starting at the provided virtual offset.
    fn read_bytes(&mut self, voffset: VirtualOffset, length: u64) -> io::Result<Vec<u8>> {
        let mut result = Vec::new();
        let mut bytes_read = 0;
        let mut current_voffset = voffset;

        while bytes_read < length {
            self.seek_virtual(current_voffset)?;

            let block_offset = current_voffset.block_offset() as usize;
            if block_offset >= self.current_block.len() {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "block offset beyond block size",
                ));
            }

            let available = self.current_block.len() - block_offset;
            let to_read = ((length - bytes_read) as usize).min(available);

            result.extend_from_slice(&self.current_block[block_offset..block_offset + to_read]);
            bytes_read += to_read as u64;

            if bytes_read < length {
                // Get current block's compressed size from the header
                self.inner
                    .seek(SeekFrom::Start(current_voffset.file_offset() + 4))?; // Skip uncomp_size
                let mut comp_size_bytes = [0u8; 4];
                self.inner.read_exact(&mut comp_size_bytes)?;
                let comp_size = u32::from_le_bytes(comp_size_bytes);

                // Calculate next block's file offset
                current_voffset =
                    VirtualOffset::new(current_voffset.file_offset() + 8 + comp_size as u64, 0);
            }
        }

        Ok(result)
    }

    /// Read a single record at the provided virtual offset (even if it spans multiple
    /// compressed blocks).
    pub fn read_record(&mut self, voffset: VirtualOffset) -> Result<T, HgIndexError>
    where
        T: for<'de> serde::Deserialize<'de>,
    {
        // First read length prefix
        let mut length_bytes = [0u8; LENGTH_PREFIX_SIZE];
        let prefix_bytes = self.read_bytes(voffset, LENGTH_PREFIX_SIZE as u64)?;
        length_bytes.copy_from_slice(&prefix_bytes);
        let record_length = u64::from_le_bytes(length_bytes);

        // Read the actual record data
        let record_bytes = self.read_bytes(
            VirtualOffset::new(
                voffset.file_offset(),
                voffset.block_offset() + LENGTH_PREFIX_SIZE as u16,
            ),
            record_length,
        )?;

        // Deserialize
        bincode::deserialize(&record_bytes)
            .map_err(|e| HgIndexError::SerializationError(e.to_string()))
    }
}

#[cfg(test)]
mod tests {
    use crate::test_utils::test_utils::TestDir;

    use super::*;
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
}
