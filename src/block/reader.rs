/// block/reader.rs
use super::CompressedBlock;
use crate::SerdeType;
use std::io::{self, Read, Seek, SeekFrom};

pub struct BlockReader<R: Read + Seek> {
    /// The underlying reader
    reader: R,
}

impl<R: Read + Seek> BlockReader<R> {
    /// Create a new BlockReader
    pub fn new(reader: R) -> Self {
        Self { reader }
    }

    /// Read a block header from the current position
    pub fn read_header(&mut self) -> io::Result<CompressedBlock> {
        // Read compressed size
        let mut size_buf = [0u8; 4];
        self.reader.read_exact(&mut size_buf)?;
        let compressed_size = u32::from_le_bytes(size_buf);

        // Read uncompressed size
        self.reader.read_exact(&mut size_buf)?;
        let uncompressed_size = u32::from_le_bytes(size_buf);

        // Read n_records
        self.reader.read_exact(&mut size_buf)?;
        let n_records = u32::from_le_bytes(size_buf);

        // Read start position
        self.reader.read_exact(&mut size_buf)?;
        let start = u32::from_le_bytes(size_buf);

        // Read end position
        self.reader.read_exact(&mut size_buf)?;
        let end = u32::from_le_bytes(size_buf);

        // Read compressed data
        let mut data = vec![0u8; compressed_size as usize];
        self.reader.read_exact(&mut data)?;

        Ok(CompressedBlock {
            compressed_size,
            uncompressed_size,
            n_records,
            data,
            start,
            end,
        })
    }

    /// Read and decompress records from a block.
    ///
    /// Returns a `Vec` of tuples containing the start position,
    /// end position, and the decompressed record data in its original
    /// type `T`.
    pub fn read_records<T: SerdeType>(
        &mut self,
        header: &CompressedBlock,
    ) -> io::Result<Vec<(u32, u32, T)>> {
        // Decompress the data
        let mut decompressed = Vec::new();
        zstd::stream::copy_decode(header.data.as_slice(), &mut decompressed)?;

        // Deserialize records
        let records = bincode::deserialize(&decompressed)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

        Ok(records)
    }

    /// Seed to a specific position in the file.
    pub fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.reader.seek(pos)
    }
}

#[cfg(test)]
mod tests {
    use crate::block::{write_block, BlockConfig, BlockWriter};

    use super::*;
    use serde::{Deserialize, Serialize};
    use std::io::Cursor;

    #[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
    struct TestRecord {
        value: i32,
    }

    #[test]
    fn test_block_read_write() {
        // First create some test data using BlockWriter
        let mut writer = BlockWriter::new(BlockConfig {
            max_records: 2,
            level: 3,
        })
        .unwrap();

        let test_records = vec![
            (1000u32, 2000u32, TestRecord { value: 42 }),
            (1500u32, 2500u32, TestRecord { value: 43 }),
        ];

        // Add records to get a block
        let block = {
            let mut block = None;
            for (start, end, record) in test_records.clone().into_iter() {
                if let Some(b) = writer.add_record(start, end, record).unwrap() {
                    block = Some(b);
                    break;
                }
            }
            block.unwrap()
        };

        // Now test reading it back
        let mut buffer: Vec<u8> = Vec::new();

        // Write header fields
        buffer.extend_from_slice(&block.compressed_size.to_le_bytes());
        buffer.extend_from_slice(&block.uncompressed_size.to_le_bytes());
        buffer.extend_from_slice(&block.n_records.to_le_bytes());
        buffer.extend_from_slice(&block.start.to_le_bytes());
        buffer.extend_from_slice(&block.end.to_le_bytes());
        buffer.extend_from_slice(&block.data);

        // Create a cursor to the LE byte array.
        let cursor = Cursor::new(buffer);
        let mut reader = BlockReader::new(cursor);

        // Read header
        let header = reader.read_header().unwrap();
        assert_eq!(header.compressed_size, block.compressed_size);
        assert_eq!(header.uncompressed_size, block.uncompressed_size);
        assert_eq!(header.n_records, block.n_records);
        assert_eq!(header.start, block.start);
        assert_eq!(header.end, block.end);

        // Read records
        let records: Vec<(u32, u32, TestRecord)> = reader.read_records(&header).unwrap();
        assert_eq!(records, test_records);
    }

    #[test]
    fn test_block_round_trip() {
        use std::io::Cursor;

        // Create test record
        let mut writer = BlockWriter::new(BlockConfig {
            max_records: 1, // Force flush after 1 record
            level: 3,
        })
        .unwrap();

        let test_record = TestRecord { value: 42 };

        // First attempt to get block via add_record
        let block = writer
            .add_record(5000000, 5000100, test_record.clone())
            .unwrap();

        // If None, try to finish the writer to get the block
        let block = if let Some(block) = block {
            block
        } else {
            writer
                .finish()
                .unwrap()
                .expect("Should get block from finish")
        };

        // println!("Got block with {} records", block.n_records);

        // Verify block metadata
        assert_eq!(block.start, 5000000);
        assert_eq!(block.end, 5000100);
        assert_eq!(block.n_records, 1);

        // Write block to buffer
        let mut buffer = Vec::new();
        write_block(&mut buffer, &block).unwrap();

        // Read back using BlockReader
        let cursor = Cursor::new(buffer);
        let mut reader = BlockReader::new(cursor);

        // Read and verify header
        let header = reader.read_header().unwrap();
        assert_eq!(header.start, 5000000);
        assert_eq!(header.end, 5000100);
        assert_eq!(header.n_records, 1);

        // Read and verify records
        let records: Vec<(u32, u32, TestRecord)> = reader.read_records(&header).unwrap();
        assert_eq!(records.len(), 1);
        let (start, end, record) = &records[0];
        assert_eq!(*start, 5000000);
        assert_eq!(*end, 5000100);
        assert_eq!(record, &test_record);
    }
}
