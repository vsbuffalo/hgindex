/// block/writer.rs
use super::{BlockBuffer, BlockConfig, CompressedBlock};
use crate::SerdeType;
use core::fmt;
use std::io;

/// Block-based compressed writer.
pub struct BlockWriter<T: SerdeType> {
    buffer: BlockBuffer<T>,
    config: BlockConfig,
}

impl<T: SerdeType> fmt::Debug for BlockWriter<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("BlockWriter")
            .field("buffer_size", &self.buffer.records.len())
            .field("buffer_start", &self.buffer.start)
            .field("buffer_end", &self.buffer.end)
            .field("config", &self.config)
            .finish()
    }
}

impl<T: SerdeType> BlockWriter<T> {
    /// Create a new [`BlockWriter`].
    pub fn new(config: BlockConfig) -> io::Result<Self> {
        Ok(Self {
            buffer: BlockBuffer {
                records: Vec::new(),
                start: u32::MAX,
                end: 0,
            },
            config,
        })
    }

    /// Check if the [`BlockWriter`] should be flushed (e.g. the
    /// current block serialized to binary, compressed, and written).
    /// The conditions for flush are that the buffer size is >= max_records.
    fn should_flush(&self) -> bool {
        if self.buffer.records.is_empty() {
            return false;
        }

        self.buffer.records.len() >= self.config.max_records
    }

    /// Flush the current [`BlockWriter`], by serializing
    /// the records in the block to binary and compressing them.
    pub fn flush(&mut self) -> io::Result<CompressedBlock> {
        // Serialize records
        let records_data = bincode::serialize(&self.buffer.records)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

        // Compress into a new Vec
        let mut compressed_data = Vec::new();
        zstd::stream::copy_encode(
            records_data.as_slice(),
            &mut compressed_data,
            self.config.level,
        )?;

        // Create block header with the compressed data
        let block = CompressedBlock {
            compressed_size: compressed_data.len() as u32,
            uncompressed_size: records_data.len() as u32,
            n_records: self.buffer.records.len() as u32,
            data: compressed_data,
            start: self.buffer.start,
            end: self.buffer.end,
        };

        // Clear buffer
        self.buffer.records.clear();
        self.buffer.start = u32::MAX;
        self.buffer.end = 0;

        Ok(block)
    }

    /// Add a record to the [`BlockWriter`] buffer.
    pub fn add_record(
        &mut self,
        start: u32,
        end: u32,
        record: T,
    ) -> io::Result<Option<CompressedBlock>> {
        self.buffer.start = self.buffer.start.min(start);
        self.buffer.end = self.buffer.end.max(end);
        self.buffer.records.push((start, end, record));

        if self.should_flush() {
            // No need to pass start
            Ok(Some(self.flush()?))
        } else {
            Ok(None)
        }
    }

    /// Flush remaining records if any exist
    pub fn finish(&mut self) -> io::Result<Option<CompressedBlock>> {
        if self.buffer.records.is_empty() {
            Ok(None)
        } else {
            Ok(Some(self.flush()?))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde::{Deserialize, Serialize};

    #[derive(Debug, Serialize, Deserialize, PartialEq)]
    struct TestRecord {
        value: i32,
    }

    #[test]
    fn test_block_flush() {
        let mut writer = BlockWriter::new(BlockConfig {
            max_records: 2,
            level: 3,
        })
        .unwrap();

        // println!("Adding first record");
        // First record shouldn't trigger flush
        let result = writer
            .add_record(1000, 2000, TestRecord { value: 1 })
            .unwrap();
        assert!(result.is_none());

        // println!("Adding second record");
        // Second record should trigger flush
        let result = writer
            .add_record(2000, 3000, TestRecord { value: 2 })
            .unwrap();
        assert!(result.is_some());

        let block = result.unwrap();
        assert_eq!(block.n_records, 2);
        assert_eq!(block.start, 1000);
        assert_eq!(block.end, 3000);
    }

    #[test]
    fn test_large_block() {
        let max_records = 5000;
        let mut writer = BlockWriter::new(BlockConfig {
            max_records,
            level: 3,
        })
        .unwrap();

        // Add records until we trigger a flush
        for i in 0..max_records {
            // println!("Adding record {}", i);
            let start: u32 = i as u32;
            let end: u32 = i as u32 + 1000;
            let result = writer
                .add_record(start, end, TestRecord { value: i as i32 })
                .unwrap();

            if i < max_records - 1 {
                assert!(
                    result.is_none(),
                    "Shouldn't flush before reaching max_records"
                );
            } else {
                assert!(result.is_some(), "Should flush on last record");
                let block = result.unwrap();
                assert_eq!(block.n_records, max_records as u32);
                assert_eq!(block.start, 0);
                assert_eq!(block.end, i as u32 + 1000);
            }
        }
    }
}
