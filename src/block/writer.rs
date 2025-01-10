// block/writer.rs

use std::{
    fs::File,
    io::{self, BufWriter, Write},
};

use crate::error::HgIndexError;
use serde::Serialize;

use super::{VirtualOffset, BLOCK_SIZE, COMPRESSION_LEVEL};

#[derive(Debug)]
pub struct BlockWriter {
    inner: BufWriter<File>,
    current_block: Vec<u8>,
    current_file_offset: u64,
}

impl BlockWriter {
    pub fn new(file: File) -> io::Result<Self> {
        Ok(Self {
            inner: BufWriter::new(file),
            current_block: Vec::with_capacity(BLOCK_SIZE),
            current_file_offset: 0,
        })
    }

    pub fn add_record<T>(&mut self, record: &T) -> Result<VirtualOffset, HgIndexError>
    where
        T: Serialize + std::fmt::Debug,
    {
        // eprintln!("Serializing record: {:?}", record);
        let record_data = bincode::serialize(record)
            .map_err(|e| HgIndexError::SerializationError(e.to_string()))?;
        // eprintln!(
        // "Serialized data length: {}, first few bytes: {:?}",
        // record_data.len(),
        // &record_data[..std::cmp::min(16, record_data.len())]
        // );

        // Calculate the total space needed including the length prefix
        let total_size = record_data.len() + 8; // 8 bytes for length prefix

        // If this record won't fit in the current block, flush the block first
        if self.current_block.len() + total_size > BLOCK_SIZE {
            self.flush_block()?;
        }

        // If the record is larger than block size, write it in its own block
        if total_size > BLOCK_SIZE {
            // Flush any existing data
            self.flush_block()?;

            // Create special block just for this record
            let start_vo = VirtualOffset::new(self.current_file_offset, 0);
            let len = record_data.len() as u64;

            // Write length and data
            let mut block = Vec::with_capacity(total_size);
            block.extend_from_slice(&len.to_le_bytes());
            block.extend_from_slice(&record_data);

            // Debug logging for block writing
            //eprintln!(
            //    "Writing large record block. Total size: {}, Record size: {}",
            //    block.len(),
            //    record_data.len()
            //);
            //eprintln!(
            //    "First few bytes of block: {:?}",
            //    &block[..std::cmp::min(16, block.len())]
            //);

            // Compress and write the block
            let compressed = zstd::encode_all(&*block, COMPRESSION_LEVEL)
                .map_err(|e| HgIndexError::IOError(io::Error::new(io::ErrorKind::Other, e)))?;

            // Write header
            let uncomp_size = block.len() as u32;
            let comp_size = compressed.len() as u32;
            self.inner.write_all(&uncomp_size.to_le_bytes())?;
            self.inner.write_all(&comp_size.to_le_bytes())?;
            self.inner.write_all(&compressed)?;

            self.current_file_offset += (8 + compressed.len()) as u64;
            return Ok(start_vo);
        }

        // Normal case: record fits in a block
        let block_offset = self.current_block.len() as u16;
        let start_vo = VirtualOffset::new(self.current_file_offset, block_offset);

        // Write length prefix and record data
        let len = record_data.len() as u64;
        self.current_block.extend_from_slice(&len.to_le_bytes());
        self.current_block.extend_from_slice(&record_data);

        Ok(start_vo)
    }

    pub fn flush_block(&mut self) -> Result<(), HgIndexError> {
        if self.current_block.is_empty() {
            return Ok(());
        }

        let compressed = zstd::encode_all(&*self.current_block, COMPRESSION_LEVEL)
            .map_err(|e| HgIndexError::IOError(io::Error::new(io::ErrorKind::Other, e)))?;

        // Write block header
        let uncomp_size = self.current_block.len() as u32;
        let comp_size = compressed.len() as u32;
        self.inner.write_all(&uncomp_size.to_le_bytes())?;
        self.inner.write_all(&comp_size.to_le_bytes())?;
        self.inner.write_all(&compressed)?;

        self.current_file_offset += (8 + compressed.len()) as u64;
        self.current_block.clear();
        Ok(())
    }

    pub fn finish(mut self) -> Result<(), HgIndexError> {
        self.flush_block()?;
        Ok(self.inner.flush()?)
    }
}
