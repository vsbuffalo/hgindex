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

        // Check if we have enough room for the ENTIRE length prefix
        if self.current_block_size + super::LENGTH_PREFIX_SIZE > BLOCK_SIZE {
            self.flush_block()?;
        }

        // Save starting position after potential flush but before writing
        let start_vo = VirtualOffset::new(self.current_file_offset, self.current_block_size as u16);

        // Now we know the full length prefix will fit in this block
        let length_bytes = length.to_le_bytes();
        self.buffer.extend_from_slice(&length_bytes);
        self.current_block_size += length_bytes.len();

        // Write record data (this can span blocks as needed)
        let mut bytes_written = 0;
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
