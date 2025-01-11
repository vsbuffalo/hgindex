// block/reader.rs
//
// This processes records in batches, which seems to
// slightly improve performance.

use std::{
    fs::File,
    io::{self, BufReader, Read, Seek, SeekFrom},
    marker::PhantomData,
};

use crate::{error::HgIndexError, GenomicCoordinates};
use serde::Deserialize;

use super::{VirtualOffset, BATCH_SIZE, BLOCK_SIZE};

#[derive(Debug)]
pub struct BlockReader<T>
where
    T: GenomicCoordinates + std::fmt::Debug,
{
    inner: BufReader<File>,
    decomp_buf: Vec<u8>,
    record_batch: Vec<Vec<u8>>, // Pre-allocated batch buffer
    _phantom: PhantomData<T>,
}

impl<T> BlockReader<T>
where
    T: GenomicCoordinates + for<'de> Deserialize<'de> + std::fmt::Debug,
{
    pub fn new(file: File) -> io::Result<Self> {
        Ok(Self {
            inner: BufReader::with_capacity(BLOCK_SIZE, file),
            decomp_buf: Vec::with_capacity(BLOCK_SIZE),
            record_batch: Vec::with_capacity(BATCH_SIZE),
            _phantom: PhantomData,
        })
    }

    pub fn read_records_between(
        &mut self,
        min_offset: VirtualOffset,
        max_offset: VirtualOffset,
        query_start: u32,
        query_end: u32,
    ) -> Result<Vec<T>, HgIndexError> {
        let mut results = Vec::new();
        self.inner.seek(SeekFrom::Start(min_offset.file_offset()))?;
        let mut current_pos = min_offset.file_offset();

        // Create a reusable decoder to avoid allocation
        let mut decoder = zstd::bulk::Decompressor::new().unwrap();
        let mut compressed = Vec::with_capacity(BLOCK_SIZE);

        while current_pos <= max_offset.file_offset() {
            // Read block header
            let mut header = [0u8; 8];
            match self.inner.read_exact(&mut header) {
                Ok(()) => (),
                Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(HgIndexError::IOError(e)),
            }

            let uncompressed_size = u32::from_le_bytes(header[0..4].try_into().unwrap());
            let compressed_size = u32::from_le_bytes(header[4..8].try_into().unwrap());

            // Reuse compressed buffer
            compressed.resize(compressed_size as usize, 0);
            self.inner.read_exact(&mut compressed)?;

            // Reuse decompression buffer
            self.decomp_buf.resize(uncompressed_size as usize, 0);

            // Use bulk decompression instead of streaming
            decoder
                .decompress_to_buffer(&compressed, &mut self.decomp_buf)
                .map_err(|e| HgIndexError::DecompressionError(e.to_string()))?;

            let mut block_offset = 0;
            self.record_batch.clear();

            // Collect records into batch
            while block_offset + 8 <= self.decomp_buf.len() {
                let record_len = u64::from_le_bytes(
                    self.decomp_buf[block_offset..block_offset + 8]
                        .try_into()
                        .unwrap(),
                );
                block_offset += 8;

                if block_offset + record_len as usize > self.decomp_buf.len() {
                    break;
                }

                self.record_batch.push(
                    self.decomp_buf[block_offset..block_offset + record_len as usize].to_vec(),
                );

                if self.record_batch.len() >= BATCH_SIZE {
                    process_record_batch(&self.record_batch, query_start, query_end, &mut results)?;
                    self.record_batch.clear();
                }

                block_offset += record_len as usize;
            }

            // Process any remaining records in batch
            if !self.record_batch.is_empty() {
                process_record_batch(&self.record_batch, query_start, query_end, &mut results)?;
            }

            current_pos += 8 + compressed_size as u64;
        }

        Ok(results)
    }
}

// Keep the existing process_record_batch function
fn process_record_batch<T>(
    batch: &[Vec<u8>],
    query_start: u32,
    query_end: u32,
    results: &mut Vec<T>,
) -> Result<(), HgIndexError>
where
    T: GenomicCoordinates + for<'de> Deserialize<'de>,
{
    for record_data in batch {
        if let Ok(record) = bincode::deserialize(record_data) {
            let record: T = record;
            if record.start() < query_end && record.end() > query_start {
                results.push(record);
            }
        }
    }
    Ok(())
}
