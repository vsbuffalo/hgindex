// block/reader.rs

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

        // Pre-allocate batch buffers
        let mut record_batch = Vec::with_capacity(BATCH_SIZE);

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

            // Read compressed block
            let mut compressed = vec![0u8; compressed_size as usize];
            self.inner.read_exact(&mut compressed)?;

            // Set the size to the uncompressed size and zero out.
            self.decomp_buf.resize(uncompressed_size as usize, 0);
            let decompressed = &mut self.decomp_buf;

            let mut decoder = zstd::stream::read::Decoder::new(std::io::Cursor::new(&compressed))?;
            io::Read::read_exact(&mut decoder, &mut decompressed[..])?;

            let mut block_offset = 0;
            while block_offset + 8 <= decompressed.len() {
                // Read record length
                let record_len = u64::from_le_bytes(
                    decompressed[block_offset..block_offset + 8]
                        .try_into()
                        .unwrap(),
                );
                block_offset += 8;

                if block_offset + record_len as usize > decompressed.len() {
                    break;
                }

                // Add record data to batch
                record_batch
                    .push(decompressed[block_offset..block_offset + record_len as usize].to_vec());

                // Process batch when full
                if record_batch.len() >= BATCH_SIZE {
                    process_record_batch(&mut record_batch, query_start, query_end, &mut results)?;
                    record_batch.clear();
                }

                block_offset += record_len as usize;
            }

            // Process remaining records in batch
            if !record_batch.is_empty() {
                process_record_batch(&mut record_batch, query_start, query_end, &mut results)?;
                record_batch.clear();
            }

            current_pos += 8 + compressed_size as u64;
        }

        Ok(results)
    }
}

// Helper function to process a batch of records
fn process_record_batch<T>(
    batch: &mut [Vec<u8>],
    query_start: u32,
    query_end: u32,
    results: &mut Vec<T>,
) -> Result<(), HgIndexError>
where
    T: GenomicCoordinates + for<'de> Deserialize<'de>,
{
    for record_data in batch.iter() {
        let record: T = bincode::deserialize(record_data)
            .map_err(|e| HgIndexError::DeserializationError(e.to_string()))?;

        if record.start() < query_end && record.end() > query_start {
            results.push(record);
        }
    }
    Ok(())
}
