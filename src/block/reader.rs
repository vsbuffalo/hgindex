// block/reader.rs

use std::{
    fs::File,
    io::{self, BufReader, Read, Seek, SeekFrom},
    marker::PhantomData,
};

use serde::Deserialize;

use crate::error::HgIndexError;

use super::{VirtualOffset, BLOCK_SIZE, LENGTH_PREFIX_SIZE};

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
