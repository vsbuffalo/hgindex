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
pub struct BlockReader<T>
where
    T: std::fmt::Debug,
{
    inner: BufReader<File>,
    current_block: Vec<u8>,
    current_block_offset: u64,
    last_virtual_offset: Option<VirtualOffset>,
    _phantom: PhantomData<T>,
}

impl<T> BlockReader<T>
where
    T: for<'de> Deserialize<'de> + std::fmt::Debug,
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
        let mut result = Vec::with_capacity(length as usize);
        let mut bytes_remaining = length;
        let mut current_voffset = voffset;

        // Check if we're too close to block boundary
        let block_remaining = BLOCK_SIZE - current_voffset.block_offset() as usize;
        if block_remaining < LENGTH_PREFIX_SIZE && bytes_remaining > block_remaining as u64 {
            // Move to next block immediately
            self.inner
                .seek(SeekFrom::Start(current_voffset.file_offset() + 4))?;
            let mut comp_size_bytes = [0u8; 4];
            self.inner.read_exact(&mut comp_size_bytes)?;
            let comp_size = u32::from_le_bytes(comp_size_bytes);

            // Move to start of next block
            current_voffset =
                VirtualOffset::new(current_voffset.file_offset() + 8 + comp_size as u64, 0);
        }

        while bytes_remaining > 0 {
            // Seek to the appropriate block
            self.seek_virtual(current_voffset)?;

            let block_offset = current_voffset.block_offset() as usize;
            if block_offset >= self.current_block.len() {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    format!(
                        "Block offset {} beyond block size {}",
                        block_offset,
                        self.current_block.len()
                    ),
                ));
            }

            // Calculate how many bytes we can read from this block
            let available_in_block = self.current_block.len() - block_offset;
            let bytes_to_read = bytes_remaining.min(available_in_block as u64);

            // Read what we can from this block
            result.extend_from_slice(
                &self.current_block[block_offset..block_offset + bytes_to_read as usize],
            );
            bytes_remaining -= bytes_to_read;

            if bytes_remaining > 0 {
                // Need to read from next block
                self.inner
                    .seek(SeekFrom::Start(current_voffset.file_offset() + 4))?;
                let mut comp_size_bytes = [0u8; 4];
                self.inner.read_exact(&mut comp_size_bytes)?;
                let comp_size = u32::from_le_bytes(comp_size_bytes);

                // Move to start of next block
                current_voffset =
                    VirtualOffset::new(current_voffset.file_offset() + 8 + comp_size as u64, 0);
            }
        }

        Ok(result)
    }

    /// Read a single record at the provided virtual offset (even if it spans multiple
    /// compressed blocks).
    pub fn read_record(&mut self, voffset: VirtualOffset) -> Result<(u32, u32, T), HgIndexError>
    where
        T: for<'de> serde::Deserialize<'de> + std::fmt::Debug,
    {
        // Check if we're too close to block boundary
        let block_remaining = BLOCK_SIZE - voffset.block_offset() as usize;
        if block_remaining < LENGTH_PREFIX_SIZE {
            return Err(HgIndexError::InvalidOffset(format!(
                "Record starts too close to block boundary: only {} bytes remaining in block",
                block_remaining
            )));
        }

        // First read length prefix
        let prefix_bytes = self.read_bytes(voffset, LENGTH_PREFIX_SIZE as u64)?;
        let record_length = u64::from_le_bytes(prefix_bytes.try_into().map_err(|_| {
            HgIndexError::DeserializationError("Failed to read length prefix".to_string())
        })?);

        // Now read the actual record data starting immediately after the length prefix
        let record_voffset = VirtualOffset::new(
            voffset.file_offset(),
            voffset.block_offset() + LENGTH_PREFIX_SIZE as u16,
        );

        let record_bytes = self.read_bytes(record_voffset, record_length)?;

        // // Debug logging - TODO remove
        // if record_bytes.len() >= 4 && record_bytes[0] != 4 {
        //     eprintln!("WARNING: Unexpected byte pattern at start of record:");
        //     eprintln!(
        //         "Virtual offset: file={}, block={}",
        //         voffset.file_offset(),
        //         voffset.block_offset()
        //     );
        //     eprintln!("Record length: {}", record_length);
        //     eprintln!(
        //         "First 8 bytes: {:?}",
        //         &record_bytes[..8.min(record_bytes.len())]
        //     );
        //     eprintln!("All bytes: {:?}", record_bytes);
        //     if let Ok(s) = std::str::from_utf8(&record_bytes) {
        //         eprintln!("As string: {:?}", s);
        //     }
        // }

        // Deserialize
        bincode::deserialize(&record_bytes).map_err(|e| {
            HgIndexError::DeserializationError(format!(
                "Failed to deserialize record: {}. Virtual offset: file={}, block={}",
                e,
                voffset.file_offset(),
                voffset.block_offset()
            ))
        })
    }
}
