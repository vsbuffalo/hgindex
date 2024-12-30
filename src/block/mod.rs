/// block/mod.rs
use std::io::{self, Write};

/// A compressed block of records, complete with metadata.
#[derive(Debug)]
pub struct CompressedBlock {
    /// Size of compressed data in bytes
    pub compressed_size: u32,
    /// Size of uncompressed data in bytes
    pub uncompressed_size: u32,
    /// Number of records in block
    pub n_records: u32,
    /// The compressed data bytes
    pub data: Vec<u8>,
    /// Leftmost start position of records
    pub start: u32,
    /// Rightmost end position of records  
    pub end: u32,
}

/// Configuration for block parameters
#[derive(Debug, Clone)]
pub struct BlockConfig {
    /// Maximum records per block
    pub max_records: usize,
    /// zstd compression level
    pub level: i32,
}

impl Default for BlockConfig {
    fn default() -> Self {
        Self {
            max_records: 10_000,
            level: 3,
        }
    }
}

/// Write a block to the given writer
pub(crate) fn write_block<W: Write>(writer: &mut W, block: &CompressedBlock) -> io::Result<()> {
    // Write header fields in little-endian format
    writer.write_all(&block.compressed_size.to_le_bytes())?;
    writer.write_all(&block.uncompressed_size.to_le_bytes())?;
    writer.write_all(&block.n_records.to_le_bytes())?;
    writer.write_all(&block.start.to_le_bytes())?;
    writer.write_all(&block.end.to_le_bytes())?;

    // Write compressed data
    writer.write_all(&block.data)?;

    Ok(())
}

/// Buffer for accumulating records before compression
#[derive(Debug)]
pub struct BlockBuffer<T> {
    /// Records stored as (start, end, data) tuples
    pub records: Vec<(u32, u32, T)>,
    /// Leftmost start position of buffered records
    pub start: u32,
    /// Rightmost end position of buffered records
    pub end: u32,
}

pub mod reader;
pub mod writer;

pub use reader::BlockReader;
pub use writer::BlockWriter;
