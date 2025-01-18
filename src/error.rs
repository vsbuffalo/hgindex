// error.rs

#[cfg(feature = "cli")]
use indicatif::style::TemplateError;
use std::num::ParseIntError;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum HgIndexError {
    #[error("Invalid interval: end ({end}) must be greater than start ({start})")]
    InvalidInterval { start: u32, end: u32 },

    #[error("GenomicDataStore has already been finalized.")]
    AlreadyFinalized,

    #[error("Unsorted features in bin {bin_id} of sequence {chrom}. Found position {current} after {previous}")]
    UnsortedFeatures {
        chrom: String,
        bin_id: u32,
        previous: u32,
        current: u32,
    },

    #[error("IO error: {0}")]
    IOError(#[from] std::io::Error),

    #[error("Invalid record: zero-length range [{0}, {1})")]
    ZeroLengthFeature(u32, u32),

    #[error("Serialization error: {0}")]
    SerializationError(String),

    #[error("Decompression error: {0}")]
    DecompressionError(String),

    #[error("Deserialization error: {0}")]
    DeserializationError(String),

    #[error("Invalid offset error: {0}")]
    InvalidOffset(String),

    #[error("Parse integer error: {0}")]
    ParseIntError(#[from] ParseIntError),

    #[error("Internal error: {0}")]
    BoxError(#[from] Box<dyn std::error::Error>),

    #[error("{0}")]
    StringError(String),

    #[cfg(feature = "cli")]
    #[error("CSV error: {0}")]
    CsvError(#[from] csv::Error),

    #[cfg(feature = "cli")]
    #[error("Template error: {0}")]
    TemplateError(#[from] TemplateError),
}

// Add a convenience implementation for &str errors
impl From<&str> for HgIndexError {
    fn from(error: &str) -> Self {
        HgIndexError::StringError(error.to_string())
    }
}

// Add a convenience implementation for string errors
impl From<String> for HgIndexError {
    fn from(error: String) -> Self {
        HgIndexError::StringError(error)
    }
}

#[cfg(feature = "cli")]
impl From<crate::io::IoError> for HgIndexError {
    fn from(error: crate::io::IoError) -> Self {
        match error {
            crate::io::IoError::IoError(e) => HgIndexError::IOError(e),
            crate::io::IoError::InvalidGzipHeader => {
                HgIndexError::StringError("Invalid gzip header".into())
            }
        }
    }
}
