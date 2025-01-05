// error.rs

use thiserror::Error;

#[derive(Error, Debug)]
pub enum HgIndexError {
    #[error("Invalid interval: end ({end}) must be greater than start ({start})")]
    InvalidInterval { start: u32, end: u32 },

    #[error("IO error: {0}")]
    IOError(#[from] std::io::Error),

    #[error("Serialization error: {0}")]
    SerializationError(String),

    #[error("Deserialization error: {0}")]
    DeserializationError(String),
}
