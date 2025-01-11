use serde::{Deserialize, Serialize};

pub mod error;
pub mod index;
#[cfg(feature = "cli")]
pub mod io;
pub mod store;

pub use index::{BinningIndex, BinningSchema, Feature, HierarchicalBins, SequenceIndex};
#[cfg(feature = "cli")]
pub use io::*;
pub use store::GenomicDataStore;

#[cfg(test)]
pub(crate) mod test_utils;

pub trait SerdeType: Serialize + for<'de> Deserialize<'de> {}
impl<T: Serialize + for<'de> Deserialize<'de>> SerdeType for T {}

/// Trait for types that have genomic coordinates
pub trait GenomicCoordinates {
    /// Get the start coordinate (0-based, inclusive)
    fn start(&self) -> u32;

    /// Get the end coordinate (0-based, exclusive)
    fn end(&self) -> u32;
}
