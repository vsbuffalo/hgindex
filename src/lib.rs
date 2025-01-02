use serde::{Deserialize, Serialize};

pub mod error;
pub mod index;
pub mod store;

pub use index::{BinningIndex, Feature, HierarchicalBins, SequenceIndex};
pub use store::GenomicDataStore;

#[cfg(test)]
pub(crate) mod test_utils;

pub trait SerdeType: Serialize + for<'de> Deserialize<'de> {}
impl<T: Serialize + for<'de> Deserialize<'de>> SerdeType for T {}
