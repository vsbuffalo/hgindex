use serde::{Deserialize, Serialize};

pub mod container;
pub mod index;

pub use container::GenomicDataStore;
pub use index::{BinningIndex, Feature, HierarchicalBins, SequenceIndex};

#[cfg(test)]
pub(crate) mod test_utils;

pub trait SerdeType: Serialize + for<'de> Deserialize<'de> {}
impl<T: Serialize + for<'de> Deserialize<'de>> SerdeType for T {}
