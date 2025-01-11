// index/mod.rs
pub mod binning;
mod binning_index;

pub use binning::{BinningSchema, HierarchicalBins};
pub use binning_index::{BinningIndex, Feature, SequenceIndex};
