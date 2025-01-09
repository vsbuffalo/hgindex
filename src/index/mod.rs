// index/mod.rs
mod binning;
mod binning_index;

pub use binning::{BinningSchema, HierarchicalBins};
pub use binning_index::{BinningIndex, Chunk, Feature, SequenceIndex};
