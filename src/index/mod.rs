// index/mod.rs
mod binning;
mod binning_index;

pub use binning::HierarchicalBins;
pub use binning_index::{BinningIndex, Feature, SequenceIndex};
