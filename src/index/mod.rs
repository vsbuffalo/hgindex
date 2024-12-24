// index/mod.rs
mod binning;
mod index;

pub use binning::HierarchicalBins;
pub use index::{BinningIndex, Feature, SequenceIndex};

