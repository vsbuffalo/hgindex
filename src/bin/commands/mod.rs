// bin/commands/mod.rs

#[cfg(feature = "cli")]
pub mod pack;
#[cfg(feature = "cli")]
pub mod query;
#[cfg(all(feature = "cli", feature = "dev"))]
pub mod random_bed;

use hgindex::GenomicCoordinates;
use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Debug, Serialize, Deserialize)]
pub struct BedRecord {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub rest: String, // Store remaining fields as a single string
}

impl GenomicCoordinates for BedRecord {
    fn start(&self) -> u32 {
        self.start
    }

    fn end(&self) -> u32 {
        self.end
    }
}

impl fmt::Display for BedRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Format the record as a TSV (tab-separated values) line
        write!(
            f,
            "{}\t{}\t{}\t{}",
            self.chrom, self.start, self.end, self.rest
        )
    }
}
