// bin/commands/mod.rs

#[cfg(feature = "cli")]
pub mod pack;
#[cfg(feature = "cli")]
pub mod query;
#[cfg(all(feature = "cli", feature = "dev"))]
pub mod random_bed;

use hgindex::DataRecord;
use hgindex::GenomicCoordinates;
use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, PartialEq, Deserialize)]
pub struct BedRecord {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub rest: String, // Store remaining fields as a single string
}

impl DataRecord for BedRecord {
    fn write_tsv_line<W: std::io::Write>(&self, writer: &mut W) -> Result<(), std::io::Error> {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}",
            self.chrom, self.start, self.end, self.rest
        )
    }
}

impl GenomicCoordinates for BedRecord {
    fn start(&self) -> u32 {
        self.start
    }

    fn end(&self) -> u32 {
        self.end
    }
}
