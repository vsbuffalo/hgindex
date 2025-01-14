// bin/commands/mod.rs

#[cfg(feature = "cli")]
pub mod pack;
#[cfg(feature = "cli")]
pub mod query;
#[cfg(all(feature = "cli", feature = "dev"))]
pub mod random_bed;
