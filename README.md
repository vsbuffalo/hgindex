[![CI](https://github.com/vsbuffalo/hgindex/actions/workflows/ci.yml/badge.svg)](https://github.com/vsbuffalo/hgindex/actions/workflows/ci.yml)

# hgindex — hierarchical genome index (& binary store!)

A *beta* flexible genomic range binning index with optional indexed binary
serialization of arbitrary types (using
[bincode](https://github.com/bincode-org/bincode)) implementation in Rust.

This is based on the UCSC genome browser's hierarchical binning scheme (also
used by Heng Li's [tabix](https://pmc.ncbi.nlm.nih.gov/articles/PMC3042176/))
but allows for custom binning configurations.

**Warning**: the API and name of this project may change as it is developed.
Please give feedback, suggestions, and reach out if you'd like to contribute.

## Features

- UCSC-compatible hierarchical binning scheme and fast index-based overlap
  lookups.
- Configurable bin sizes and levels
- Memory-efficient linear index for fast region queries
- Binary serializable indices for persistent storage
- `GenomicDataStore` for serializing any data type with bincode and
  [serde](https://serde.rs)
- Support for optional index-level metadata storage (e.g. for storing mapping
  between chromosome names and indices, etc)

## Implementation Details

The index uses a hierarchical binning strategy where genomic ranges are
assigned to bins based on their size. The default configuration matches UCSC:

- 5 levels with 8x scaling between levels
- Base bin size of 128kb
- Level bin counts: 4096, 512, 64, 8, 1

Custom configurations available:
```rust
// Default UCSC scheme
let bins = HierarchicalBins::ucsc();  // 128kb base, 8x scaling, 5 levels

// Higher resolution scheme
let dense = HierarchicalBins::dense();  // 16kb base, 4x scaling, 6 levels

// Custom scheme
let custom = HierarchicalBins::new(
    17,  // base_shift (128kb bins)
    3,   // level_shift (8x scaling)
    5    // num_levels
);
```

## Usage

Basic usage with a custom record type:

```rust
use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize)]
struct Record {
    name: String,
    score: f64,
}

// Create new store
let mut store = GenomicDataStore::<Record, ()>::create("data/", None)?;

// Add records
store.add_record("chr1", 1000, 2000, &Record {
    name: "feature1".into(),
    score: 0.5,
})?;

// Finalize and save
store.finalize()?;

// Query overlapping features
let mut store = GenomicDataStore::<Record, ()>::open("data/", None)?;
let results = store.get_overlapping("chr1", 1500, 2500);
```

With metadata:

```rust
#[derive(Serialize, Deserialize)]
struct Metadata {
    genome_build: String,
    creation_date: String,
}

let mut store = GenomicDataStore::<Record, Metadata>::create("data/", None)?;
store.set_metadata(Metadata {
    genome_build: "GRCh38".into(),
    creation_date: "2024-03-20".into(),
});
```

## Performance

Lookups are quite fast, but because bincode isn't as efficient as packed C
structures, there is some disk overhead.

