[![CI](https://github.com/vsbuffalo/hgindex/actions/workflows/ci.yml/badge.svg)](https://github.com/vsbuffalo/hgindex/actions/workflows/ci.yml)

# hgindex — hierarchical genome index (& binary store!)

A *beta* flexible genomic range binning index with *uncompressed* indexed
binary serialization of arbitrary types (using
[bincode](https://github.com/bincode-org/bincode)) implementation in Rust.
Tabix-like block compression is the branch `block-compression`, but is
currently quite slow for the particular task I've needed this for, and still
experimental.

This binning is based on the UCSC genome browser's hierarchical binning scheme
(also used by Heng Li's
[tabix](https://pmc.ncbi.nlm.nih.gov/articles/PMC3042176/)) but allows for
custom binning configurations. However, the full index (i.e. all range, rather
than ranges of virtual offsets) is currently used now (the `block-compression`
uses the tabix-like offset range though).

The current version is not block compressed, which allows for direct memory
mapping. This is on the other side of the disk-space / query-speed tradeoff
than tabix. Here are some preliminary query benchmarks, using the UCSC repeat
masker track as the database, and the UCSC RefGene track as the query:

```
$ hyperfine --warmup 3 \
  './target/release/hgidx  query --regions tests/data/refgene.bed --input tests/data/repeat_masker.hgidx > /dev/null' \
  'tabix tests/data/repeat_masker.bed.bgz --regions tests/data/refgene.bed > /dev/null'
    Finished `release` profile [optimized] target(s) in 0.19s
Benchmark 1: ./target/release/hgidx  query --regions tests/data/refgene.bed --input tests/data/repeat_masker.hgidx > /dev/null
  Time (mean ± σ):      2.158 s ±  0.015 s    [User: 2.047 s, System: 0.082 s]
  Range (min … max):    2.140 s …  2.181 s    10 runs

Benchmark 2: tabix tests/data/repeat_masker.bed.bgz --regions tests/data/refgene.bed > /dev/null
  Time (mean ± σ):      3.920 s ±  0.035 s    [User: 3.785 s, System: 0.091 s]
  Range (min … max):    3.876 s …  3.991 s    10 runs

Summary
  ./target/release/hgidx  query --regions tests/data/refgene.bed --input tests/data/repeat_masker.hgidx > /dev/null ran
    1.82 ± 0.02 times faster than tabix tests/data/repeat_masker.bed.bgz --regions tests/data/refgene.bed > /dev/null

$ du -h tests/data/repeat_masker.hgidx
622M tests/data/repeat_masker.hgidx

$ du -h tests/data/repeat_masker.bed.bgz
144M     tests/data/repeat_masker.bed.bgz
```

So presently hgindex is about 1.8x faster than tabix, with about 4.3x the disk
usage.

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
- Simple command line tool for working with BED-like files (mostly for testing/benchmarks)

## Command line tool

There is a simple, optional command line tool (use `--features=cli` when
installing to enable), primarily for integration/validation tests and
benchmarks against tabix's results.

```
$ hgidx pack tests/data/repeat_masker.bed        # serialie the BED records and index
$ hgidx  query --regions tests/data/refgene.bed    # query a set of ranges
$ hgidx  query chr2:4131233-4131233                # query a single range
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


