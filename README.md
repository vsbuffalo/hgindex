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
$  bash scripts/benchmark_query.sh
Benchmark 1: ./target/release/hgidx query \
    --regions tests/data/refgene.bed.gz \
    --input tests/data/repeat_masker_autosomes.hgidx \
    > /dev/null
  Time (mean ± σ):     486.6 ms ±   3.8 ms    [User: 444.6 ms, System: 36.7 ms]
  Range (min … max):   480.8 ms … 493.5 ms    20 runs

Benchmark 2: tabix \
    tests/data/repeat_masker_autosomes.bed.bgz \
    --regions tests/data/refgene.bed.gz \
    > /dev/null
  Time (mean ± σ):      2.708 s ±  0.028 s    [User: 2.627 s, System: 0.059 s]
  Range (min … max):    2.679 s …  2.762 s    20 runs

Summary
    [hgidx] ran 5.56 ± 0.07 times faster than tabix

$ bash scripts/compare_sizes.sh
comparison of data sizes (ratios included):
-------------------------------------------
tabix file (.bgz):              144M     1x
hgidx bins:                     358M	 2.47x
tabix index (.tbi):             560K     1x
hgidx index (index.bin):        118M     215.96x

```

So presently hgindex is about 5.6x faster than tabix, with about 2.5x the disk
for the binary data. The index, since all ranges are stored in it (rather than
just block offset ranges), is quite a bit larger (216x) but still small
enough to load into memory. While index files can get quite large, they are still
small enough to load in memory (e.g. 100M records of start/end as u32, 8 bytes
each, would be roughly ≈800MB).

**Warning**: the API and name of this project may change as it is developed.
Please give feedback, suggestions, and reach out if you'd like to contribute.

## Features


- Zero cost deserialization and memory-mapped, so very optimized queries.
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


