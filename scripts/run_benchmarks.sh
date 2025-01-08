#!/bin/bash
set -e

mkdir -p data

# Ensure we have our test data
if [ ! -f "data/test.bed" ]; then
    echo "Generating test data..."
    cargo run --release --features=cli,dev -- random-bed -o data/test.bed -n 1000000
fi

# Create both compression versions
echo "Compressing files..."
if [ ! -f "data/test.bed.gz" ] || [ "data/test.bed" -nt "data/test.bed.gz" ]; then
    # Force standard gzip with -n flag
    cat data/test.bed | gzip > data/test.bed.gz
fi

if [ ! -f "data/test.bed.bgz" ] || [ "data/test.bed" -nt "data/test.bed.bgz" ]; then
    bgzip -c data/test.bed > data/test.bed.bgz
fi

# Index with both tools
echo "Indexing with tabix..."
rm -f data/test.bed.bgz.tbi
tabix -p bed -c '#' data/test.bed.bgz

echo "Indexing with hgindex..."
cargo run --release --features=cli -- pack -f data/test.bed -o data/test.bed.hgidx

# Update benchmark to use correct files
echo "Running benchmarks..."
# Make sure to update the benchmark code to use test.bed.bgz for tabix
cargo bench --features=cli
