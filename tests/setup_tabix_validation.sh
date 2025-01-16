#!/bin/bash
# This is a setup like that in tabix_validation.rs for manual comparison
# (i.e. if tests fail).

set -e  # Exit on error

# Check required tools
bgzip --version || { echo "bgzip not found"; exit 1; }
tabix --version || { echo "tabix not found"; exit 1; }

# Set up directories
TEST_DIR="target/test_files"
mkdir -p "$TEST_DIR"
TEST_BED="$TEST_DIR/test.bed"
BGZIPPED="$TEST_DIR/test.bed.bgz"
HGINDEX="$TEST_DIR/test.hgidx"

echo "Test files location:"
echo "BED file: $TEST_BED"
echo "BGZip file: $BGZIPPED"
echo "HGIndex file: $HGINDEX"

# Generate random BED file
echo "Generating test BED file..."
cargo run --release --features=cli,dev -- random-bed -o "$TEST_BED" -n 1000000

# Compress with BGZip
echo "Running BGZip compression..."
rm -f "$BGZIPPED"  # Remove if exists
bgzip -c "$TEST_BED" > "$BGZIPPED"

# Create tabix index
echo "Running Tabix indexing..."
tabix -p bed -c "#" "$BGZIPPED"

# Create HGIndex
echo "Running HGIndex pack..."
cargo run --release --features=cli,dev -- pack --output "$HGINDEX" --force  $TEST_BED

# Example regions to test
echo "Now you can test queries like:"
echo "tabix $BGZIPPED chr1:15000-17000"
echo "cargo run --release --features=cli,dev -- query -i $HGINDEX chr1:15000-17000"
