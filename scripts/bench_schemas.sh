#!/bin/bash

set -euo pipefail

# Define variables
INPUT_BED="tests/data/repeat_masker_autosomes.bed"
REGIONS_BED="tests/data/refgene.bed"
BASE_DIR="tests/data/bench_schema"
SCHEMAS=("dense" "tabix" "ucsc")
HGIDX_BINARY="./target/release/hgidx"
TABIX_BINARY="tabix"
RESULTS_DIR="benchmark_results"
WARMUP=10
RUNS=20

# Ensure required directories exist
mkdir -p "$BASE_DIR" "$RESULTS_DIR"

# Pre-build the hgidx binary
cargo build --release --features=cli,dev

# Iterate over schemas and benchmark each one
for SCHEMA in "${SCHEMAS[@]}"; do
    SCHEMA_DIR="$BASE_DIR/$SCHEMA"
    HGIDX_FILE="$SCHEMA_DIR/repeat_masker_autosomes.hgidx"

    echo "Benchmarking schema: $SCHEMA"

    # Remove old index and create a new one
    rm -rf "$SCHEMA_DIR"
    mkdir -p "$SCHEMA_DIR"
    "$HGIDX_BINARY" pack "$INPUT_BED" --force --schema "$SCHEMA" --output "$HGIDX_FILE"

    # Run benchmarks and save results
    hyperfine --warmup "$WARMUP" --min-runs "$RUNS" \
        "$HGIDX_BINARY query --regions $REGIONS_BED --input $HGIDX_FILE > /dev/null" \
        "$TABIX_BINARY ${INPUT_BED}.bgz --regions $REGIONS_BED > /dev/null" \
        --export-json "$RESULTS_DIR/benchmark_${SCHEMA}.json" \
        --export-markdown "$RESULTS_DIR/benchmark_${SCHEMA}.md"
done

echo "Benchmarking complete. Results are saved in $RESULTS_DIR."

