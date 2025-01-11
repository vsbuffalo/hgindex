#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BENCH_DIR="${SCRIPT_DIR}/../benches/results"
DATE=$(date +%Y%m%d_%H%M%S)
# Colors for output
GREEN='\033[0;32m'
NC='\033[0m' # No Color
# Create benchmark results directory
mkdir -p "${BENCH_DIR}"

# Make sure we have our test data
if [ ! -f "${SCRIPT_DIR}/../tests/data/repeat_masker.bed" ]; then
    echo "Downloading, bgzipping, and tabix indexing repeat masker data..."
    "${SCRIPT_DIR}/download_repeat_masker.sh"
fi
if [ ! -f "${SCRIPT_DIR}/../tests/data/refgene.bed" ]; then
    echo "Downloading refgene data..."
    "${SCRIPT_DIR}/download_refgene.sh"
fi

# Create chr2-only file if it doesn't exist
if [ ! -f "${SCRIPT_DIR}/../tests/data/refgene_chr2.bed" ]; then
    echo "Creating chr2-only refgene file..."
    grep "^chr2\t" "${SCRIPT_DIR}/../tests/data/refgene.bed" > "${SCRIPT_DIR}/../tests/data/refgene_chr2.bed"
fi

# Build hgidx index
echo "Building hgidx index..."
cargo build --release --features=cli
"${SCRIPT_DIR}/../target/release/hgidx" pack --force "${SCRIPT_DIR}/../tests/data/repeat_masker.bed" -o "${SCRIPT_DIR}/../tests/data/repeat_masker.hgidx"

# Run benchmarks for full data
echo "Running benchmarks for full data..."
RESULTS_FILE="${BENCH_DIR}/benchmark_full_${DATE}.txt"
echo "Benchmark Results (Full Data) - ${DATE}" > "${RESULTS_FILE}"
echo "----------------------" >> "${RESULTS_FILE}"

hyperfine --warmup 3 \
    --export-markdown "${BENCH_DIR}/benchmark_full_${DATE}.md" \
    --export-json "${BENCH_DIR}/benchmark_full_${DATE}.json" \
    "tabix ${SCRIPT_DIR}/../tests/data/repeat_masker.bed.bgz --regions ${SCRIPT_DIR}/../tests/data/refgene.bed > /dev/null" \
    "${SCRIPT_DIR}/../target/release/hgidx query --regions ${SCRIPT_DIR}/../tests/data/refgene.bed --input ${SCRIPT_DIR}/../tests/data/repeat_masker.hgidx > /dev/null"

# Run benchmarks for chr2 data
echo "Running benchmarks for chr2 data..."
RESULTS_FILE="${BENCH_DIR}/benchmark_chr2_${DATE}.txt"
echo "Benchmark Results (Chr2 Only) - ${DATE}" > "${RESULTS_FILE}"
echo "----------------------" >> "${RESULTS_FILE}"

hyperfine --warmup 3 \
    --export-markdown "${BENCH_DIR}/benchmark_chr2_${DATE}.md" \
    --export-json "${BENCH_DIR}/benchmark_chr2_${DATE}.json" \
    "tabix ${SCRIPT_DIR}/../tests/data/repeat_masker.bed.bgz --regions ${SCRIPT_DIR}/../tests/data/refgene_chr2.bed > /dev/null" \
    "${SCRIPT_DIR}/../target/release/hgidx query --regions ${SCRIPT_DIR}/../tests/data/refgene_chr2.bed --input ${SCRIPT_DIR}/../tests/data/repeat_masker.hgidx > /dev/null"

# Append file sizes to both results
for type in "full" "chr2"; do
    RESULTS_FILE="${BENCH_DIR}/benchmark_${type}_${DATE}.txt"
    echo -e "\nFile Sizes:" >> "${RESULTS_FILE}"
    echo "------------" >> "${RESULTS_FILE}"
    ls -lh "${SCRIPT_DIR}/../tests/data/repeat_masker.bed.bgz" >> "${RESULTS_FILE}"
    ls -lh "${SCRIPT_DIR}/../tests/data/repeat_masker.bed.bgz.tbi" >> "${RESULTS_FILE}"
    ls -lh "${SCRIPT_DIR}/../tests/data/repeat_masker.hgidx" >> "${RESULTS_FILE}"
    ls -lh "${SCRIPT_DIR}/../tests/data/refgene.bed" >> "${RESULTS_FILE}"
    ls -lh "${SCRIPT_DIR}/../tests/data/refgene_chr2.bed" >> "${RESULTS_FILE}"
done

echo -e "${GREEN}Benchmarks completed! Results saved to:${NC}"
echo "Full data results:"
echo "- ${BENCH_DIR}/benchmark_full_${DATE}.txt"
echo "- ${BENCH_DIR}/benchmark_full_${DATE}.md"
echo "- ${BENCH_DIR}/benchmark_full_${DATE}.json"
echo "Chr2 data results:"
echo "- ${BENCH_DIR}/benchmark_chr2_${DATE}.txt"
echo "- ${BENCH_DIR}/benchmark_chr2_${DATE}.md"
echo "- ${BENCH_DIR}/benchmark_chr2_${DATE}.json"
