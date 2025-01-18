#!/bin/bash

# Paths
BGZ_FILE="tests/data/repeat_masker.bed.bgz"
TBI_FILE="tests/data/repeat_masker.bed.bgz.tbi"
HGIDX_DIR="tests/data/repeat_masker_autosomes.hgidx"

# Function to get file size in bytes (portable for macOS and Linux)
get_size_in_bytes() {
    if [[ "$OSTYPE" == "darwin"* ]]; then
        stat -f%z "$1"
    else
        stat -c%s "$1"
    fi
}

# Get sizes in bytes for calculations
BGZ_BYTES=$(get_size_in_bytes "$BGZ_FILE")
TBI_BYTES=$(get_size_in_bytes "$TBI_FILE")
HGIDX_BINS_BYTES=$(find "$HGIDX_DIR" -name 'chr*bin' -exec stat -f%z {} + | awk '{sum+=$1} END {print sum}')
HGIDX_INDEX_BYTES=$(get_size_in_bytes "$HGIDX_DIR/index.bin")

# Compute human-readable sizes for display
BGZ_SIZE=$(du -h "$BGZ_FILE" | awk '{print $1}')
TBI_SIZE=$(du -h "$TBI_FILE" | awk '{print $1}')
HGIDX_BINS_SIZE=$(du -ch "$HGIDX_DIR"/chr*bin | tail -n 1 | sed 's/total//')
HGIDX_INDEX_SIZE=$(du -h "$HGIDX_DIR/index.bin" | awk '{print $1}')

# Compute ratios
BINS_TO_TABIX_RATIO=$(echo "scale=2; $HGIDX_BINS_BYTES / $BGZ_BYTES" | bc)
INDEX_TO_TABIX_RATIO=$(echo "scale=2; $HGIDX_INDEX_BYTES / $TBI_BYTES" | bc)

# Output results
echo "comparison of data sizes (ratios included):"
echo "-------------------------------------------"
printf "%-25s %10s %6s\n" "tabix file (.bgz):" "$BGZ_SIZE" "1x"
printf "%-25s %11s %2s\n" "hgidx bins:" "$HGIDX_BINS_SIZE" "${BINS_TO_TABIX_RATIO}x"
printf "%-25s %10s %6s\n" "tabix index (.tbi):" "$TBI_SIZE" "1x"
printf "%-25s %10s %11s\n" "hgidx index (index.bin):" "$HGIDX_INDEX_SIZE" "${INDEX_TO_TABIX_RATIO}x"

