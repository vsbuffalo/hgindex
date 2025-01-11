#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Source the install script to get DuckDB if needed
source "${SCRIPT_DIR}/install_duckdb.sh"

# Create data directory
mkdir -p "${SCRIPT_DIR}/../tests/data"

# Run the SQL script using the local DuckDB installation
"${SCRIPT_DIR}/duckdb/duckdb" < "${SCRIPT_DIR}/repeat_masker.sql"

# Compress repeat masker with bgzip
if ! command -v bgzip &> /dev/null; then
    echo "bgzip not found. Please install htslib."
    exit 1
fi
echo "Preparing test data..."
bgzip -c "${SCRIPT_DIR}/../tests/data/repeat_masker.bed" > "${SCRIPT_DIR}/../tests/data/repeat_masker.bed.bgz"
tabix -p bed "${SCRIPT_DIR}/../tests/data/repeat_masker.bed.bgz"


echo -e "${GREEN}Repeat masker data downloaded successfully${NC}"
