#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Source the install script to get DuckDB if needed
source "${SCRIPT_DIR}/install_duckdb.sh"

# Create data directory
mkdir -p "${SCRIPT_DIR}/../tests/data"

# Run the SQL script using the local DuckDB installation
# This creates the compressed TSV files from UCSC SQL query.
"${SCRIPT_DIR}/duckdb/duckdb" < "${SCRIPT_DIR}/repeat_masker.sql"

# Since we're already outputting gzipped files, we just need to create bgzip versions
echo "Preparing test data..."
gunzip -c "${SCRIPT_DIR}/../tests/data/repeat_masker.bed.gz" | bgzip -c > "${SCRIPT_DIR}/../tests/data/repeat_masker.bed.bgz"
tabix -p bed "${SCRIPT_DIR}/../tests/data/repeat_masker.bed.bgz"

echo "Repeat masker data downloaded successfully"
