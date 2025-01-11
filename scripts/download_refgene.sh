#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Source the install script to get DuckDB if needed
source "${SCRIPT_DIR}/install_duckdb.sh"

# Create data directory
mkdir -p "${SCRIPT_DIR}/../tests/data"

# Run the SQL script using the local DuckDB installation
"${SCRIPT_DIR}/duckdb/duckdb" < "${SCRIPT_DIR}/download_refgene.sql"

echo -e "${GREEN}RefGene data downloaded successfully${NC}"
