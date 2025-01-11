#!/usr/bin/env bash
# Exit on error, undefined vars are errors, pipeline fails on first error
set -euo pipefail

# Get script directory for relative paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Version variable
DUCKDB_VERSION="v1.1.3"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

# Parse command line arguments
FORCE_INSTALL=0
while [[ $# -gt 0 ]]; do
    case $1 in
        -f|--force)
            FORCE_INSTALL=1
            shift
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            echo "Usage: $0 [-f|--force]"
            exit 1
            ;;
    esac
done

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to install DuckDB based on platform
install_duckdb() {
    echo "Installing DuckDB ${DUCKDB_VERSION}..."

    # Create duckdb directory in scripts if it doesn't exist
    mkdir -p "${SCRIPT_DIR}/duckdb"

    # Detect OS and architecture
    OS=$(uname -s | tr '[:upper:]' '[:lower:]')

    if [[ "$OS" == "darwin" ]]; then
        DUCKDB_URL="https://github.com/duckdb/duckdb/releases/download/${DUCKDB_VERSION}/duckdb_cli-osx-universal.zip"
    elif [[ "$OS" == "linux" ]]; then
        DUCKDB_URL="https://github.com/duckdb/duckdb/releases/download/${DUCKDB_VERSION}/duckdb_cli-linux-amd64.zip"
    else
        echo -e "${RED}Unsupported operating system${NC}"
        exit 1
    fi

    curl -L "$DUCKDB_URL" -o "${SCRIPT_DIR}/duckdb/duckdb.zip"
    unzip -o "${SCRIPT_DIR}/duckdb/duckdb.zip" -d "${SCRIPT_DIR}/duckdb"
    rm "${SCRIPT_DIR}/duckdb/duckdb.zip"
    chmod +x "${SCRIPT_DIR}/duckdb/duckdb"

    echo -e "${GREEN}DuckDB ${DUCKDB_VERSION} installed successfully${NC}"
}

# Main installation logic
if [ $FORCE_INSTALL -eq 1 ]; then
    echo "Force installing DuckDB..."
    install_duckdb
    export PATH="${SCRIPT_DIR}/duckdb:$PATH"
elif ! command_exists "${SCRIPT_DIR}/duckdb/duckdb"; then
    if [[ -x "${SCRIPT_DIR}/duckdb/duckdb" ]]; then
        export PATH="${SCRIPT_DIR}/duckdb:$PATH"
    else
        install_duckdb
        export PATH="${SCRIPT_DIR}/duckdb:$PATH"
    fi
fi
