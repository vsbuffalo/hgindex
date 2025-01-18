# Makefile for Rust genomic test data setup
SHELL := /bin/bash
SCRIPT_DIR := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
TEST_DATA_DIR := $(SCRIPT_DIR)/../tests/data
DUCKDB_DIR := $(SCRIPT_DIR)/duckdb
DUCKDB := $(DUCKDB_DIR)/duckdb
SQL_DIR := $(SCRIPT_DIR)/sql
HGIDX := $(SCRIPT_DIR)/../target/release/hgidx

# Color codes for output
GREEN := \033[0;32m
NC := \033[0m

# Define the data types we're working with
DATA_TYPES := refgene repeat_masker repeat_masker_autosomes

# Preserve intermediate files
.SECONDARY:

# Main targets
.PHONY: all clean clean-derived test-data check-test-data
all: test-data

# Test data preparation targets
test-data: genomic-test-data index-data

.PHONY: genomic-test-data index-data
genomic-test-data: $(foreach type,$(DATA_TYPES),$(TEST_DATA_DIR)/$(type).bed.bgz.tbi)
	@echo -e "$(GREEN)All genomic test data prepared successfully$(NC)"

# Build hgidx if needed
$(HGIDX):
	cd $(SCRIPT_DIR)/.. && cargo build --release --features=cli,dev

# Index creation rules
index-data: $(foreach type,$(DATA_TYPES),$(TEST_DATA_DIR)/$(type).hgidx)
	@echo -e "$(GREEN)All indices created successfully$(NC)"

$(TEST_DATA_DIR)/%.hgidx: $(TEST_DATA_DIR)/%.bed.gz $(HGIDX)
	$(HGIDX) pack $< --force

# Directory creation
$(TEST_DATA_DIR):
	mkdir -p $@

# Special rule for repeat_masker files since they're created together
$(TEST_DATA_DIR)/repeat_masker.bed.gz $(TEST_DATA_DIR)/repeat_masker_autosomes.bed.gz: $(SQL_DIR)/download_repeatmasker.sql $(DUCKDB) | $(TEST_DATA_DIR)
	cd $(SCRIPT_DIR)/.. && $(DUCKDB) < $<

# Generic pattern rules for other data types
$(TEST_DATA_DIR)/%.bed.gz: $(SQL_DIR)/download_%.sql $(DUCKDB) | $(TEST_DATA_DIR)
	cd $(SCRIPT_DIR)/.. && $(DUCKDB) < $<

$(TEST_DATA_DIR)/%.bed.bgz: $(TEST_DATA_DIR)/%.bed.gz
	gunzip -c $< | bgzip -c > $@

$(TEST_DATA_DIR)/%.bed.bgz.tbi: $(TEST_DATA_DIR)/%.bed.bgz
	tabix -p bed $<

# Clean targets
.PHONY: clean clean-derived clean-index clean-bgzip clean-all

# Clean only derived files (bgzip, tabix, and indices)
clean-derived: clean-index clean-bgzip
	@echo "Cleaned derived files (keeping SQL downloaded data)"

# Clean everything including SQL downloads
clean-all: clean-derived
	rm -f $(foreach type,$(DATA_TYPES),$(TEST_DATA_DIR)/$(type).bed.gz)
	@echo "Cleaned all data including SQL downloads"

# For backwards compatibility
clean: clean-derived

clean-bgzip:
	rm -f $(foreach type,$(DATA_TYPES),$(TEST_DATA_DIR)/$(type).bed.bgz)
	rm -f $(foreach type,$(DATA_TYPES),$(TEST_DATA_DIR)/$(type).bed.bgz.tbi)

clean-index:
	rm -rf $(foreach type,$(DATA_TYPES),$(TEST_DATA_DIR)/$(type).hgidx)
