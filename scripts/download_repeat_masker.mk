.PHONY: all clean download compress split

SCRIPT_DIR := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DATA_DIR := $(SCRIPT_DIR)/../tests/data
DUCKDB_DIR := $(SCRIPT_DIR)/duckdb
DUCKDB := $(DUCKDB_DIR)/duckdb
RM_BED := $(DATA_DIR)/repeat_masker.bed
RM_BGZ := $(DATA_DIR)/repeat_masker.bed.bgz
RM_TBI := $(DATA_DIR)/repeat_masker.bed.bgz.tbi
CHROM_DIR := $(DATA_DIR)/repeat_masker_by_chrom
CHROMS := $(patsubst %,$(CHROM_DIR)/chr%.bed,$(shell seq 1 22))
SHELL := /bin/bash

$(DUCKDB_DIR):
	mkdir -p $@

$(DUCKDB): | $(DUCKDB_DIR)
	bash $(SCRIPT_DIR)/install_duckdb.sh

all: $(RM_TBI) split

$(DATA_DIR) $(CHROM_DIR):
	mkdir -p $@

$(RM_BED): $(DUCKDB) | $(DATA_DIR)
	$< < $(SCRIPT_DIR)/repeat_masker.sql

$(RM_BGZ): $(RM_BED)
	bgzip -c $< > $@

$(RM_TBI): $(RM_BGZ)
	tabix -p bed $<

# Rule to create individual chromosome files in parallel
$(CHROM_DIR)/chr%.bed: $(RM_BED) | $(CHROM_DIR)
	awk -v chr="chr$*" '$$1 == chr' $< > $@

# Split target that can be run in parallel with -j
split: $(CHROMS)

clean:
	rm -f $(RM_BED) $(RM_BGZ) $(RM_TBI)
	rm -rf $(CHROM_DIR)
