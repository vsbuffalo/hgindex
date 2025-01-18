#!/bin/bash

hyperfine --warmup 10 --min-runs 20  './target/release/hgidx query \
    --regions tests/data/refgene.bed.gz \
    --input tests/data/repeat_masker_autosomes.hgidx \
    > /dev/null' \
\
  'tabix \
    tests/data/repeat_masker_autosomes.bed.bgz \
    --regions tests/data/refgene.bed.gz \
    > /dev/null'

