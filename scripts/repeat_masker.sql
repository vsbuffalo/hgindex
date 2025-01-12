-- First, install and attach MySQL
INSTALL mysql;
ATTACH 'host=genome-mysql.soe.ucsc.edu user=genome' AS ucsc (TYPE mysql);

-- Create a temporary table to store the data
CREATE TEMP TABLE rmsk_temp AS
SELECT
    genoName AS chrom,
    genoStart,
    genoEnd,
    repName,
    swScore,
    strand,
    repClass,
    repFamily,
    milliDiv,
    milliDel,
    milliIns,
    genoLeft,
    repStart,
    repEnd,
    repLeft,
    id
FROM ucsc.hg38.rmsk;

-- Write all chromosomes to first file
COPY (
    SELECT *
    FROM rmsk_temp
    ORDER BY chrom, genoStart
)
TO 'tests/data/repeat_masker.bed'
WITH (HEADER false, DELIMITER E'\t');

-- Write only chr1-chr22 to second file
-- This is for back-compatibility over a bug
-- with branch block-compression
COPY (
    SELECT *
    FROM rmsk_temp
    WHERE chrom IN (
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
        'chr20', 'chr21', 'chr22'
    )
    ORDER BY chrom, genoStart
)
TO 'tests/data/repeat_masker_autosomes.bed'
WITH (HEADER false, DELIMITER E'\t');
