INSTALL mysql;

ATTACH 'host=genome-mysql.soe.ucsc.edu user=genome' AS ucsc (TYPE mysql);

COPY (
    SELECT
        chrom,           -- BED standard cols first
        txStart,         -- BED start
        txEnd,           -- BED end
        -- Additional useful fields
        name,           -- gene name/ID
        name2,          -- alternative name
        strand,
        cdsStart,
        cdsEnd,
        exonCount,
        exonStarts,
        exonEnds,
        score
    FROM ucsc.hg38.refGene
    ORDER BY chrom, txStart  -- Sort by chromosome and start position
)
TO 'tests/data/refgene.bed'
WITH (HEADER false, DELIMITER E'\t');
