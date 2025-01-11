INSTALL mysql;

ATTACH 'host=genome-mysql.soe.ucsc.edu user=genome' AS ucsc (TYPE mysql);

COPY (
    SELECT
        genoName AS chrom,  -- BED standard cols first
        genoStart,          -- BED start
        genoEnd,            -- BED end
        -- Additional fields after BED3
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
    FROM ucsc.hg38.rmsk
    ORDER BY genoName, genoStart  -- Sort by chromosome and start position
)
TO 'tests/data/repeat_masker.bed'
WITH (HEADER false, DELIMITER E'\t');
