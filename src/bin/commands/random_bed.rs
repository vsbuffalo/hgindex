// bin/commands/random_bed.rs

use clap::Args;
use hgindex::error::HgIndexError;
use hgindex::io::OutputStream;
use hgindex::BedRecord;
use rand::{seq::SliceRandom, Rng, SeedableRng};
use std::fmt::Write;
use std::path::PathBuf;

#[derive(Args)]
pub struct RandomBedArgs {
    /// Output file path (.bed or .bed.gz)
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Number of records to generate
    #[arg(short = 'n', long, default_value = "1000000")]
    pub num_records: usize,

    /// Optional seed for random number generation
    #[arg(short, long)]
    pub seed: Option<u64>,
}

pub fn run(args: RandomBedArgs) -> Result<(), HgIndexError> {
    eprintln!(
        "Generating {} random BED records to {}",
        args.num_records,
        args.output
            .as_ref()
            .map_or("<stdout>".to_string(), |v| v.to_string_lossy().to_string())
    );

    let output = OutputStream::new(args.output);
    let mut output_writer = output.writer()?;

    //// Write header comments
    //writeln!(
    //    writer,
    //    "# Random BED file for testing/benchmarking purposes\n\
    //     # Generated records: {}",
    //    args.num_records
    //)?;

    // Generate and write records
    let records = generate_random_bed_records(args.num_records, args.seed);

    // Re-usable line buffer
    let mut line_buffer = String::new();

    for (chrom, record) in records {
        line_buffer.clear(); // Clear the buffer for the next record
        write!(
            line_buffer,
            "{}\t{}\t{}\t{}",
            chrom, record.start, record.end, record.rest
        )
        .unwrap();
        writeln!(output_writer, "{}", line_buffer)?; // Write the record
    }

    eprintln!("Done!");
    Ok(())
}

fn generate_random_bed_records(num_records: usize, seed: Option<u64>) -> Vec<(String, BedRecord)> {
    const CHROMS: &[&str] = &["chr1", "chr2", "chr3", "chr4", "chr5", "chrX", "chrY"];
    const FEATURE_TYPES: &[&str] = &[
        "gene",
        "exon",
        "promoter",
        "enhancer",
        "UTR",
        "intron",
        "repeat",
        "peak",
        "binding_site",
        "methylation",
    ];

    let mut rng = match seed {
        Some(s) => rand::rngs::StdRng::seed_from_u64(s),
        None => rand::rngs::StdRng::from_entropy(),
    };

    let mut records = Vec::with_capacity(num_records);

    for _ in 0..num_records {
        records.push(generate_single_record(&mut rng, CHROMS, FEATURE_TYPES));
    }

    records.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then(a.1.start.cmp(&b.1.start))
            .then(a.1.end.cmp(&b.1.end))
    });

    records
}

fn generate_single_record<R: Rng>(
    rng: &mut R,
    chroms: &[&str],
    feature_types: &[&str],
) -> (String, BedRecord) {
    let chrom = *chroms.choose(rng).unwrap();
    let start = rng.gen_range(0..1_000_000);
    let length = rng.gen_range(100..10_000);
    let end = start + length;

    let num_extra_fields = rng.gen_range(0..6);
    let rest = (0..num_extra_fields)
        .map(|_| generate_extra_field(rng, feature_types))
        .collect::<Vec<_>>()
        .join("\t");

    (chrom.to_string(), BedRecord { start, end, rest })
}

fn generate_extra_field<R: Rng>(rng: &mut R, feature_types: &[&str]) -> String {
    match rng.gen_range(0..4) {
        0 => feature_types.choose(rng).unwrap().to_string(),
        1 => rng.gen_range(0..1000).to_string(),
        2 => if rng.gen_bool(0.5) { "+" } else { "-" }.to_string(),
        3 => {
            let key = feature_types.choose(rng).unwrap();
            let value = rng.gen_range(0..100);
            format!("{}={}", key, value)
        }
        _ => unreachable!(),
    }
}

#[allow(dead_code)]
fn compare_chromosomes(a: &str, b: &str) -> std::cmp::Ordering {
    fn chrom_value(chrom: &str) -> u32 {
        chrom
            .strip_prefix("chr")
            .and_then(|s| match s {
                "X" => Some(u32::MAX - 1),
                "Y" => Some(u32::MAX),
                num => num.parse().ok(),
            })
            .unwrap_or(0)
    }

    chrom_value(a).cmp(&chrom_value(b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::NamedTempFile;

    #[test]
    fn test_reproducible_generation() {
        let seed = 42;
        let records1 = generate_random_bed_records(100, Some(seed));
        let records2 = generate_random_bed_records(100, Some(seed));
        assert_eq!(records1, records2);
    }

    #[test]
    fn test_chromosome_ordering() {
        assert!(compare_chromosomes("chr1", "chr2").is_lt());
        assert!(compare_chromosomes("chr2", "chr1").is_gt());
        assert!(compare_chromosomes("chrX", "chrY").is_lt());
        assert!(compare_chromosomes("chr1", "chrX").is_lt());
    }

    #[test]
    fn test_output_file_creation() -> Result<(), HgIndexError> {
        let test_file = NamedTempFile::new().unwrap();
        let args = RandomBedArgs {
            output: Some(test_file.path().to_path_buf()),
            num_records: 10,
            seed: Some(42),
        };

        run(args)?;

        let mut content = String::new();
        let mut file = std::fs::File::open(test_file.path())?;
        file.read_to_string(&mut content)?;

        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 10); // 2 comment lines + 10 records

        Ok(())
    }
}
