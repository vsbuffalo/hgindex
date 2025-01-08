// bin/commands/random_bed.rs

#[cfg(all(feature = "cli", feature = "dev"))]
pub mod random_bed {
    use crate::commands::BedRecord;
    use clap::Args;
    use hgindex::error::HgIndexError;
    use hgindex::io::io::OutputFile;
    use rand::{seq::SliceRandom, Rng};
    use std::io::Write;
    use std::path::PathBuf;

    #[derive(Args)]
    pub struct RandomBedArgs {
        /// Output file path (.bed or .bed.gz)
        #[arg(short, long)]
        pub output: PathBuf,

        /// Number of records to generate
        #[arg(short = 'n', long, default_value = "1000000")]
        pub num_records: usize,
    }

    pub fn run(args: RandomBedArgs) -> Result<(), HgIndexError> {
        eprintln!(
            "Generating {} random BED records to {}",
            args.num_records,
            args.output.display()
        );

        let output = OutputFile::new(&args.output);
        let mut writer = output.writer()?;

        //// Write header comments
        //writeln!(
        //    writer,
        //    "# Random BED file for testing/benchmarking purposes\n\
        //     # Generated records: {}",
        //    args.num_records
        //)?;

        // Generate and write records
        let records = generate_random_bed_records(args.num_records);
        for record in records {
            writeln!(writer, "{}", record.to_string())?;
        }

        eprintln!("Done!");
        Ok(())
    }

    fn generate_random_bed_records(num_records: usize) -> Vec<BedRecord> {
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

        let mut rng = rand::thread_rng();
        let mut records = Vec::with_capacity(num_records);

        for _ in 0..num_records {
            records.push(generate_single_record(&mut rng, CHROMS, FEATURE_TYPES));
        }

        records.sort_by(|a, b| {
            let chrom_cmp = compare_chromosomes(&a.chrom, &b.chrom);
            chrom_cmp.then(a.start.cmp(&b.start))
        });

        records
    }

    fn generate_single_record<R: Rng>(
        rng: &mut R,
        chroms: &[&str],
        feature_types: &[&str],
    ) -> BedRecord {
        let chrom = *chroms.choose(rng).unwrap();
        let start = rng.gen_range(0..1_000_000);
        let length = rng.gen_range(100..10_000);
        let end = start + length;

        let num_extra_fields = rng.gen_range(0..6);
        let rest = (0..num_extra_fields)
            .map(|_| generate_extra_field(rng, feature_types))
            .collect::<Vec<_>>()
            .join("\t");

        BedRecord {
            chrom: chrom.to_string(),
            start,
            end,
            rest,
        }
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
        fn test_bed_record_generation() {
            let mut rng = rand::thread_rng();
            let record = generate_single_record(&mut rng, &["chr1", "chr2"], &["gene", "exon"]);

            assert!(record.chrom.starts_with("chr"));
            assert!(record.end > record.start);
        }

        #[test]
        fn test_chromosome_ordering() {
            assert!(compare_chromosomes("chr1", "chr2").is_lt());
            assert!(compare_chromosomes("chr2", "chr1").is_gt());
            assert!(compare_chromosomes("chrX", "chrY").is_lt());
            assert!(compare_chromosomes("chr1", "chrX").is_lt());
        }

        #[test]
        fn test_output_file_creation() -> Result<(), RandomBedError> {
            let test_file = NamedTempFile::new().unwrap();
            let args = RandomBedArgs {
                output: test_file.path().to_path_buf(),
                num_records: 10,
            };

            run(args)?;

            let mut content = String::new();
            let mut file = std::fs::File::open(test_file.path())?;
            file.read_to_string(&mut content)?;

            let lines: Vec<&str> = content.lines().collect();
            assert!(lines.len() >= 12); // 2 comment lines + 10 records
            assert!(lines[0].starts_with('#'));

            Ok(())
        }
    }
}
