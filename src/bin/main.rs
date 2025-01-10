#[cfg(feature = "cli")]
mod commands;

#[cfg(feature = "cli")]
mod cli {
    #[cfg(all(feature = "dev"))]
    use super::commands::random_bed::random_bed;
    //#[cfg(all(feature = "dev"))]
    //use crate::commands::analyze;
    use crate::commands::pack::pack;
    use crate::commands::query::query;
    use clap::Parser;
    use hgindex::error::HgIndexError;

    #[derive(Parser)]
    #[command(author, version, about, long_about = None)]
    pub struct Cli {
        #[command(subcommand)]
        command: Commands,
    }

    #[derive(clap::Subcommand)]
    enum Commands {
        //#[cfg(feature = "dev")]
        ///// Analyze index structure and performance metrics
        //Analyze(analyze::AnalyzeArgs),
        /// Block-compress and index a file.
        Pack(pack::PackArgs),
        Query(query::QueryArgs),
        #[cfg(feature = "dev")]
        /// Generate a random BED file for benchmarking (only with dev feature)
        RandomBed(random_bed::RandomBedArgs),
    }

    pub fn run() -> Result<(), HgIndexError> {
        let cli = Cli::parse();
        match cli.command {
            //#[cfg(feature = "dev")]
            //Commands::Analyze(args) => analyze::run(args),
            Commands::Pack(args) => pack::run(args),
            Commands::Query(args) => query::run(args),
            #[cfg(feature = "dev")]
            Commands::RandomBed(args) => random_bed::run(args),
        }
    }
}

fn main() {
    #[cfg(feature = "cli")]
    if let Err(e) = cli::run() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }

    #[cfg(not(feature = "cli"))]
    {
        eprintln!("CLI feature not enabled. Please rebuild with --features cli");
        std::process::exit(1);
    }
}
