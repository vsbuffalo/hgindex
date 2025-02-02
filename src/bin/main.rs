#[cfg(feature = "cli")]
mod commands;

#[cfg(all(feature = "cli", feature = "dev"))]
use crate::commands::random_bed;
//#[cfg(all(feature = "dev"))]
//use crate::commands::analyze;
use crate::commands::pack;
use crate::commands::query;
use crate::commands::stats;
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
    #[cfg(all(feature = "cli", feature = "dev"))]
    /// Generate a random BED file for benchmarking (only with dev feature)
    RandomBed(random_bed::RandomBedArgs),
    Stats(stats::StatsArgs),
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
        Commands::Stats(args) => stats::run(args),
    }
}

fn main() {
    #[cfg(feature = "cli")]
    if let Err(e) = run() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }

    #[cfg(not(feature = "cli"))]
    {
        eprintln!("CLI feature not enabled. Please rebuild with --features cli");
        std::process::exit(1);
    }
}
