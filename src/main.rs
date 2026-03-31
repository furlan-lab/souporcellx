mod cli;
mod manifest;
mod paths;
mod pipeline;
mod slurm;
mod toolchain;

use anyhow::Result;
use clap::Parser;
use cli::{Cli, Commands, ToolCommands};

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Tools { command } => match command {
            ToolCommands::Fetch => toolchain::fetch(),
            ToolCommands::Update => toolchain::update(),
            ToolCommands::Bootstrap => toolchain::bootstrap(),
            ToolCommands::Show => toolchain::show(),
        },
        Commands::Validate {
            sample_manifest,
            vcf_manifest,
        } => pipeline::validate_manifests(&sample_manifest, &vcf_manifest),
        Commands::Run(args) => pipeline::run(args),
    }
}
