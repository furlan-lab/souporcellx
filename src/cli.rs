use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(name = "souporcellx")]
#[command(about = "Souporcell workflow orchestrator for Slurm clusters")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Manage pinned locally-built Rust tool binaries.
    Tools {
        #[command(subcommand)]
        command: ToolCommands,
    },
    /// Validate sample and VCF manifests.
    Validate {
        #[arg(long)]
        sample_manifest: PathBuf,
        #[arg(long)]
        vcf_manifest: PathBuf,
    },
    /// Submit or print the workflow plan.
    Run(RunArgs),
}

#[derive(Subcommand, Debug)]
pub enum ToolCommands {
    /// Clone upstream vartrix and souporcell repos into vendor/.
    Fetch,
    /// Pull latest upstream sources and rebuild.
    Update,
    /// Build vendored vartrix, souporcell, and troublet from source.
    Bootstrap,
    /// Print discovered managed tool paths.
    Show,
}

#[derive(clap::Args, Debug, Clone)]
pub struct RunArgs {
    #[arg(long)]
    pub sample_manifest: PathBuf,
    #[arg(long)]
    pub vcf_manifest: PathBuf,
    #[arg(long)]
    pub r#ref: PathBuf,
    #[arg(long)]
    pub workdir: PathBuf,
    #[arg(long, default_value = "1,2,3,4,5,6,7,8")]
    pub ks: String,
    #[arg(long, default_value = "campus-new")]
    pub partition: String,
    #[arg(long, default_value_t = 24)]
    pub vartrix_threads: u32,
    #[arg(long, default_value_t = 35)]
    pub souporcell_threads: u32,
    #[arg(long, default_value_t = 16000)]
    pub heavy_mem_mb: u32,
    #[arg(long, default_value_t = 8000)]
    pub light_mem_mb: u32,
    #[arg(long, default_value_t = false)]
    pub submit: bool,
}
