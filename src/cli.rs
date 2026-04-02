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
    /// Generate a sample manifest from Cell Ranger output directories.
    Manifest {
        /// One or more Cell Ranger output directories.
        #[arg(long, required = true, num_args = 1..)]
        cellranger_dirs: Vec<PathBuf>,
        /// Group ID for all samples (defaults to "group1").
        #[arg(long, default_value = "group1")]
        group_id: String,
        /// Barcode prefixes (one per sample). If omitted, library_id is used.
        #[arg(long, num_args = 1..)]
        prefixes: Option<Vec<String>>,
        /// Write output to a file instead of stdout.
        #[arg(long)]
        output: Option<PathBuf>,
    },
    /// Combine vartrix output matrices across samples (used internally by the pipeline).
    Combine {
        /// Vartrix output directories to combine.
        #[arg(long, required = true, num_args = 1..)]
        inputs: Vec<PathBuf>,
        /// Label for each input (used as barcode prefix). Must match number of inputs.
        #[arg(long, required = true, num_args = 1..)]
        labels: Vec<String>,
        /// Output directory for combined matrices.
        #[arg(long)]
        output: PathBuf,
    },
    /// Filter a VCF to retain only variants with sufficient BAM coverage (used internally by the pipeline).
    FilterVcf {
        /// Input VCF file (plain or gzipped).
        #[arg(long)]
        vcf: PathBuf,
        /// BAM files to compute coverage from (must be indexed).
        #[arg(long, required = true, num_args = 1..)]
        bams: Vec<PathBuf>,
        /// Minimum combined read depth to retain a variant.
        #[arg(long, default_value_t = 20)]
        min_cov: u32,
        /// Output filtered VCF path.
        #[arg(long)]
        output: PathBuf,
    },
    /// Run freebayes variant calling on BAMs (used internally by the pipeline).
    Freebayes {
        /// Input BAM files (will be merged if multiple).
        #[arg(long, required = true, num_args = 1..)]
        bams: Vec<PathBuf>,
        /// Reference FASTA (must have .fai index).
        #[arg(long)]
        r#ref: PathBuf,
        /// Number of parallel freebayes processes.
        #[arg(long, default_value_t = 24)]
        threads: u32,
        /// Minimum coverage for freebayes.
        #[arg(long, default_value_t = 20)]
        min_cov: u32,
        /// Output directory (final VCF will be variants.vcf.gz).
        #[arg(long)]
        output_dir: PathBuf,
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
    /// VCF manifest. Required unless --remap is used (freebayes will discover variants).
    #[arg(long)]
    pub vcf_manifest: Option<PathBuf>,
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
    /// Minimum alt read depth for coverage filtering and souporcell clustering.
    #[arg(long, default_value_t = 10)]
    pub min_alt: u32,
    /// Minimum ref read depth for coverage filtering and souporcell clustering.
    #[arg(long, default_value_t = 10)]
    pub min_ref: u32,
    /// Skip VCF coverage filtering (pass raw VCFs directly to vartrix).
    #[arg(long, default_value_t = false)]
    pub skip_coverage_filter: bool,
    /// Enable remapping stage (renamer → minimap2 → retag → sort/index) before vartrix.
    #[arg(long, default_value_t = false)]
    pub remap: bool,
    /// Number of threads for minimap2 and samtools during remapping.
    #[arg(long, default_value_t = 24)]
    pub remap_threads: u32,
    /// BAM tag for UMIs (passed to renamer.py / retag.py).
    #[arg(long, default_value = "UB")]
    pub umi_tag: String,
    /// BAM tag for cell barcodes (passed to renamer.py / retag.py).
    #[arg(long, default_value = "CB")]
    pub cell_tag: String,
    /// Set if BAM files lack UMI tags.
    #[arg(long, default_value_t = false)]
    pub no_umi: bool,
}
