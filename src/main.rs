mod cli;
mod combine;
mod filter_vcf;
mod freebayes;
mod manifest;
mod paths;
mod pipeline;
mod renamer;
mod retag;
mod slurm;
mod stage;
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
        Commands::Manifest {
            cellranger_dirs,
            group_id,
            prefixes,
            output,
        } => manifest::generate_manifest_from_cellranger(
            &cellranger_dirs,
            &group_id,
            prefixes.as_deref(),
            output.as_deref(),
        ),
        Commands::FilterVcf {
            vcf,
            bams,
            min_cov,
            output,
        } => filter_vcf::filter_vcf(&vcf, &bams, min_cov, &output),
        Commands::Rename {
            bam,
            barcodes,
            output,
            umi_tag,
            cell_tag,
            no_umi,
        } => renamer::rename_bam_to_fastq(&bam, &barcodes, &output, no_umi, &umi_tag, &cell_tag),
        Commands::Retag {
            bam,
            output,
            umi_tag,
            cell_tag,
            no_umi,
        } => retag::retag_bam(&bam, &output, no_umi, &umi_tag, &cell_tag),
        Commands::Freebayes {
            bams,
            r#ref,
            threads,
            min_cov,
            output_dir,
        } => freebayes::run_freebayes(&bams, &r#ref, threads, min_cov, &output_dir),
        Commands::Combine {
            inputs,
            labels,
            output,
        } => combine::combine_matrices(&inputs, &labels, &output),
        Commands::Validate {
            sample_manifest,
            vcf_manifest,
        } => pipeline::validate_manifests(&sample_manifest, &vcf_manifest),
        Commands::Stage {
            cellranger_dirs,
            dest,
            souporcell_dirs,
            include_vdj,
            dry_run,
        } => stage::run(stage::StageArgs {
            cellranger_dirs,
            dest,
            souporcell_dirs,
            include_vdj,
            dry_run,
        }),
        Commands::Run(args) => pipeline::run(args),
    }
}
