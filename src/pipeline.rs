use crate::cli::RunArgs;
use crate::manifest::{read_sample_manifest, read_vcf_manifest, validate_sample_rows, validate_vcf_rows};
use crate::slurm::{format_command, submit, JobSpec};
use crate::toolchain;
use anyhow::{Context, Result};
use std::fs;
use std::path::Path;
use which::which;

pub fn validate_manifests(sample_manifest: &Path, vcf_manifest: &Path) -> Result<()> {
    let sample_rows = read_sample_manifest(sample_manifest)?;
    validate_sample_rows(&sample_rows)?;

    let vcf_rows = read_vcf_manifest(vcf_manifest)?;
    validate_vcf_rows(&vcf_rows)?;

    println!(
        "Validated {} sample rows and {} VCF rows.",
        sample_rows.len(),
        vcf_rows.len()
    );
    Ok(())
}

pub fn run(args: RunArgs) -> Result<()> {
    validate_manifests(&args.sample_manifest, &args.vcf_manifest)?;
    let _ = which("minimap2").context("minimap2 not found on PATH")?;
    let _ = which("freebayes").context("freebayes not found on PATH")?;

    let vartrix = toolchain::resolve_binary("vartrix")?;
    let souporcell = toolchain::resolve_binary("souporcell")?;
    let troublet = toolchain::resolve_binary("troublet")?;

    let sample_rows = read_sample_manifest(&args.sample_manifest)?;
    let vcf_rows = read_vcf_manifest(&args.vcf_manifest)?;
    let ks = parse_ks(&args.ks)?;

    fs::create_dir_all(&args.workdir)
        .with_context(|| format!("failed to create {}", args.workdir.display()))?;

    println!("Plan:");
    println!("  vartrix  = {}", vartrix.display());
    println!("  souporcell = {}", souporcell.display());
    println!("  troublet = {}", troublet.display());
    println!("  groups   = {}", unique_count(sample_rows.iter().map(|r| &r.group_id)));
    println!("  libraries= {}", unique_count(sample_rows.iter().map(|r| &r.library_id)));
    println!("  vcfs     = {}", vcf_rows.len());
    println!("  ks       = {:?}", ks);

    let dry_run = !args.submit;
    let mut job_counter = 0usize;
    let mut total = 0usize;

    // Helper: either submit the job or print the sbatch command (dry run).
    // Returns a job ID string (real from sbatch, or a placeholder for dry run).
    let mut dispatch = |spec: &JobSpec| -> Result<String> {
        if dry_run {
            println!("{}\n", format_command(spec));
            job_counter += 1;
            Ok(format!("PLACEHOLDER_{}", job_counter))
        } else {
            let jobid = submit(spec)?;
            println!("submitted {} -> {}", spec.job_name, jobid);
            Ok(jobid)
        }
    };

    for vcf in &vcf_rows {
        let vcf_base = args.workdir.join(&vcf.vcf_id);
        let logdir = vcf_base.join("logs");
        if !dry_run {
            fs::create_dir_all(&logdir)?;
        }

        let mut vartrix_jobids = Vec::new();
        for row in &sample_rows {
            let bam_base = row
                .bam
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("bam");
            let row_out = vcf_base
                .join("raw_vartrix")
                .join(format!("{}__{}", row.library_id, bam_base));
            if !dry_run {
                fs::create_dir_all(&row_out)?;
            }

            let cmd = format!(
                "{} --bam {} --cell-barcodes {} --vcf {} --fasta {} --out-matrix {} --ref-matrix {} --out-barcodes {}",
                shell_escape(vartrix.display().to_string()),
                shell_escape(row.bam.display().to_string()),
                shell_escape(row.barcodes.display().to_string()),
                shell_escape(vcf.vcf_path.display().to_string()),
                shell_escape(args.r#ref.display().to_string()),
                shell_escape(row_out.join("alt.mtx").display().to_string()),
                shell_escape(row_out.join("ref.mtx").display().to_string()),
                shell_escape(row_out.join("barcodes.tsv").display().to_string()),
            );
            let spec = JobSpec {
                job_name: format!("vtx_{}_{}_{}", vcf.vcf_id, row.library_id, bam_base),
                partition: args.partition.clone(),
                cpus: args.vartrix_threads,
                mem_mb: args.heavy_mem_mb,
                dependency: None,
                command: cmd,
                log_path: logdir
                    .join(format!("vtx_{}_{}_{}.log", vcf.vcf_id, row.library_id, bam_base))
                    .display()
                    .to_string(),
            };
            let jobid = dispatch(&spec)?;
            vartrix_jobids.push(jobid);
            total += 1;
        }

        let dep = Some(vartrix_jobids.join(":"));
        let combine_spec = JobSpec {
            job_name: format!("combine_{}", vcf.vcf_id),
            partition: args.partition.clone(),
            cpus: 4,
            mem_mb: args.light_mem_mb,
            dependency: dep,
            command: format!(
                "echo TODO: collapse and combine grouped matrices for {}",
                shell_escape(vcf.vcf_id.clone())
            ),
            log_path: logdir
                .join(format!("combine_{}.log", vcf.vcf_id))
                .display()
                .to_string(),
        };
        let combine_jobid = dispatch(&combine_spec)?;
        total += 1;

        for k in &ks {
            let outdir = vcf_base.join(format!("souporcell_{}", k));
            if !dry_run {
                fs::create_dir_all(&outdir)?;
            }
            let cluster_spec = JobSpec {
                job_name: format!("soup_{}_k{}", vcf.vcf_id, k),
                partition: args.partition.clone(),
                cpus: args.souporcell_threads,
                mem_mb: args.heavy_mem_mb,
                dependency: Some(combine_jobid.clone()),
                command: format!(
                    "{} -a {} -r {} -b {} -k {} -t {} > {} 2> {}",
                    shell_escape(souporcell.display().to_string()),
                    shell_escape(vcf_base.join("combined").join("GROUP").join("alt.mtx").display().to_string()),
                    shell_escape(vcf_base.join("combined").join("GROUP").join("ref.mtx").display().to_string()),
                    shell_escape(vcf_base.join("combined").join("GROUP").join("barcodes.tsv").display().to_string()),
                    k,
                    args.souporcell_threads,
                    shell_escape(outdir.join("clusters_tmp.tsv").display().to_string()),
                    shell_escape(outdir.join("log.tsv").display().to_string()),
                ),
                log_path: logdir
                    .join(format!("soup_{}_k{}.log", vcf.vcf_id, k))
                    .display()
                    .to_string(),
            };
            let cluster_jobid = dispatch(&cluster_spec)?;
            total += 1;

            let troublet_spec = JobSpec {
                job_name: format!("troublet_{}_k{}", vcf.vcf_id, k),
                partition: args.partition.clone(),
                cpus: 1,
                mem_mb: args.light_mem_mb,
                dependency: Some(cluster_jobid),
                command: format!(
                    "{} -a {} -r {} --clusters {} > {}",
                    shell_escape(troublet.display().to_string()),
                    shell_escape(vcf_base.join("combined").join("GROUP").join("alt.mtx").display().to_string()),
                    shell_escape(vcf_base.join("combined").join("GROUP").join("ref.mtx").display().to_string()),
                    shell_escape(outdir.join("clusters_tmp.tsv").display().to_string()),
                    shell_escape(outdir.join("clusters.tsv").display().to_string()),
                ),
                log_path: logdir
                    .join(format!("troublet_{}_k{}.log", vcf.vcf_id, k))
                    .display()
                    .to_string(),
            };
            dispatch(&troublet_spec)?;
            total += 1;
        }
    }

    if dry_run {
        println!("Dry run: {} job(s) would be submitted. Re-run with --submit to submit.", total);
    } else {
        println!("Submitted {} jobs total.", total);
    }
    Ok(())
}

fn parse_ks(s: &str) -> Result<Vec<u32>> {
    let mut out = Vec::new();
    for part in s.split(',') {
        let k = part.trim().parse::<u32>().with_context(|| format!("invalid K value: {}", part))?;
        out.push(k);
    }
    Ok(out)
}

fn unique_count<'a>(iter: impl Iterator<Item = &'a String>) -> usize {
    use std::collections::BTreeSet;
    let set: BTreeSet<&String> = iter.collect();
    set.len()
}

fn shell_escape(s: String) -> String {
    format!("'{}'", s.replace('\'', "'\\''"))
}
