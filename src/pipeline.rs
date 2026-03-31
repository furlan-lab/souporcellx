use crate::cli::RunArgs;
use crate::manifest::{read_sample_manifest, read_vcf_manifest, validate_sample_rows, validate_vcf_rows};
use crate::slurm::{format_command, submit, JobSpec};
use crate::toolchain;
use anyhow::{Context, Result};
use std::collections::BTreeMap;
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
    let self_bin = std::env::current_exe().context("cannot determine souporcellx binary path")?;

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

        // Track vartrix jobs and output dirs per group.
        // Each entry: (jobid, output_dir, barcode_prefix)
        let mut groups: BTreeMap<&str, Vec<(String, String, String)>> = BTreeMap::new();

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
            groups
                .entry(row.group_id.as_str())
                .or_default()
                .push((
                    jobid,
                    row_out.display().to_string(),
                    row.barcode_prefix().to_string(),
                ));
            total += 1;
        }

        // For each group: combine, then souporcell + troublet per K.
        for (group_id, entries) in &groups {
            let combined_dir = vcf_base.join("combined").join(group_id);

            let vartrix_jobids: Vec<&str> = entries.iter().map(|(id, _, _)| id.as_str()).collect();
            let input_args: String = entries
                .iter()
                .map(|(_, dir, _)| shell_escape(dir.clone()))
                .collect::<Vec<_>>()
                .join(" ");
            let label_args: String = entries
                .iter()
                .map(|(_, _, prefix)| shell_escape(prefix.clone()))
                .collect::<Vec<_>>()
                .join(" ");

            let combine_cmd = format!(
                "{} combine --inputs {} --labels {} --output {}",
                shell_escape(self_bin.display().to_string()),
                input_args,
                label_args,
                shell_escape(combined_dir.display().to_string()),
            );

            let dep = Some(vartrix_jobids.join(":"));
            let combine_spec = JobSpec {
                job_name: format!("combine_{}_{}", vcf.vcf_id, group_id),
                partition: args.partition.clone(),
                cpus: 4,
                mem_mb: args.light_mem_mb,
                dependency: dep,
                command: combine_cmd,
                log_path: logdir
                    .join(format!("combine_{}_{}.log", vcf.vcf_id, group_id))
                    .display()
                    .to_string(),
            };
            let combine_jobid = dispatch(&combine_spec)?;
            total += 1;

            for k in &ks {
                let outdir = vcf_base.join(format!("souporcell_{}", k)).join(group_id);
                if !dry_run {
                    fs::create_dir_all(&outdir)?;
                }
                let cluster_spec = JobSpec {
                    job_name: format!("soup_{}_{}_k{}", vcf.vcf_id, group_id, k),
                    partition: args.partition.clone(),
                    cpus: args.souporcell_threads,
                    mem_mb: args.heavy_mem_mb,
                    dependency: Some(combine_jobid.clone()),
                    command: format!(
                        "{} -a {} -r {} -b {} -k {} -t {} > {} 2> {}",
                        shell_escape(souporcell.display().to_string()),
                        shell_escape(combined_dir.join("alt.mtx").display().to_string()),
                        shell_escape(combined_dir.join("ref.mtx").display().to_string()),
                        shell_escape(combined_dir.join("barcodes.tsv").display().to_string()),
                        k,
                        args.souporcell_threads,
                        shell_escape(outdir.join("clusters_tmp.tsv").display().to_string()),
                        shell_escape(outdir.join("log.tsv").display().to_string()),
                    ),
                    log_path: logdir
                        .join(format!("soup_{}_{}_k{}.log", vcf.vcf_id, group_id, k))
                        .display()
                        .to_string(),
                };
                let cluster_jobid = dispatch(&cluster_spec)?;
                total += 1;

                let troublet_spec = JobSpec {
                    job_name: format!("troublet_{}_{}_k{}", vcf.vcf_id, group_id, k),
                    partition: args.partition.clone(),
                    cpus: 1,
                    mem_mb: args.light_mem_mb,
                    dependency: Some(cluster_jobid),
                    command: format!(
                        "{} -a {} -r {} --clusters {} > {}",
                        shell_escape(troublet.display().to_string()),
                        shell_escape(combined_dir.join("alt.mtx").display().to_string()),
                        shell_escape(combined_dir.join("ref.mtx").display().to_string()),
                        shell_escape(outdir.join("clusters_tmp.tsv").display().to_string()),
                        shell_escape(outdir.join("clusters.tsv").display().to_string()),
                    ),
                    log_path: logdir
                        .join(format!("troublet_{}_{}_k{}.log", vcf.vcf_id, group_id, k))
                        .display()
                        .to_string(),
                };
                dispatch(&troublet_spec)?;
                total += 1;
            }
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
