use crate::cli::RunArgs;
use crate::manifest::{read_sample_manifest, read_vcf_manifest, validate_sample_rows, validate_vcf_rows};
use crate::slurm::{format_command, submit, JobSpec};
use crate::toolchain;
use anyhow::{Context, Result};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
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

    if args.remap {
        let _ = which("samtools").context("samtools not found on PATH (required for --remap)")?;
        let _ = resolve_python().context("python3/python not found on PATH (required for --remap)")?;
    }

    let vartrix = toolchain::resolve_binary("vartrix")?;
    let souporcell = toolchain::resolve_binary("souporcell")?;
    let troublet = toolchain::resolve_binary("troublet")?;
    let self_bin = std::env::current_exe().context("cannot determine souporcellx binary path")?;

    let (renamer, retagger) = if args.remap {
        (
            Some(toolchain::resolve_script("renamer.py")?),
            Some(toolchain::resolve_script("retag.py")?),
        )
    } else {
        (None, None)
    };

    let sample_rows = read_sample_manifest(&args.sample_manifest)?;
    let vcf_rows = read_vcf_manifest(&args.vcf_manifest)?;
    let ks = parse_ks(&args.ks)?;

    fs::create_dir_all(&args.workdir)
        .with_context(|| format!("failed to create {}", args.workdir.display()))?;

    println!("Plan:");
    println!("  vartrix  = {}", vartrix.display());
    println!("  souporcell = {}", souporcell.display());
    println!("  troublet = {}", troublet.display());
    if args.remap {
        println!("  renamer  = {}", renamer.as_ref().unwrap().display());
        println!("  retag    = {}", retagger.as_ref().unwrap().display());
        println!("  remap    = enabled (threads={})", args.remap_threads);
    }
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

    let use_coverage_filter = !args.skip_coverage_filter;
    let min_cov = args.min_alt + args.min_ref;

    // Pre-collect group -> sample rows for filter job generation.
    let mut group_samples: BTreeMap<&str, Vec<&crate::manifest::SampleRow>> = BTreeMap::new();
    for row in &sample_rows {
        group_samples.entry(row.group_id.as_str()).or_default().push(row);
    }

    let python = resolve_python().unwrap_or_default();
    let no_umi_str = if args.no_umi { "True" } else { "False" };

    // --- Remap stage (top-level, once per sample, shared across VCFs) ---
    // Maps (library_id, bam_stem) -> (remap_jobid, remapped_bam_path)
    let mut remap_jobs: BTreeMap<String, (String, PathBuf)> = BTreeMap::new();

    if args.remap {
        let remap_base = args.workdir.join("remap");
        let remap_logdir = remap_base.join("logs");
        if !dry_run {
            fs::create_dir_all(&remap_logdir)?;
        }

        let renamer = renamer.as_ref().unwrap();
        let retagger = retagger.as_ref().unwrap();

        for row in &sample_rows {
            let bam_base = row
                .bam
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("bam");
            let remap_key = format!("{}__{}", row.library_id, bam_base);
            let remap_dir = remap_base.join(&remap_key);
            if !dry_run {
                fs::create_dir_all(&remap_dir)?;
            }

            let remapped_bam = remap_dir.join("retagged_sorted.bam");

            let remap_cmd = format!(
                "{python} {renamer} --bam {bam} --barcodes {barcodes} --out {fq} \
                 --no_umi {no_umi} --umi_tag {umi_tag} --cell_tag {cell_tag} && \
                 gzip {fq} && \
                 minimap2 -ax splice -t {threads} -G50k -k 21 \
                 -w 11 --sr -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 \
                 -n2 -m20 -s40 -g2000 -2K50m --secondary=no \
                 {ref_fa} {fqgz} 2> {mm2log} | samtools view -bS -o {remapped_bam} && \
                 {python} {retagger} --sam {remapped_bam} --out {retagged} \
                 --no_umi {no_umi} --umi_tag {umi_tag} --cell_tag {cell_tag} && \
                 samtools sort -@ {threads} {retagged} -o {sorted} && \
                 samtools index -@ {threads} {sorted}",
                python = shell_escape(python.clone()),
                renamer = shell_escape(renamer.display().to_string()),
                bam = shell_escape(row.bam.display().to_string()),
                barcodes = shell_escape(row.barcodes.display().to_string()),
                fq = shell_escape(remap_dir.join("renamed.fq").display().to_string()),
                fqgz = shell_escape(remap_dir.join("renamed.fq.gz").display().to_string()),
                no_umi = no_umi_str,
                umi_tag = &args.umi_tag,
                cell_tag = &args.cell_tag,
                threads = args.remap_threads,
                ref_fa = shell_escape(args.r#ref.display().to_string()),
                remapped_bam = shell_escape(remap_dir.join("remapped.bam").display().to_string()),
                mm2log = shell_escape(remap_dir.join("minimap2.log").display().to_string()),
                retagger = shell_escape(retagger.display().to_string()),
                retagged = shell_escape(remap_dir.join("retagged.bam").display().to_string()),
                sorted = shell_escape(remapped_bam.display().to_string()),
            );

            let remap_spec = JobSpec {
                job_name: format!("remap_{}_{}", row.library_id, bam_base),
                partition: args.partition.clone(),
                cpus: args.remap_threads,
                mem_mb: args.heavy_mem_mb,
                dependency: None,
                command: remap_cmd,
                log_path: remap_logdir
                    .join(format!("remap_{}_{}.log", row.library_id, bam_base))
                    .display()
                    .to_string(),
            };
            let remap_jobid = dispatch(&remap_spec)?;
            remap_jobs.insert(remap_key, (remap_jobid, remapped_bam));
            total += 1;
        }
    }

    for vcf in &vcf_rows {
        let vcf_base = args.workdir.join(&vcf.vcf_id);
        let logdir = vcf_base.join("logs");
        if !dry_run {
            fs::create_dir_all(&logdir)?;
        }

        // Submit coverage filter jobs per group (before vartrix).
        // Maps group_id -> (filter_jobid, filtered_vcf_path)
        let mut filter_jobs: BTreeMap<&str, (String, String)> = BTreeMap::new();
        if use_coverage_filter {
            for (group_id, rows) in &group_samples {
                let filter_dir = vcf_base.join("filtered").join(group_id);
                if !dry_run {
                    fs::create_dir_all(&filter_dir)?;
                }
                let filtered_vcf = filter_dir.join("filtered.vcf.gz");

                // Use remapped BAMs for coverage filtering when --remap is enabled.
                let bam_args: String = rows
                    .iter()
                    .map(|r| {
                        let bam_path = resolve_bam_path(r, &remap_jobs);
                        shell_escape(bam_path)
                    })
                    .collect::<Vec<_>>()
                    .join(" ");

                // Covfilt depends on remap jobs for all samples in this group.
                let covfilt_dep = if args.remap {
                    let dep_ids: Vec<&str> = rows
                        .iter()
                        .filter_map(|r| {
                            let key = remap_key_for_row(r);
                            remap_jobs.get(&key).map(|(id, _)| id.as_str())
                        })
                        .collect();
                    if dep_ids.is_empty() { None } else { Some(dep_ids.join(":")) }
                } else {
                    None
                };

                let filter_cmd = format!(
                    "{} filter-vcf --vcf {} --bams {} --min-cov {} --output {}",
                    shell_escape(self_bin.display().to_string()),
                    shell_escape(vcf.vcf_path.display().to_string()),
                    bam_args,
                    min_cov,
                    shell_escape(filtered_vcf.display().to_string()),
                );

                let filter_spec = JobSpec {
                    job_name: format!("covfilt_{}_{}", vcf.vcf_id, group_id),
                    partition: args.partition.clone(),
                    cpus: 4,
                    mem_mb: args.light_mem_mb,
                    dependency: covfilt_dep,
                    command: filter_cmd,
                    log_path: logdir
                        .join(format!("covfilt_{}_{}.log", vcf.vcf_id, group_id))
                        .display()
                        .to_string(),
                };
                let filter_jobid = dispatch(&filter_spec)?;
                filter_jobs.insert(group_id, (filter_jobid, filtered_vcf.display().to_string()));
                total += 1;
            }
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

            // Use filtered VCF if coverage filtering is enabled, otherwise raw VCF.
            let vcf_path_for_vartrix = if let Some((_, ref filtered_path)) =
                filter_jobs.get(row.group_id.as_str())
            {
                filtered_path.clone()
            } else {
                vcf.vcf_path.display().to_string()
            };

            // Use remapped BAM when --remap is enabled.
            let bam_for_vartrix = resolve_bam_path(row, &remap_jobs);

            let cmd = format!(
                "{} --bam {} --cell-barcodes {} --vcf {} --fasta {} --scoring-method coverage --out-matrix {} --ref-matrix {} --out-barcodes {}",
                shell_escape(vartrix.display().to_string()),
                shell_escape(bam_for_vartrix),
                shell_escape(row.barcodes.display().to_string()),
                shell_escape(vcf_path_for_vartrix),
                shell_escape(args.r#ref.display().to_string()),
                shell_escape(row_out.join("alt.mtx").display().to_string()),
                shell_escape(row_out.join("ref.mtx").display().to_string()),
                shell_escape(row_out.join("barcodes.tsv").display().to_string()),
            );

            // Vartrix depends on covfilt (if present), and also on remap (if present).
            let vartrix_dep = {
                let mut deps = Vec::new();
                if let Some((filter_jobid, _)) = filter_jobs.get(row.group_id.as_str()) {
                    deps.push(filter_jobid.clone());
                } else if let Some((remap_jobid, _)) = remap_jobs.get(&remap_key_for_row(row)) {
                    // Only add remap dep directly if there's no covfilt (covfilt already depends on remap).
                    deps.push(remap_jobid.clone());
                }
                if deps.is_empty() { None } else { Some(deps.join(":")) }
            };

            let spec = JobSpec {
                job_name: format!("vtx_{}_{}_{}", vcf.vcf_id, row.library_id, bam_base),
                partition: args.partition.clone(),
                cpus: args.vartrix_threads,
                mem_mb: args.heavy_mem_mb,
                dependency: vartrix_dep,
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
                        "{} -a {} -r {} -b {} -k {} -t {} --min_alt {} --min_ref {} > {} 2> {}",
                        shell_escape(souporcell.display().to_string()),
                        shell_escape(combined_dir.join("alt.mtx").display().to_string()),
                        shell_escape(combined_dir.join("ref.mtx").display().to_string()),
                        shell_escape(combined_dir.join("barcodes.tsv").display().to_string()),
                        k,
                        args.souporcell_threads,
                        args.min_alt,
                        args.min_ref,
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

/// Find python3 or python on PATH.
fn resolve_python() -> Result<String> {
    if which("python3").is_ok() {
        Ok("python3".to_string())
    } else if which("python").is_ok() {
        Ok("python".to_string())
    } else {
        anyhow::bail!("neither python3 nor python found on PATH")
    }
}

/// Build the remap key for a sample row (matches the remap output directory name).
fn remap_key_for_row(row: &crate::manifest::SampleRow) -> String {
    let bam_base = row
        .bam
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("bam");
    format!("{}__{}", row.library_id, bam_base)
}

/// Return the BAM path to use for a sample: remapped BAM if available, otherwise original.
fn resolve_bam_path(
    row: &crate::manifest::SampleRow,
    remap_jobs: &BTreeMap<String, (String, std::path::PathBuf)>,
) -> String {
    let key = remap_key_for_row(row);
    if let Some((_, ref remapped_bam)) = remap_jobs.get(&key) {
        remapped_bam.display().to_string()
    } else {
        row.bam.display().to_string()
    }
}
