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
    // Validate manifests.
    let sample_rows = read_sample_manifest(&args.sample_manifest)?;
    validate_sample_rows(&sample_rows)?;

    let use_freebayes = args.remap && args.vcf_manifest.is_none();

    let vcf_rows = if let Some(ref vcf_path) = args.vcf_manifest {
        let rows = read_vcf_manifest(vcf_path)?;
        validate_vcf_rows(&rows)?;
        rows
    } else if !args.remap {
        anyhow::bail!("--vcf-manifest is required unless --remap is used (freebayes will discover variants)");
    } else {
        Vec::new()
    };

    // Check external tool availability.
    let _ = which("minimap2").context("minimap2 not found on PATH")?;
    if !use_freebayes {
        // freebayes only needed when explicitly using VCF manifest (legacy check).
        // When doing de novo calling, the freebayes subcommand will check for it.
    }
    if use_freebayes {
        let _ = which("freebayes").context("freebayes not found on PATH (required for de novo variant calling)")?;
        let _ = which("bcftools").context("bcftools not found on PATH (required for de novo variant calling)")?;
        let _ = which("bgzip").context("bgzip not found on PATH (required for de novo variant calling)")?;
        let _ = which("tabix").context("tabix not found on PATH (required for de novo variant calling)")?;
    }
    if args.remap {
        let _ = which("samtools").context("samtools not found on PATH (required for --remap)")?;
    }

    let vartrix = toolchain::resolve_binary("vartrix")?;
    let souporcell = toolchain::resolve_binary("souporcell")?;
    let troublet = toolchain::resolve_binary("troublet")?;
    let self_bin = std::env::current_exe().context("cannot determine souporcellx binary path")?;

    let ks = parse_ks(&args.ks)?;

    fs::create_dir_all(&args.workdir)
        .with_context(|| format!("failed to create {}", args.workdir.display()))?;

    println!("Plan:");
    println!("  vartrix  = {}", vartrix.display());
    println!("  souporcell = {}", souporcell.display());
    println!("  troublet = {}", troublet.display());
    if args.remap {
        println!("  remap    = enabled (threads={})", args.remap_threads);
    }
    if use_freebayes {
        println!("  variants = de novo (freebayes)");
    } else {
        println!("  vcfs     = {}", vcf_rows.len());
    }
    println!("  groups   = {}", unique_count(sample_rows.iter().map(|r| &r.group_id)));
    println!("  libraries= {}", unique_count(sample_rows.iter().map(|r| &r.library_id)));
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

    let use_coverage_filter = !args.skip_coverage_filter && !use_freebayes;
    let min_cov = args.min_alt + args.min_ref;

    // Pre-collect group -> sample rows.
    let mut group_samples: BTreeMap<&str, Vec<&crate::manifest::SampleRow>> = BTreeMap::new();
    for row in &sample_rows {
        group_samples.entry(row.group_id.as_str()).or_default().push(row);
    }

    // --- Remap stage (top-level, once per sample, shared across VCFs) ---
    let mut remap_jobs: BTreeMap<String, (String, PathBuf)> = BTreeMap::new();

    if args.remap {
        let remap_base = args.workdir.join("remap");
        let remap_logdir = remap_base.join("logs");
        if !dry_run {
            fs::create_dir_all(&remap_logdir)?;
        }

        let no_umi_flag = if args.no_umi { " --no-umi" } else { "" };

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
                "{self_bin} rename --bam {bam} --barcodes {barcodes} --output {fqgz} \
                 --umi-tag {umi_tag} --cell-tag {cell_tag}{no_umi_flag} && \
                 minimap2 -ax splice -t {threads} -G50k -k 21 \
                 -w 11 --sr -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 \
                 -n2 -m20 -s40 -g2000 -2K50m --secondary=no \
                 {ref_fa} {fqgz} 2> {mm2log} | samtools view -bS -o {remapped_bam} && \
                 {self_bin} retag --bam {remapped_bam} --output {retagged} \
                 --umi-tag {umi_tag} --cell-tag {cell_tag}{no_umi_flag} && \
                 samtools sort -@ {threads} {retagged} -o {sorted} && \
                 samtools index -@ {threads} {sorted}",
                self_bin = shell_escape(self_bin.display().to_string()),
                bam = shell_escape(row.bam.display().to_string()),
                barcodes = shell_escape(row.barcodes.display().to_string()),
                fqgz = shell_escape(remap_dir.join("renamed.fq.gz").display().to_string()),
                umi_tag = &args.umi_tag,
                cell_tag = &args.cell_tag,
                no_umi_flag = no_umi_flag,
                threads = args.remap_threads,
                ref_fa = shell_escape(args.r#ref.display().to_string()),
                remapped_bam = shell_escape(remap_dir.join("remapped.bam").display().to_string()),
                mm2log = shell_escape(remap_dir.join("minimap2.log").display().to_string()),
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

    // --- Build VCF panels ---
    // Each panel has a vcf_id and maps group_id -> (vcf_path, upstream_dep).
    // For manifest entries: one panel per VCF row.
    // For freebayes: one panel ("denovo") with per-group discovered VCFs.
    struct VcfPanel {
        vcf_id: String,
        /// group_id -> (vcf_path_string, dependency_jobid)
        group_vcfs: BTreeMap<String, (String, Option<String>)>,
    }

    let mut panels: Vec<VcfPanel> = Vec::new();

    if use_freebayes {
        // --- Freebayes: submit one job per group, producing a per-group VCF ---
        let fb_base = args.workdir.join("denovo");
        let fb_logdir = fb_base.join("logs");
        if !dry_run {
            fs::create_dir_all(&fb_logdir)?;
        }

        let mut group_vcfs: BTreeMap<String, (String, Option<String>)> = BTreeMap::new();

        for (group_id, rows) in &group_samples {
            let fb_dir = fb_base.join("variants").join(group_id);
            if !dry_run {
                fs::create_dir_all(&fb_dir)?;
            }
            let fb_vcf = fb_dir.join("variants.vcf.gz");

            let bam_args: String = rows
                .iter()
                .map(|r| {
                    let path = resolve_bam_path(r, &remap_jobs);
                    format!("--bams {}", shell_escape(path))
                })
                .collect::<Vec<_>>()
                .join(" ");

            let fb_cmd = format!(
                "{} freebayes {} --ref {} --threads {} --min-cov {} --output-dir {}",
                shell_escape(self_bin.display().to_string()),
                bam_args,
                shell_escape(args.r#ref.display().to_string()),
                args.remap_threads,
                min_cov,
                shell_escape(fb_dir.display().to_string()),
            );

            // Freebayes depends on all remap jobs for this group.
            let dep_ids: Vec<&str> = rows
                .iter()
                .filter_map(|r| {
                    let key = remap_key_for_row(r);
                    remap_jobs.get(&key).map(|(id, _)| id.as_str())
                })
                .collect();
            let fb_dep = if dep_ids.is_empty() { None } else { Some(dep_ids.join(":")) };

            let fb_spec = JobSpec {
                job_name: format!("freebayes_{}", group_id),
                partition: args.partition.clone(),
                cpus: args.remap_threads,
                mem_mb: args.heavy_mem_mb,
                dependency: fb_dep,
                command: fb_cmd,
                log_path: fb_logdir
                    .join(format!("freebayes_{}.log", group_id))
                    .display()
                    .to_string(),
            };
            let fb_jobid = dispatch(&fb_spec)?;
            group_vcfs.insert(
                group_id.to_string(),
                (fb_vcf.display().to_string(), Some(fb_jobid)),
            );
            total += 1;
        }

        panels.push(VcfPanel {
            vcf_id: "denovo".to_string(),
            group_vcfs,
        });
    } else {
        // --- VCF manifest: one panel per VCF row ---
        for vcf in &vcf_rows {
            let vcf_base = args.workdir.join(&vcf.vcf_id);
            let logdir = vcf_base.join("logs");
            if !dry_run {
                fs::create_dir_all(&logdir)?;
            }

            let mut group_vcfs: BTreeMap<String, (String, Option<String>)> = BTreeMap::new();

            if use_coverage_filter {
                for (group_id, rows) in &group_samples {
                    let filter_dir = vcf_base.join("filtered").join(group_id);
                    if !dry_run {
                        fs::create_dir_all(&filter_dir)?;
                    }
                    let filtered_vcf = filter_dir.join("filtered.vcf.gz");

                    let bam_args: String = rows
                        .iter()
                        .map(|r| {
                            let bam_path = resolve_bam_path(r, &remap_jobs);
                            shell_escape(bam_path)
                        })
                        .collect::<Vec<_>>()
                        .join(" ");

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
                    group_vcfs.insert(
                        group_id.to_string(),
                        (filtered_vcf.display().to_string(), Some(filter_jobid)),
                    );
                    total += 1;
                }
            } else {
                // No coverage filter: all groups use the raw VCF.
                for group_id in group_samples.keys() {
                    let remap_dep = if args.remap {
                        let rows = &group_samples[group_id];
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
                    group_vcfs.insert(
                        group_id.to_string(),
                        (vcf.vcf_path.display().to_string(), remap_dep),
                    );
                }
            }

            panels.push(VcfPanel {
                vcf_id: vcf.vcf_id.clone(),
                group_vcfs,
            });
        }
    }

    // --- Process each VCF panel: vartrix → combine → souporcell → troublet ---
    for panel in &panels {
        let vcf_base = args.workdir.join(&panel.vcf_id);
        let logdir = vcf_base.join("logs");
        if !dry_run {
            fs::create_dir_all(&logdir)?;
        }

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

            let (vcf_path, upstream_dep) = panel
                .group_vcfs
                .get(row.group_id.as_str())
                .expect("group missing from panel");

            let bam_for_vartrix = resolve_bam_path(row, &remap_jobs);

            let cmd = format!(
                "{} --bam {} --cell-barcodes {} --vcf {} --fasta {} --scoring-method coverage --out-matrix {} --ref-matrix {} --out-barcodes {}",
                shell_escape(vartrix.display().to_string()),
                shell_escape(bam_for_vartrix),
                shell_escape(row.barcodes.display().to_string()),
                shell_escape(vcf_path.clone()),
                shell_escape(args.r#ref.display().to_string()),
                shell_escape(row_out.join("alt.mtx").display().to_string()),
                shell_escape(row_out.join("ref.mtx").display().to_string()),
                shell_escape(row_out.join("barcodes.tsv").display().to_string()),
            );

            let spec = JobSpec {
                job_name: format!("vtx_{}_{}_{}", panel.vcf_id, row.library_id, bam_base),
                partition: args.partition.clone(),
                cpus: args.vartrix_threads,
                mem_mb: args.heavy_mem_mb,
                dependency: upstream_dep.clone(),
                command: cmd,
                log_path: logdir
                    .join(format!("vtx_{}_{}_{}.log", panel.vcf_id, row.library_id, bam_base))
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
                job_name: format!("combine_{}_{}", panel.vcf_id, group_id),
                partition: args.partition.clone(),
                cpus: 4,
                mem_mb: args.light_mem_mb,
                dependency: dep,
                command: combine_cmd,
                log_path: logdir
                    .join(format!("combine_{}_{}.log", panel.vcf_id, group_id))
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
                    job_name: format!("soup_{}_{}_k{}", panel.vcf_id, group_id, k),
                    partition: args.partition.clone(),
                    cpus: args.souporcell_threads,
                    mem_mb: args.heavy_mem_mb,
                    dependency: Some(combine_jobid.clone()),
                    command: format!(
                        "{} -a {} -r {} -b {} -k {} -t {} --min_alt {} --min_ref {}{}{} > {} 2> {}",
                        shell_escape(souporcell.display().to_string()),
                        shell_escape(combined_dir.join("alt.mtx").display().to_string()),
                        shell_escape(combined_dir.join("ref.mtx").display().to_string()),
                        shell_escape(combined_dir.join("barcodes.tsv").display().to_string()),
                        k,
                        args.souporcell_threads,
                        args.min_alt,
                        args.min_ref,
                        if args.souporcell3 { " -s true" } else { "" },
                        match args.clustering_method.as_deref() {
                            Some(m) => format!(" -m {}", m),
                            None => String::new(),
                        },
                        shell_escape(outdir.join("clusters_tmp.tsv").display().to_string()),
                        shell_escape(outdir.join("log.tsv").display().to_string()),
                    ),
                    log_path: logdir
                        .join(format!("soup_{}_{}_k{}.log", panel.vcf_id, group_id, k))
                        .display()
                        .to_string(),
                };
                let cluster_jobid = dispatch(&cluster_spec)?;
                total += 1;

                let troublet_spec = JobSpec {
                    job_name: format!("troublet_{}_{}_k{}", panel.vcf_id, group_id, k),
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
                        .join(format!("troublet_{}_{}_k{}.log", panel.vcf_id, group_id, k))
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
