use anyhow::{bail, Context, Result};
use std::collections::VecDeque;
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::Duration;

#[derive(Debug, Clone)]
struct Region {
    chrom: String,
    start: u64,
    end: u64,
}

impl Region {
    fn freebayes_arg(&self) -> String {
        format!("{}:{}-{}", self.chrom, self.start, self.end)
    }
}

/// Run freebayes variant calling on one or more BAMs.
///
/// 1. Merge BAMs (if multiple) with samtools merge
/// 2. Read FASTA index for genome regions
/// 3. Run freebayes in parallel across regions
/// 4. Merge per-region VCFs with bcftools
///
/// The final VCF is written to `output_dir/variants.vcf.gz`.
pub fn run_freebayes(
    bams: &[PathBuf],
    reference: &Path,
    threads: u32,
    min_cov: u32,
    output_dir: &Path,
) -> Result<()> {
    fs::create_dir_all(output_dir)?;
    let logs_dir = output_dir.join("logs");
    fs::create_dir_all(&logs_dir)?;

    // Step 1: Merge BAMs if multiple, otherwise use the single BAM directly.
    let working_bam = if bams.len() > 1 {
        let merged = output_dir.join("merged.bam");
        eprintln!("merging {} BAMs...", bams.len());
        let mut cmd = Command::new("samtools");
        cmd.arg("merge")
            .arg("-@")
            .arg(threads.to_string())
            .arg(&merged);
        for bam in bams {
            cmd.arg(bam);
        }
        let status = cmd.status().context("failed to run samtools merge")?;
        if !status.success() {
            bail!("samtools merge failed");
        }
        let status = Command::new("samtools")
            .arg("index")
            .arg("-@")
            .arg(threads.to_string())
            .arg(&merged)
            .status()
            .context("failed to run samtools index")?;
        if !status.success() {
            bail!("samtools index failed on merged BAM");
        }
        merged
    } else {
        bams[0].clone()
    };

    // Step 2: Read FASTA index for regions.
    let fai_path = PathBuf::from(format!("{}.fai", reference.display()));
    if !fai_path.exists() {
        bail!(
            "FASTA index not found: {}. Run `samtools faidx` first.",
            fai_path.display()
        );
    }
    let region_groups = split_fasta_regions(&fai_path, threads as usize)?;
    let mut all_regions: Vec<(usize, Region)> = Vec::new();
    let mut idx = 0;
    for group in &region_groups {
        for region in group {
            all_regions.push((idx, region.clone()));
            idx += 1;
        }
    }
    eprintln!(
        "running freebayes across {} regions with {} threads...",
        all_regions.len(),
        threads
    );

    // Step 3: Run freebayes in parallel across regions.
    let mut pending: VecDeque<(usize, Region)> = all_regions.into_iter().collect();
    let mut running: Vec<(usize, std::process::Child)> = Vec::new();
    let mut completed: Vec<PathBuf> = Vec::new();

    while !pending.is_empty() || !running.is_empty() {
        // Harvest completed processes.
        let mut i = 0;
        while i < running.len() {
            match running[i].1.try_wait()? {
                Some(status) => {
                    let (rid, _) = running.remove(i);
                    if !status.success() {
                        bail!(
                            "freebayes failed for region {} (see {}/freebayes_{}.vcf.err)",
                            rid,
                            output_dir.display(),
                            rid
                        );
                    }
                    completed.push(output_dir.join(format!("freebayes_{}.vcf", rid)));
                }
                None => {
                    i += 1;
                }
            }
        }

        // Fill available slots with new processes.
        while running.len() < threads as usize {
            if let Some((rid, region)) = pending.pop_front() {
                let vcf_path = output_dir.join(format!("freebayes_{}.vcf", rid));
                let err_path = output_dir.join(format!("freebayes_{}.vcf.err", rid));
                let vcf_file = File::create(&vcf_path)
                    .with_context(|| format!("create {}", vcf_path.display()))?;
                let err_file = File::create(&err_path)
                    .with_context(|| format!("create {}", err_path.display()))?;

                let child = Command::new("freebayes")
                    .arg("-f")
                    .arg(reference)
                    .args(["-iXu", "-C", "2", "-q", "20", "-n", "3", "-E", "1", "-m", "30"])
                    .arg("--min-coverage")
                    .arg(min_cov.to_string())
                    .args(["--pooled-continuous", "--skip-coverage", "100000"])
                    .arg("-r")
                    .arg(region.freebayes_arg())
                    .arg(&working_bam)
                    .stdout(Stdio::from(vcf_file))
                    .stderr(Stdio::from(err_file))
                    .spawn()
                    .context("failed to spawn freebayes")?;

                running.push((rid, child));
            } else {
                break;
            }
        }

        if !running.is_empty() {
            std::thread::sleep(Duration::from_millis(500));
        }
    }

    // Step 4: bgzip and index each per-region VCF.
    eprintln!("compressing and indexing {} per-region VCFs...", completed.len());
    for vcf in &completed {
        run_cmd("bgzip", &[vcf.as_os_str()])?;
    }

    let gz_files: Vec<PathBuf> = completed
        .iter()
        .map(|f| PathBuf::from(format!("{}.gz", f.display())))
        .collect();

    for gz in &gz_files {
        run_cmd("bcftools", &["index".as_ref(), gz.as_os_str()])?;
    }

    // Step 5: Merge with bcftools concat, sort, bgzip, and tabix.
    eprintln!("merging VCFs...");
    let merged_vcf = output_dir.join("merged.vcf");
    {
        let out_file = File::create(&merged_vcf)?;
        let mut cmd = Command::new("bcftools");
        cmd.arg("concat").arg("-a");
        for gz in &gz_files {
            cmd.arg(gz);
        }
        let status = cmd
            .stdout(Stdio::from(out_file))
            .status()
            .context("bcftools concat")?;
        if !status.success() {
            bail!("bcftools concat failed");
        }
    }

    let sorted_vcf = output_dir.join("variants.vcf");
    {
        let out_file = File::create(&sorted_vcf)?;
        let status = Command::new("bcftools")
            .arg("sort")
            .arg(&merged_vcf)
            .stdout(Stdio::from(out_file))
            .status()
            .context("bcftools sort")?;
        if !status.success() {
            bail!("bcftools sort failed");
        }
    }

    run_cmd("bgzip", &[sorted_vcf.as_os_str()])?;

    let final_vcf = output_dir.join("variants.vcf.gz");
    run_cmd(
        "tabix",
        &["-p".as_ref(), "vcf".as_ref(), final_vcf.as_os_str()],
    )?;

    // Cleanup temp files.
    for gz in &gz_files {
        let _ = fs::remove_file(gz);
        let _ = fs::remove_file(format!("{}.csi", gz.display()));
    }
    let _ = fs::remove_file(&merged_vcf);
    for vcf in &completed {
        let _ = fs::remove_file(format!("{}.err", vcf.display()));
    }
    if bams.len() > 1 {
        let _ = fs::remove_file(output_dir.join("merged.bam"));
        let _ = fs::remove_file(output_dir.join("merged.bam.bai"));
    }

    eprintln!("freebayes complete: {}", final_vcf.display());
    Ok(())
}

/// Split the reference genome into roughly equal-sized region groups for parallel processing.
/// Mirrors the logic in souporcell_pipeline.py `get_fasta_regions`.
fn split_fasta_regions(fai_path: &Path, threads: usize) -> Result<Vec<Vec<Region>>> {
    let file =
        File::open(fai_path).with_context(|| format!("open {}", fai_path.display()))?;
    let reader = BufReader::new(file);

    let mut chroms: Vec<(String, u64)> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 {
            continue;
        }
        let name = fields[0].to_string();
        let length: u64 = fields[1]
            .parse()
            .with_context(|| format!("invalid length in .fai: {}", fields[1]))?;
        if length >= 250_000 {
            chroms.push((name, length));
        }
    }

    if chroms.is_empty() {
        bail!("no chromosomes >= 250 kb found in {}", fai_path.display());
    }

    let total_length: u64 = chroms.iter().map(|(_, l)| *l).sum();
    let step_length = (total_length as f64 / threads as f64).ceil() as u64;

    let mut regions: Vec<Vec<Region>> = Vec::new();
    let mut current_group: Vec<Region> = Vec::new();
    let mut group_so_far: u64 = 0;

    for (chrom, chrom_length) in &chroms {
        let mut chrom_so_far: u64 = 0;
        loop {
            let remaining_in_chrom = chrom_length - chrom_so_far;
            if group_so_far + remaining_in_chrom <= step_length {
                current_group.push(Region {
                    chrom: chrom.clone(),
                    start: chrom_so_far,
                    end: *chrom_length,
                });
                group_so_far += remaining_in_chrom;
                break;
            } else {
                let end = chrom_so_far + step_length - group_so_far;
                current_group.push(Region {
                    chrom: chrom.clone(),
                    start: chrom_so_far,
                    end,
                });
                regions.push(current_group);
                current_group = Vec::new();
                chrom_so_far = end;
                group_so_far = 0;
            }
        }
    }

    if !current_group.is_empty() {
        if regions.len() == threads {
            regions.last_mut().unwrap().extend(current_group);
        } else {
            regions.push(current_group);
        }
    }

    Ok(regions)
}

/// Run a simple command, bail on failure.
fn run_cmd(program: &str, args: &[&std::ffi::OsStr]) -> Result<()> {
    let status = Command::new(program)
        .args(args)
        .status()
        .with_context(|| format!("failed to run {}", program))?;
    if !status.success() {
        bail!("{} failed", program);
    }
    Ok(())
}
