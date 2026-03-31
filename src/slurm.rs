use anyhow::{bail, Context, Result};
use std::process::Command;

#[derive(Debug, Clone)]
pub struct JobSpec {
    pub job_name: String,
    pub partition: String,
    pub cpus: u32,
    pub mem_mb: u32,
    pub dependency: Option<String>,
    pub command: String,
    pub log_path: String,
}

/// Format a JobSpec as the sbatch command line that would be executed.
pub fn format_command(spec: &JobSpec) -> String {
    let mut parts = vec![
        "sbatch".to_string(),
        "-n 1".to_string(),
        format!("-c {}", spec.cpus),
        format!("-p {}", spec.partition),
        format!("--mem-per-cpu {}MB", spec.mem_mb),
        format!("--job-name {}", spec.job_name),
        format!("--output {}", spec.log_path),
    ];

    if let Some(dep) = &spec.dependency {
        parts.push(format!("--dependency afterok:{}", dep));
    }

    parts.push(format!("--wrap '{}'", spec.command.replace('\'', "'\\''")));
    parts.join(" \\\n  ")
}

pub fn submit(spec: &JobSpec) -> Result<String> {
    let mut cmd = Command::new("sbatch");
    cmd.arg("-n")
        .arg("1")
        .arg("-c")
        .arg(spec.cpus.to_string())
        .arg("-p")
        .arg(&spec.partition)
        .arg("--mem-per-cpu")
        .arg(format!("{}MB", spec.mem_mb))
        .arg("--job-name")
        .arg(&spec.job_name)
        .arg("--output")
        .arg(&spec.log_path);

    if let Some(dep) = &spec.dependency {
        cmd.arg("--dependency").arg(format!("afterok:{}", dep));
    }

    cmd.arg("--wrap").arg(&spec.command);

    let out = cmd.output().context("failed to execute sbatch")?;
    if !out.status.success() {
        bail!("sbatch failed: {}", String::from_utf8_lossy(&out.stderr));
    }
    let stdout = String::from_utf8_lossy(&out.stdout);
    let job_id = stdout
        .split_whitespace()
        .find(|tok| tok.chars().all(|c| c.is_ascii_digit()))
        .ok_or_else(|| anyhow::anyhow!("could not parse sbatch output: {}", stdout.trim()))?;
    Ok(job_id.to_string())
}
