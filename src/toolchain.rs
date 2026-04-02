use crate::paths;
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ToolRecord {
    name: String,
    binary_path: String,
    source_path: String,
}

const VARTRIX_REPO: &str = "https://github.com/10XGenomics/vartrix.git";
const SOUPORCELL_REPO: &str = "https://github.com/wheaton5/souporcell.git";

fn vendor_dir() -> Result<PathBuf> {
    let repo_root = std::env::current_dir().context("failed to read current directory")?;
    Ok(repo_root.join("vendor"))
}

fn git_clone(url: &str, dest: &Path) -> Result<()> {
    println!("Cloning {} -> {}", url, dest.display());
    let status = Command::new("git")
        .args(["clone", url])
        .arg(dest)
        .status()
        .with_context(|| format!("failed to run git clone {}", url))?;
    if !status.success() {
        bail!("git clone failed for {}", url);
    }
    Ok(())
}

fn git_pull(repo: &Path) -> Result<()> {
    println!("Pulling latest in {}", repo.display());
    let status = Command::new("git")
        .args(["pull"])
        .current_dir(repo)
        .status()
        .with_context(|| format!("failed to run git pull in {}", repo.display()))?;
    if !status.success() {
        bail!("git pull failed in {}", repo.display());
    }
    Ok(())
}

pub fn fetch() -> Result<()> {
    let vendor = vendor_dir()?;
    fs::create_dir_all(&vendor)
        .with_context(|| format!("failed to create {}", vendor.display()))?;

    let vartrix_src = vendor.join("vartrix");
    let soup_src = vendor.join("souporcell");

    if vartrix_src.exists() {
        println!("vendor/vartrix already exists, skipping (use `tools update` to pull latest)");
    } else {
        git_clone(VARTRIX_REPO, &vartrix_src)?;
    }

    if soup_src.exists() {
        println!("vendor/souporcell already exists, skipping (use `tools update` to pull latest)");
    } else {
        git_clone(SOUPORCELL_REPO, &soup_src)?;
    }

    println!("\nFetch complete. Run `souporcellx tools bootstrap` to build.");
    Ok(())
}

pub fn update() -> Result<()> {
    let vendor = vendor_dir()?;
    let vartrix_src = vendor.join("vartrix");
    let soup_src = vendor.join("souporcell");

    if !vartrix_src.exists() || !soup_src.exists() {
        bail!("vendor sources not found. Run `souporcellx tools fetch` first.");
    }

    git_pull(&vartrix_src)?;
    git_pull(&soup_src)?;

    println!("\nSources updated. Rebuilding...\n");
    bootstrap()
}

pub fn bootstrap() -> Result<()> {
    let vendor = vendor_dir()?;
    let vartrix_src = vendor.join("vartrix");
    let soup_src = vendor.join("souporcell");
    let souporcell_rust = soup_src.join("souporcell");
    let troublet_rust = soup_src.join("troublet");

    for p in [&vartrix_src, &souporcell_rust, &troublet_rust] {
        if !p.exists() {
            bail!(
                "missing vendored source directory: {}. Run `souporcellx tools fetch` first.",
                p.display()
            );
        }
    }

    let root = paths::tool_root()?;
    let vartrix_bin = build_cargo_bin(&vartrix_src, "vartrix", &root)?;
    let souporcell_bin = build_cargo_bin(&souporcell_rust, "souporcell", &root)?;
    let troublet_bin = build_cargo_bin(&troublet_rust, "troublet", &root)?;

    let records = vec![
        ToolRecord {
            name: "vartrix".into(),
            binary_path: vartrix_bin.display().to_string(),
            source_path: vartrix_src.display().to_string(),
        },
        ToolRecord {
            name: "souporcell".into(),
            binary_path: souporcell_bin.display().to_string(),
            source_path: souporcell_rust.display().to_string(),
        },
        ToolRecord {
            name: "troublet".into(),
            binary_path: troublet_bin.display().to_string(),
            source_path: troublet_rust.display().to_string(),
        },
    ];
    write_registry(&records)?;

    println!("Built managed tools:");
    show()?;
    Ok(())
}

pub fn show() -> Result<()> {
    let registry = read_registry()?;
    if registry.is_empty() {
        println!("No managed tools registered. Run: souporcellx tools bootstrap");
        return Ok(());
    }
    for rec in registry {
        println!("{}  {}  {}", rec.name, rec.binary_path, rec.source_path);
    }
    Ok(())
}

pub fn resolve_binary(name: &str) -> Result<PathBuf> {
    let registry = read_registry()?;
    registry
        .into_iter()
        .find(|r| r.name == name)
        .map(|r| PathBuf::from(r.binary_path))
        .ok_or_else(|| anyhow::anyhow!("managed tool not found: {}. Run tools bootstrap.", name))
}

fn build_cargo_bin(src: &Path, bin_name: &str, root: &Path) -> Result<PathBuf> {
    println!("Building {} in {}", bin_name, src.display());
    let status = Command::new("cargo")
        .arg("build")
        .arg("--release")
        .current_dir(src)
        .status()
        .with_context(|| format!("failed to run cargo build in {}", src.display()))?;
    if !status.success() {
        bail!("cargo build failed in {}", src.display());
    }

    let built = src.join("target").join("release").join(bin_name);
    if !built.exists() {
        bail!("built binary not found: {}", built.display());
    }

    let dest_dir = root.join(bin_name).join("bin");
    fs::create_dir_all(&dest_dir)
        .with_context(|| format!("failed to create {}", dest_dir.display()))?;
    let dest = dest_dir.join(bin_name);
    fs::copy(&built, &dest)
        .with_context(|| format!("failed to copy {} to {}", built.display(), dest.display()))?;
    Ok(dest)
}

fn write_registry(records: &[ToolRecord]) -> Result<()> {
    let path = paths::registry_file()?;
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b',')
        .from_path(&path)
        .with_context(|| format!("failed to open {} for writing", path.display()))?;
    for rec in records {
        wtr.serialize(rec)?;
    }
    wtr.flush()?;
    Ok(())
}

fn read_registry() -> Result<Vec<ToolRecord>> {
    let path = paths::registry_file()?;
    if !path.exists() {
        return Ok(Vec::new());
    }
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b',')
        .from_path(&path)
        .with_context(|| format!("failed to open registry {}", path.display()))?;
    let mut out = Vec::new();
    for rec in rdr.deserialize() {
        out.push(rec?);
    }
    Ok(out)
}
