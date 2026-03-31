use anyhow::{Context, Result};
use dirs::cache_dir;
use std::fs;
use std::path::PathBuf;

pub fn cache_root() -> Result<PathBuf> {
    let base = cache_dir().context("failed to determine user cache directory")?;
    let root = base.join("souporcellx");
    fs::create_dir_all(&root).with_context(|| format!("failed to create {}", root.display()))?;
    Ok(root)
}

pub fn tool_root() -> Result<PathBuf> {
    let p = cache_root()?.join("tools");
    fs::create_dir_all(&p).with_context(|| format!("failed to create {}", p.display()))?;
    Ok(p)
}

pub fn registry_file() -> Result<PathBuf> {
    Ok(tool_root()?.join("registry.tsv"))
}
