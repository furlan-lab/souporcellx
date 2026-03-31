use anyhow::{bail, Context, Result};
use serde::Deserialize;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum CombineMode {
    Sum,
    Concat,
}

#[derive(Debug, Clone, Deserialize)]
pub struct SampleRow {
    pub group_id: String,
    pub library_id: String,
    pub mode: CombineMode,
    pub bam: PathBuf,
    pub barcodes: PathBuf,
}

#[derive(Debug, Clone, Deserialize)]
pub struct VcfRow {
    pub vcf_id: String,
    pub vcf_path: PathBuf,
}

pub fn read_sample_manifest(path: &Path) -> Result<Vec<SampleRow>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to open sample manifest: {}", path.display()))?;
    let mut rows = Vec::new();
    for rec in rdr.deserialize() {
        let row: SampleRow = rec.with_context(|| format!("failed to parse row in {}", path.display()))?;
        rows.push(row);
    }
    if rows.is_empty() {
        bail!("sample manifest is empty: {}", path.display());
    }
    Ok(rows)
}

pub fn read_vcf_manifest(path: &Path) -> Result<Vec<VcfRow>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to open VCF manifest: {}", path.display()))?;
    let mut rows = Vec::new();
    for rec in rdr.deserialize() {
        let row: VcfRow = rec.with_context(|| format!("failed to parse row in {}", path.display()))?;
        rows.push(row);
    }
    if rows.is_empty() {
        bail!("VCF manifest is empty: {}", path.display());
    }
    Ok(rows)
}

pub fn validate_sample_rows(rows: &[SampleRow]) -> Result<()> {
    let mut group_seen = BTreeSet::new();
    let mut repeated: BTreeMap<&str, Vec<&SampleRow>> = BTreeMap::new();

    for row in rows {
        if row.group_id.trim().is_empty() || row.library_id.trim().is_empty() {
            bail!("group_id and library_id must not be empty");
        }
        if !row.bam.exists() {
            bail!("missing BAM: {}", row.bam.display());
        }
        if !row.barcodes.exists() {
            bail!("missing barcode file: {}", row.barcodes.display());
        }
        group_seen.insert(row.group_id.as_str());
        repeated.entry(row.library_id.as_str()).or_default().push(row);
    }

    for (lib, lib_rows) in repeated {
        if lib_rows.len() > 1 && lib_rows.iter().any(|r| r.mode != CombineMode::Sum) {
            bail!(
                "library_id '{}' appears multiple times but not all rows are mode=sum",
                lib
            );
        }
    }

    if group_seen.is_empty() {
        bail!("no group_id values found");
    }
    Ok(())
}

pub fn validate_vcf_rows(rows: &[VcfRow]) -> Result<()> {
    let mut ids = BTreeSet::new();
    for row in rows {
        if row.vcf_id.trim().is_empty() {
            bail!("vcf_id must not be empty");
        }
        if !ids.insert(row.vcf_id.clone()) {
            bail!("duplicate vcf_id: {}", row.vcf_id);
        }
        if !row.vcf_path.exists() {
            bail!("missing VCF: {}", row.vcf_path.display());
        }
    }
    Ok(())
}
