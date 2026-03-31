use anyhow::{bail, Context, Result};
use serde::Deserialize;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Deserialize)]
pub struct SampleRow {
    pub group_id: String,
    pub library_id: String,
    pub bam: PathBuf,
    pub barcodes: PathBuf,
    /// Optional barcode prefix. Falls back to library_id if not set.
    #[serde(default)]
    pub prefix: Option<String>,
}

impl SampleRow {
    pub fn barcode_prefix(&self) -> &str {
        self.prefix.as_deref().unwrap_or(&self.library_id)
    }
}

#[derive(Debug, Clone, Deserialize)]
pub struct VcfRow {
    pub vcf_id: String,
    pub vcf_path: PathBuf,
}

pub fn read_sample_manifest(path: &Path) -> Result<Vec<SampleRow>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b',')
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
        .delimiter(b',')
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
    }

    if group_seen.is_empty() {
        bail!("no group_id values found");
    }
    Ok(())
}

pub fn generate_manifest_from_cellranger(
    dirs: &[PathBuf],
    group_id: &str,
    prefixes: Option<&[String]>,
    output: Option<&Path>,
) -> Result<()> {
    let mut rows: Vec<(String, String, PathBuf, PathBuf)> = Vec::new();

    for dir in dirs {
        let canon = dir
            .canonicalize()
            .with_context(|| format!("cannot resolve path: {}", dir.display()))?;

        // Try cellranger count layout first
        let count_bam = canon.join("outs/possorted_genome_bam.bam");
        let count_barcodes = canon.join("outs/filtered_feature_bc_matrix/barcodes.tsv.gz");

        if count_bam.exists() && count_barcodes.exists() {
            let library_id = canon
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_else(|| "unknown".to_string());
            rows.push((group_id.to_string(), library_id, count_bam, count_barcodes));
            continue;
        }

        // Try cellranger multi layout
        let per_sample = canon.join("outs/per_sample_outs");
        if per_sample.is_dir() {
            let mut found = false;
            let mut entries: Vec<_> = std::fs::read_dir(&per_sample)
                .with_context(|| {
                    format!("cannot read per_sample_outs: {}", per_sample.display())
                })?
                .collect::<std::result::Result<Vec<_>, _>>()?;
            entries.sort_by_key(|e| e.file_name());
            for entry in entries {
                if !entry.path().is_dir() {
                    continue;
                }
                let sample_dir = entry.path().join("count");
                let multi_bam = sample_dir.join("sample_alignments.bam");
                let multi_barcodes = sample_dir
                    .join("sample_filtered_feature_bc_matrix/barcodes.tsv.gz");
                if multi_bam.exists() && multi_barcodes.exists() {
                    let library_id = entry.file_name().to_string_lossy().into_owned();
                    rows.push((
                        group_id.to_string(),
                        library_id,
                        multi_bam,
                        multi_barcodes,
                    ));
                    found = true;
                }
            }
            if found {
                continue;
            }
        }

        bail!(
            "no Cell Ranger output found in {}. Expected either \
             outs/possorted_genome_bam.bam (cellranger count) or \
             outs/per_sample_outs/*/count/sample_alignments.bam (cellranger multi)",
            canon.display()
        );
    }

    if rows.is_empty() {
        bail!("no Cell Ranger directories provided");
    }

    if let Some(pfx) = prefixes {
        if pfx.len() != rows.len() {
            bail!(
                "number of prefixes ({}) does not match number of samples found ({})",
                pfx.len(),
                rows.len()
            );
        }
    }

    let writer: Box<dyn std::io::Write> = match output {
        Some(p) => Box::new(
            std::fs::File::create(p)
                .with_context(|| format!("cannot create output file: {}", p.display()))?,
        ),
        None => Box::new(std::io::stdout()),
    };

    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b',')
        .from_writer(writer);

    let has_prefix = prefixes.is_some();
    if has_prefix {
        wtr.write_record(["group_id", "library_id", "bam", "barcodes", "prefix"])?;
    } else {
        wtr.write_record(["group_id", "library_id", "bam", "barcodes"])?;
    }

    for (i, (group, lib, bam, barcodes)) in rows.iter().enumerate() {
        let mut record = vec![
            group.as_str().to_string(),
            lib.clone(),
            bam.display().to_string(),
            barcodes.display().to_string(),
        ];
        if let Some(pfx) = prefixes {
            record.push(pfx[i].clone());
        }
        wtr.write_record(&record)?;
    }

    wtr.flush()?;

    if let Some(p) = output {
        eprintln!("wrote {} rows to {}", rows.len(), p.display());
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
