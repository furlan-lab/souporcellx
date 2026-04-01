use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use rust_htslib::bam::{self, Read as BamRead};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

const MAX_DEPTH: u32 = 100_000;

/// Count the total read depth at a single genomic position across all BAMs.
///
/// Skips unmapped, secondary, supplementary, duplicate, and QC-failed reads
/// (matching `samtools depth` default behaviour).
fn depth_at_position(
    readers: &mut [bam::IndexedReader],
    chrom: &str,
    pos: i64, // 0-based
) -> Result<u32> {
    let mut total: u32 = 0;
    for reader in readers.iter_mut() {
        let tid = match reader.header().tid(chrom.as_bytes()) {
            Some(t) => t,
            None => continue, // contig not in this BAM — depth is 0
        };
        reader
            .fetch((tid, pos, pos + 1))
            .with_context(|| format!("failed to fetch {}:{}", chrom, pos))?;
        for record in reader.records() {
            let record = record?;
            let flags = record.flags();
            // skip unmapped (0x4), secondary (0x100), duplicate (0x400),
            // qcfail (0x200), supplementary (0x800)
            if flags & 0x704 != 0 || flags & 0x800 != 0 {
                continue;
            }
            total += 1;
        }
    }
    Ok(total)
}

/// Open a VCF (plain or gzipped) and return a boxed buffered reader.
fn open_vcf(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let is_gz = path
        .extension()
        .is_some_and(|ext| ext == "gz" || ext == "bgz");
    if is_gz {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Filter a VCF to only include variants with sufficient read depth in the
/// provided BAM files.
///
/// For each variant position the combined depth across all BAMs is computed
/// using indexed random access.  Variants with depth >= `min_cov` and
/// depth < 100 000 are retained.  VCF header lines are always preserved.
pub fn filter_vcf(
    vcf: &Path,
    bams: &[PathBuf],
    min_cov: u32,
    output: &Path,
) -> Result<()> {
    if bams.is_empty() {
        anyhow::bail!("no BAM files provided");
    }

    // Open indexed BAM readers
    let mut readers: Vec<bam::IndexedReader> = Vec::with_capacity(bams.len());
    for bam_path in bams {
        let reader = bam::IndexedReader::from_path(bam_path)
            .with_context(|| format!("cannot open indexed BAM: {}", bam_path.display()))?;
        readers.push(reader);
    }

    // Ensure output directory exists
    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("cannot create output dir: {}", parent.display()))?;
    }

    let vcf_reader = open_vcf(vcf)?;
    let out_file =
        File::create(output).with_context(|| format!("cannot create {}", output.display()))?;
    let mut writer = BufWriter::new(out_file);

    let mut kept = 0u64;
    let mut skipped = 0u64;

    for line in vcf_reader.lines() {
        let line = line?;

        // Header lines pass through unconditionally
        if line.starts_with('#') {
            writeln!(writer, "{}", line)?;
            continue;
        }

        // Parse CHROM and POS from the data line
        let mut cols = line.splitn(3, '\t');
        let chrom = match cols.next() {
            Some(c) if !c.is_empty() => c,
            _ => continue,
        };
        let pos_str = match cols.next() {
            Some(p) => p,
            None => continue,
        };
        let pos: i64 = pos_str.parse().with_context(|| {
            format!("invalid POS value '{}' in VCF line", pos_str)
        })?;
        // VCF POS is 1-based; rust-htslib uses 0-based coordinates
        let pos0 = pos - 1;

        let depth = depth_at_position(&mut readers, chrom, pos0)?;

        if depth >= min_cov && depth < MAX_DEPTH {
            writeln!(writer, "{}", line)?;
            kept += 1;
        } else {
            skipped += 1;
        }
    }

    writer.flush()?;
    println!(
        "coverage filter: kept {} variants, skipped {} (min_cov={}, max={})",
        kept, skipped, min_cov, MAX_DEPTH
    );
    Ok(())
}
