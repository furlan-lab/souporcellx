use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use rust_htslib::bam::{self, Read as BamRead};
use rust_htslib::bgzf;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

const MAX_DEPTH: u32 = 100_000;

/// Open a VCF (plain or gzipped) and return a boxed buffered reader.
fn open_vcf(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let is_gz = path
        .extension()
        .is_some_and(|ext| ext == "gz" || ext == "bgz");
    if is_gz {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Build a map of chrom → set of 0-based positions from a VCF file.
fn collect_vcf_positions(vcf: &Path) -> Result<HashMap<String, HashSet<i64>>> {
    let reader = open_vcf(vcf)?;
    let mut positions: HashMap<String, HashSet<i64>> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let mut cols = line.splitn(3, '\t');
        let chrom = match cols.next() {
            Some(c) if !c.is_empty() => c.to_owned(),
            _ => continue,
        };
        let pos_str = match cols.next() {
            Some(p) => p,
            None => continue,
        };
        let pos: i64 = pos_str
            .parse()
            .with_context(|| format!("invalid POS value '{}' in VCF line", pos_str))?;
        // VCF POS is 1-based; rust-htslib uses 0-based
        positions.entry(chrom).or_default().insert(pos - 1);
    }
    Ok(positions)
}

/// Walk BAM files using pileup and accumulate depth at the requested positions.
///
/// Returns a map of (chrom, 0-based pos) → total depth across all BAMs.
fn accumulate_depth(
    bams: &[PathBuf],
    positions: &HashMap<String, HashSet<i64>>,
) -> Result<HashMap<(String, i64), u32>> {
    let mut depth_map: HashMap<(String, i64), u32> = HashMap::new();

    for bam_path in bams {
        let mut reader = bam::IndexedReader::from_path(bam_path)
            .with_context(|| format!("cannot open indexed BAM: {}", bam_path.display()))?;

        // Build tid → chrom name lookup, filtered to chroms we care about
        let tid_to_chrom: HashMap<u32, String> = {
            let header = reader.header().clone();
            let mut map = HashMap::new();
            for tid in 0..header.target_count() {
                let name = String::from_utf8_lossy(header.tid2name(tid)).to_string();
                if positions.contains_key(&name) {
                    map.insert(tid, name);
                }
            }
            map
        };

        // Iterate chromosome by chromosome so we only scan relevant regions
        for (tid, chrom) in &tid_to_chrom {
            let chrom_positions = match positions.get(chrom) {
                Some(p) => p,
                None => continue,
            };

            // Fetch the entire chromosome
            reader
                .fetch(*tid)
                .with_context(|| format!("failed to fetch tid {} ({})", tid, chrom))?;

            // Pileup iterates every covered position sequentially
            for pileup in reader.pileup() {
                let pileup = pileup?;
                let pos = pileup.pos() as i64;

                if !chrom_positions.contains(&pos) {
                    continue;
                }

                // Count reads, skipping unmapped/secondary/dup/qcfail/supplementary
                let mut count: u32 = 0;
                for alignment in pileup.alignments() {
                    let flags = alignment.record().flags();
                    if flags & (0x4 | 0x100 | 0x200 | 0x400 | 0x800) != 0 {
                        continue;
                    }
                    count += 1;
                }

                *depth_map
                    .entry((chrom.clone(), pos))
                    .or_insert(0) += count;
            }
        }
    }

    Ok(depth_map)
}

/// Filter a VCF to only include variants with sufficient read depth in the
/// provided BAM files.
///
/// Uses a two-pass approach: first collects all VCF positions, then walks each
/// BAM once via pileup to accumulate depth (sequential I/O), then re-reads the
/// VCF and emits passing variants.
pub fn filter_vcf(
    vcf: &Path,
    bams: &[PathBuf],
    min_cov: u32,
    output: &Path,
) -> Result<()> {
    if bams.is_empty() {
        anyhow::bail!("no BAM files provided");
    }

    // Pass 1: collect all variant positions from VCF
    eprintln!("coverage filter: reading VCF positions...");
    let positions = collect_vcf_positions(vcf)?;
    let total_positions: usize = positions.values().map(|s| s.len()).sum();
    eprintln!(
        "coverage filter: {} positions across {} contigs",
        total_positions,
        positions.len()
    );

    // Pass 2: walk BAMs via pileup, accumulate depth at VCF positions
    eprintln!("coverage filter: scanning {} BAM file(s)...", bams.len());
    let depth_map = accumulate_depth(bams, &positions)?;

    // Ensure output directory exists
    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("cannot create output dir: {}", parent.display()))?;
    }

    // Pass 3: re-read VCF and emit variants that pass the depth filter
    let vcf_reader = open_vcf(vcf)?;
    let bgzf_writer = bgzf::Writer::from_path(output)
        .map_err(|e| anyhow::anyhow!("cannot create {}: {}", output.display(), e))?;
    let mut writer = BufWriter::new(bgzf_writer);

    let mut kept = 0u64;
    let mut skipped = 0u64;

    for line in vcf_reader.lines() {
        let line = line?;

        if line.starts_with('#') {
            writeln!(writer, "{}", line)?;
            continue;
        }

        let mut cols = line.splitn(3, '\t');
        let chrom = match cols.next() {
            Some(c) if !c.is_empty() => c,
            _ => continue,
        };
        let pos_str = match cols.next() {
            Some(p) => p,
            None => continue,
        };
        let pos: i64 = pos_str
            .parse()
            .with_context(|| format!("invalid POS value '{}' in VCF line", pos_str))?;
        let pos0 = pos - 1;

        let depth = depth_map
            .get(&(chrom.to_owned(), pos0))
            .copied()
            .unwrap_or(0);

        if depth >= min_cov && depth < MAX_DEPTH {
            writeln!(writer, "{}", line)?;
            kept += 1;
        } else {
            skipped += 1;
        }
    }

    writer.flush()?;
    eprintln!(
        "coverage filter: kept {} variants, skipped {} (min_cov={}, max={})",
        kept, skipped, min_cov, MAX_DEPTH
    );
    Ok(())
}
