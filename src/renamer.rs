use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Read barcodes from a plain-text or gzipped file into a HashSet.
fn load_barcodes(path: &Path) -> Result<HashSet<String>> {
    let file = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let reader: Box<dyn BufRead> = if path
        .extension()
        .is_some_and(|ext| ext == "gz" || ext == "bgz")
    {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut barcodes = HashSet::new();
    for line in reader.lines() {
        let line = line?;
        let token = line.split_whitespace().next();
        if let Some(bc) = token {
            barcodes.insert(bc.to_string());
        }
    }
    Ok(barcodes)
}

/// Convert a BAM to gzipped FASTQ with cell barcodes (and UMIs) encoded in read names.
///
/// Each output read name has the form `@readname;CB;UMI` (or `@readname;CB` if `no_umi`).
/// Reads are filtered to those matching the barcode whitelist, and deduplicated by
/// cell_barcode + UMI + position (or cell_barcode + position if `no_umi`).
///
/// This is a Rust reimplementation of the upstream souporcell `renamer.py`.
pub fn rename_bam_to_fastq(
    bam_path: &Path,
    barcodes_path: &Path,
    output: &Path,
    no_umi: bool,
    umi_tag: &str,
    cell_tag: &str,
) -> Result<()> {
    let cell_barcodes = load_barcodes(barcodes_path)?;
    eprintln!(
        "loaded {} barcodes from {}",
        cell_barcodes.len(),
        barcodes_path.display()
    );

    let mut reader =
        bam::Reader::from_path(bam_path).with_context(|| format!("cannot open {}", bam_path.display()))?;

    let out_file =
        File::create(output).with_context(|| format!("cannot create {}", output.display()))?;
    let mut writer = BufWriter::new(GzEncoder::new(out_file, Compression::fast()));

    let mut seen = HashSet::new();
    let mut kept = 0u64;
    let mut total = 0u64;

    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        result.context("error reading BAM record")?;
        total += 1;

        // Skip reads without cell barcode tag.
        let cb = match record.aux(cell_tag.as_bytes()) {
            Ok(bam::record::Aux::String(s)) => s.to_string(),
            _ => continue,
        };

        // Skip secondary and supplementary alignments.
        if record.is_secondary() || record.is_supplementary() {
            continue;
        }

        // Skip reads whose barcode is not in the whitelist.
        if !cell_barcodes.contains(&cb) {
            continue;
        }

        // Deduplicate by barcode + UMI + position (or barcode + position if no_umi).
        let pos = record.pos();
        let dedup_key = if no_umi {
            format!("{}{}", cb, pos)
        } else {
            match record.aux(umi_tag.as_bytes()) {
                Ok(bam::record::Aux::String(s)) => format!("{}{}{}", cb, s, pos),
                _ => continue,
            }
        };
        if !seen.insert(dedup_key) {
            continue;
        }

        // Skip reads with no sequence.
        let seq = record.seq();
        if seq.is_empty() {
            continue;
        }

        // Build FASTQ read name with encoded tags.
        let qname = std::str::from_utf8(record.qname()).unwrap_or("?");
        let qual: String = record.qual().iter().map(|&q| (q + 33) as char).collect();
        let seq_str: String = seq.as_bytes().iter().map(|&b| b as char).collect();

        if no_umi {
            writeln!(writer, "@{};{}", qname, cb)?;
        } else {
            let umi = match record.aux(umi_tag.as_bytes()) {
                Ok(bam::record::Aux::String(s)) => s,
                _ => continue,
            };
            writeln!(writer, "@{};{};{}", qname, cb, umi)?;
        }
        writeln!(writer, "{}", seq_str)?;
        writeln!(writer, "+")?;
        writeln!(writer, "{}", qual)?;

        kept += 1;
    }

    writer.flush()?;
    eprintln!(
        "wrote {} reads from {} total to {}",
        kept,
        total,
        output.display()
    );
    Ok(())
}
