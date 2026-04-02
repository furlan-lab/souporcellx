use anyhow::{bail, Context, Result};
use rust_htslib::bam::{self, Read as BamRead};
use std::path::Path;

/// Re-tag BAM reads with cell barcodes and UMIs parsed from read names.
///
/// After minimap2 realignment, reads have names of the form `readname;CB;UMI`
/// (or `readname;CB` if `no_umi`). This function splits the read name, strips
/// the encoded tags, and reattaches them as BAM auxiliary tags.
///
/// This is a Rust reimplementation of the upstream souporcell `retag.py`.
pub fn retag_bam(
    input_bam: &Path,
    output_bam: &Path,
    no_umi: bool,
    umi_tag: &str,
    cell_tag: &str,
) -> Result<()> {
    let mut reader =
        bam::Reader::from_path(input_bam).with_context(|| format!("cannot open {}", input_bam.display()))?;

    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(output_bam, &header, bam::Format::Bam)
        .with_context(|| format!("cannot create {}", output_bam.display()))?;

    let mut count = 0u64;
    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        result.context("error reading BAM record")?;

        // Clone tokens out of record before mutating it.
        let qname = std::str::from_utf8(record.qname()).unwrap_or("").to_string();
        let tokens: Vec<String> = qname.split(';').map(|s| s.to_string()).collect();

        if no_umi {
            if tokens.len() != 2 {
                bail!(
                    "expected 2 tokens in read name (no_umi mode), got {}: '{}'",
                    tokens.len(),
                    qname
                );
            }
            record.set_qname(tokens[0].as_bytes());
            record.push_aux(cell_tag.as_bytes(), bam::record::Aux::String(&tokens[1]))?;
        } else {
            if tokens.len() != 3 {
                bail!(
                    "expected 3 tokens in read name, got {}: '{}'",
                    tokens.len(),
                    qname
                );
            }
            record.set_qname(tokens[0].as_bytes());
            record.push_aux(cell_tag.as_bytes(), bam::record::Aux::String(&tokens[1]))?;
            record.push_aux(umi_tag.as_bytes(), bam::record::Aux::String(&tokens[2]))?;
        }

        writer.write(&record)?;
        count += 1;
    }

    eprintln!("retagged {} reads to {}", count, output_bam.display());
    Ok(())
}
