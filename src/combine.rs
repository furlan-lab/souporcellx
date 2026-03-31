use anyhow::{bail, Context, Result};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

struct MtxHeader {
    rows: usize,
    cols: usize,
    entries: usize,
    /// Number of lines before the first data line (banner + comments + dimension line).
    header_lines: usize,
}

fn parse_mtx_header(path: &Path) -> Result<MtxHeader> {
    let file = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut header_lines = 0;

    for line in reader.lines() {
        let line = line?;
        header_lines += 1;
        if line.starts_with('%') || line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() != 3 {
            bail!("invalid Matrix Market dimension line in {}", path.display());
        }
        return Ok(MtxHeader {
            rows: parts[0].parse()?,
            cols: parts[1].parse()?,
            entries: parts[2].parse()?,
            header_lines,
        });
    }
    bail!("empty or invalid Matrix Market file: {}", path.display())
}

/// Horizontally concatenate Matrix Market files (same rows, different columns).
fn combine_single_matrix(
    input_dirs: &[PathBuf],
    output_dir: &Path,
    name: &str,
) -> Result<()> {
    let mut headers = Vec::new();
    for dir in input_dirs {
        let path = dir.join(name);
        headers.push(parse_mtx_header(&path)?);
    }

    let num_rows = headers[0].rows;
    for (i, h) in headers.iter().enumerate() {
        if h.rows != num_rows {
            bail!(
                "{} in {} has {} rows but expected {} (from first input)",
                name,
                input_dirs[i].display(),
                h.rows,
                num_rows
            );
        }
    }

    let total_cols: usize = headers.iter().map(|h| h.cols).sum();
    let total_entries: usize = headers.iter().map(|h| h.entries).sum();

    let out_path = output_dir.join(name);
    let out_file = File::create(&out_path)
        .with_context(|| format!("cannot create {}", out_path.display()))?;
    let mut writer = BufWriter::new(out_file);

    writeln!(writer, "%%MatrixMarket matrix coordinate integer general")?;
    writeln!(writer, "%")?;
    writeln!(writer, "{} {} {}", num_rows, total_cols, total_entries)?;

    let mut col_offset: usize = 0;
    for (i, dir) in input_dirs.iter().enumerate() {
        let path = dir.join(name);
        let file = File::open(&path)?;
        let reader = BufReader::new(file);
        let skip = headers[i].header_lines;

        for (line_idx, line) in reader.lines().enumerate() {
            let line = line?;
            if line_idx < skip {
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 3 {
                let row = parts[0];
                let col: usize = parts[1].parse()?;
                writeln!(writer, "{} {} {}", row, col + col_offset, parts[2])?;
            }
        }
        col_offset += headers[i].cols;
    }

    writer.flush()?;
    Ok(())
}

/// Concatenate barcode files, prefixing each barcode with the sample label.
fn combine_barcodes(
    input_dirs: &[PathBuf],
    labels: &[String],
    output_dir: &Path,
) -> Result<()> {
    let out_path = output_dir.join("barcodes.tsv");
    let out_file = File::create(&out_path)
        .with_context(|| format!("cannot create {}", out_path.display()))?;
    let mut writer = BufWriter::new(out_file);

    for (dir, label) in input_dirs.iter().zip(labels.iter()) {
        let path = dir.join("barcodes.tsv");
        let file =
            File::open(&path).with_context(|| format!("cannot open {}", path.display()))?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let barcode = line?;
            writeln!(writer, "{}_{}", label, barcode)?;
        }
    }

    writer.flush()?;
    Ok(())
}

/// Combine vartrix outputs from multiple samples into a single set of matrices.
pub fn combine_matrices(
    input_dirs: &[PathBuf],
    labels: &[String],
    output_dir: &Path,
) -> Result<()> {
    if input_dirs.is_empty() {
        bail!("no input directories provided");
    }
    if input_dirs.len() != labels.len() {
        bail!(
            "number of inputs ({}) does not match number of labels ({})",
            input_dirs.len(),
            labels.len()
        );
    }

    fs::create_dir_all(output_dir)
        .with_context(|| format!("cannot create {}", output_dir.display()))?;

    combine_single_matrix(input_dirs, output_dir, "alt.mtx")?;
    combine_single_matrix(input_dirs, output_dir, "ref.mtx")?;
    combine_barcodes(input_dirs, labels, output_dir)?;

    println!(
        "combined {} inputs into {}",
        input_dirs.len(),
        output_dir.display()
    );
    Ok(())
}
