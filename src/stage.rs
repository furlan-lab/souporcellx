use anyhow::{Context, Result};
use std::path::{Path, PathBuf};

pub struct StageArgs {
    pub cellranger_dirs: Vec<PathBuf>,
    pub dest: PathBuf,
    pub souporcell_dirs: Option<Vec<PathBuf>>,
    pub include_vdj: bool,
    pub dry_run: bool,
}

pub fn run(args: StageArgs) -> Result<()> {
    let mut copies: Vec<(PathBuf, PathBuf)> = Vec::new();

    for cr_dir in &args.cellranger_dirs {
        let canon = cr_dir
            .canonicalize()
            .with_context(|| format!("cannot resolve path: {}", cr_dir.display()))?;
        let dir_name = canon
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| "unknown".to_string());
        let dest_root = args.dest.join(&dir_name);

        let per_sample = canon.join("outs/per_sample_outs");
        if per_sample.is_dir() {
            collect_multi(&canon, &per_sample, &dest_root, args.include_vdj, &mut copies)?;
        } else {
            collect_count(&canon, &dest_root, &mut copies)?;
        }
    }

    if let Some(soup_dirs) = &args.souporcell_dirs {
        for soup_dir in soup_dirs {
            let canon = soup_dir
                .canonicalize()
                .with_context(|| format!("cannot resolve path: {}", soup_dir.display()))?;
            let dir_name = canon
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_else(|| "unknown".to_string());
            let dest_root = args.dest.join(&dir_name);
            collect_souporcell(&canon, &dest_root, &mut copies)?;
        }
    }

    if copies.is_empty() {
        eprintln!("no files found to stage");
        return Ok(());
    }

    let mut count = 0usize;
    for (src, dst) in &copies {
        if args.dry_run {
            println!("{} -> {}", src.display(), dst.display());
        } else {
            if let Some(parent) = dst.parent() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("cannot create directory: {}", parent.display()))?;
            }
            std::fs::copy(src, dst).with_context(|| {
                format!("failed to copy {} -> {}", src.display(), dst.display())
            })?;
            count += 1;
        }
    }

    if args.dry_run {
        eprintln!("dry run: {} files would be staged to {}", copies.len(), args.dest.display());
    } else {
        eprintln!("staged {} files to {}", count, args.dest.display());
    }

    Ok(())
}

/// Copy relevant files from a cellranger count layout:
///   outs/*.html, outs/*.csv, outs/*matrix.h5
fn collect_count(cr_dir: &Path, dest_root: &Path, copies: &mut Vec<(PathBuf, PathBuf)>) -> Result<()> {
    let outs = cr_dir.join("outs");
    if !outs.is_dir() {
        return Ok(());
    }
    for entry in std::fs::read_dir(&outs)
        .with_context(|| format!("cannot read: {}", outs.display()))?
    {
        let path = entry?.path();
        if !path.is_file() {
            continue;
        }
        if is_stageable_outs_file(&path) {
            let rel = path.strip_prefix(cr_dir).expect("path must be under cr_dir").to_path_buf();
            copies.push((path, dest_root.join(rel)));
        }
    }
    Ok(())
}

/// Copy relevant files from a cellranger multi layout under per_sample_outs/:
///   <sample>/*.html, <sample>/*.csv
///   <sample>/count/*.csv, <sample>/count/*matrix.h5
///   <sample>/vdj_*/*.{csv,fasta,tsv,json}  (if include_vdj)
fn collect_multi(
    cr_dir: &Path,
    per_sample: &Path,
    dest_root: &Path,
    include_vdj: bool,
    copies: &mut Vec<(PathBuf, PathBuf)>,
) -> Result<()> {
    let mut entries: Vec<_> = std::fs::read_dir(per_sample)
        .with_context(|| format!("cannot read: {}", per_sample.display()))?
        .collect::<std::result::Result<Vec<_>, _>>()?;
    entries.sort_by_key(|e| e.file_name());

    for entry in entries {
        let sample_dir = entry.path();
        if !sample_dir.is_dir() {
            continue;
        }

        // Files directly in the sample dir
        for f in read_files(&sample_dir)? {
            if is_stageable_outs_file(&f) {
                let rel = f.strip_prefix(cr_dir).expect("must be under cr_dir").to_path_buf();
                copies.push((f, dest_root.join(rel)));
            }
        }

        // Files in count/
        let count_dir = sample_dir.join("count");
        if count_dir.is_dir() {
            for f in read_files(&count_dir)? {
                if is_stageable_count_file(&f) {
                    let rel = f.strip_prefix(cr_dir).expect("must be under cr_dir").to_path_buf();
                    copies.push((f, dest_root.join(rel)));
                }
            }
        }

        // VDJ files
        if include_vdj {
            let mut vdj_entries: Vec<_> = std::fs::read_dir(&sample_dir)
                .with_context(|| format!("cannot read: {}", sample_dir.display()))?
                .collect::<std::result::Result<Vec<_>, _>>()?;
            vdj_entries.sort_by_key(|e| e.file_name());
            for vdj_entry in vdj_entries {
                let vdj_dir = vdj_entry.path();
                if !vdj_dir.is_dir() {
                    continue;
                }
                let name = vdj_dir.file_name().unwrap_or_default().to_string_lossy();
                if !name.starts_with("vdj") {
                    continue;
                }
                for f in read_files(&vdj_dir)? {
                    if is_vdj_file(&f) {
                        let rel = f.strip_prefix(cr_dir).expect("must be under cr_dir").to_path_buf();
                        copies.push((f, dest_root.join(rel)));
                    }
                }
            }
        }
    }
    Ok(())
}

/// Copy clusters.tsv and log.tsv from souporcell_* subdirectories.
fn collect_souporcell(
    soup_dir: &Path,
    dest_root: &Path,
    copies: &mut Vec<(PathBuf, PathBuf)>,
) -> Result<()> {
    let mut entries: Vec<_> = std::fs::read_dir(soup_dir)
        .with_context(|| format!("cannot read: {}", soup_dir.display()))?
        .collect::<std::result::Result<Vec<_>, _>>()?;
    entries.sort_by_key(|e| e.file_name());

    for entry in entries {
        let sub = entry.path();
        if !sub.is_dir() {
            continue;
        }
        let name = sub.file_name().unwrap_or_default().to_string_lossy();
        if !name.starts_with("souporcell") {
            continue;
        }
        for filename in &["clusters.tsv", "log.tsv"] {
            let f = sub.join(filename);
            if f.is_file() {
                let rel = f.strip_prefix(soup_dir).expect("must be under soup_dir").to_path_buf();
                copies.push((f, dest_root.join(rel)));
            }
        }
    }
    Ok(())
}

fn read_files(dir: &Path) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    for entry in std::fs::read_dir(dir)
        .with_context(|| format!("cannot read: {}", dir.display()))?
    {
        let path = entry?.path();
        if path.is_file() {
            files.push(path);
        }
    }
    Ok(files)
}

/// Matches files to stage from outs/ (count layout) or per_sample_outs/<sample>/ root (multi layout).
fn is_stageable_outs_file(path: &Path) -> bool {
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    name.ends_with(".html") || name.ends_with(".csv")
}

/// Matches files to stage from per_sample_outs/<sample>/count/ (multi layout)
/// or outs/ for count layout when looking for matrix files.
fn is_stageable_count_file(path: &Path) -> bool {
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    name.ends_with(".csv") || name.ends_with("matrix.h5")
}

/// Matches VDJ output files.
fn is_vdj_file(path: &Path) -> bool {
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    name.ends_with(".csv")
        || name.ends_with(".fasta")
        || name.ends_with(".tsv")
        || name.ends_with(".json")
}
