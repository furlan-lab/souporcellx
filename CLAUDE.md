# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

souporcellx is a Rust CLI that orchestrates souporcell single-cell genotype demultiplexing workflows on **Slurm clusters**. It manages vendored tool builds, validates CSV input manifests, and submits multi-stage Slurm job DAGs with proper dependency chains. It requires `sbatch`, `minimap2`, and `freebayes` on the cluster PATH.

## Build & Install

```bash
cargo install --path .     # Install souporcellx binary to ~/.cargo/bin/
cargo build --release      # Or build without installing (binary at target/release/souporcellx)
cargo test                 # Run tests
cargo fmt                  # Format code
cargo clippy               # Lint
```

Key CLI subcommands:
```bash
souporcellx tools fetch        # Clone upstream vartrix & souporcell repos into vendor/
souporcellx tools bootstrap    # Build vendored tools from source
souporcellx tools update       # Pull latest upstream sources and rebuild
souporcellx tools show         # Print managed tool paths
souporcellx validate --sample-manifest samples.csv --vcf-manifest vcfs.csv
souporcellx run --sample-manifest samples.csv --vcf-manifest vcfs.csv \
    --ref genome.fa --workdir /path/to/run --ks 1,2,3,4 [--submit]
```

The `run` command is **dry-run by default**, printing every sbatch command; pass `--submit` to actually submit Slurm jobs.

## Architecture

Six modules behind a single `main.rs` entry point:

- **cli** — clap derive-based CLI definition (`Cli`, `Commands`, `RunArgs`)
- **manifest** — CSV parsing/validation for sample manifests (group_id, library_id, bam, barcodes) and VCF manifests (vcf_id, vcf_path)
- **toolchain** — Fetches, builds, and updates vendored Rust binaries from `vendor/` into `~/.cache/souporcellx/tools/`, writes a registry CSV, and resolves binary paths at runtime
- **pipeline** — Orchestrates the Slurm job DAG: validates inputs, resolves tools, then submits jobs
- **slurm** — Thin wrapper around `sbatch` for job submission
- **paths** — Cache directory and registry file path helpers

## Pipeline Job DAG

```
With VCF manifest (per VCF panel):
  [covfilt per group] → vartrix per sample → combine per group → souporcell per K → troublet

With --remap + VCF manifest:
  remap per sample (top-level, shared) → [covfilt] → vartrix → combine → souporcell → troublet

With --remap, no VCF manifest (de novo):
  remap per sample → freebayes per group → vartrix → combine → souporcell → troublet
```

Jobs are linked via Slurm `--dependency=afterok:JOBID` chains.

## Vendored Tool Layout

```
vendor/
  vartrix/                          # VarTrix repo root
  souporcell/
    souporcell/                     # souporcell Rust crate
    troublet/                       # troublet Rust crate
```

`tools fetch` clones upstream repos, `tools bootstrap` builds them, `tools update` pulls latest and rebuilds.

## External Dependencies (must be on PATH)

- `minimap2`, `freebayes` — checked at runtime before job submission
- `sbatch` — required when `--submit` is used
- `samtools`, `python3`/`python` — required when `--remap` is used
- `freebayes`, `bcftools`, `bgzip`, `tabix` — required for de novo variant calling (`--remap` without `--vcf-manifest`)
