# souporcellx

A Rust-first orchestrator for souporcell workflows on **Slurm clusters**. Supports grouped multi-BAM runs, multiple VCF panels, and pinned local builds for Rust tools.

> **Slurm required.** souporcellx generates and submits jobs via `sbatch`. It is designed exclusively for HPC environments running the Slurm workload manager.

## Requirements

The following must be available on the cluster `PATH`:
- `minimap2`
- `freebayes`
- `sbatch` (Slurm)

## Install

```bash
cargo install --path .
```

## Quick start

```bash
cargo install --path .
souporcellx tools fetch       # clone upstream vartrix & souporcell repos into vendor/
souporcellx tools bootstrap   # build vendored tools from source
souporcellx run --sample-manifest samples.tsv --vcf-manifest vcfs.tsv --ref genome.fa --workdir /path/to/run --ks 1,2,3,4
```

## Commands

```bash
souporcellx tools fetch       # clone upstream repos into vendor/
souporcellx tools bootstrap   # build vendored tools from source
souporcellx tools update      # pull latest sources and rebuild
souporcellx tools show        # print managed tool paths
souporcellx validate --sample-manifest samples.tsv --vcf-manifest vcfs.tsv
souporcellx run --sample-manifest samples.tsv --vcf-manifest vcfs.tsv --ref genome.fa --workdir /path/to/run --ks 1,2,3,4
```

By default `run` performs a dry run, printing every `sbatch` command that would be executed. Add `--submit` to actually submit jobs to Slurm.

## Sample manifest

```tsv
group_id	library_id	mode	bam	barcodes
merge_R1D2R2D1	AML_MRD_R1_D2_A1	concat	/path/sample1.bam	/path/barcodes1.tsv.gz
merge_R1D2R2D1	AML_MRD_R1_D2_A2	concat	/path/sample2.bam	/path/barcodes2.tsv.gz
merge_R1D2R2D1	AML_MRD_R1_D2_B1	sum	/path/sample3_part1.bam	/path/barcodes3.tsv.gz
merge_R1D2R2D1	AML_MRD_R1_D2_B1	sum	/path/sample3_part2.bam	/path/barcodes3.tsv.gz
```

## VCF manifest

```tsv
vcf_id	vcf_path
kg1k	/path/filtered_2p_1kgenomes_chr.vcf
common_highconf	/path/common_highconf.vcf
```

## Vendored tools

`tools fetch` clones the upstream repositories into `vendor/`:

```text
vendor/
  vartrix/                          # 10XGenomics/vartrix
  souporcell/                       # wheaton5/souporcell
    souporcell/                     #   souporcell Rust crate
    troublet/                       #   troublet Rust crate
```

`tools bootstrap` builds release binaries and registers them under the local cache directory (`~/.cache/souporcellx/tools/`).

`tools update` pulls the latest from both repos and rebuilds.
