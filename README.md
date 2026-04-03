# souporcellx

<p align="center">
  <img src="souporcellx_logo.png" alt="souporcellx logo" width="400">
</p>

A Rust-first orchestrator for [souporcell](https://github.com/wheaton5/souporcell) workflows on **Slurm clusters**. Supports grouped multi-BAM runs, multiple VCF panels, and pinned local builds for Rust tools.

> **Linux only.** souporcellx is designed exclusively for Linux HPC environments running the Slurm workload manager. It generates and submits jobs via `sbatch`.

## Installing Rust

souporcellx is built with Rust. If you don't have Rust installed, install it via [rustup](https://rustup.rs/):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Follow the prompts, then restart your shell (or run `source ~/.cargo/env`). Verify with `rustc --version`.

## Requirements

The following must be available on the cluster `PATH`:
- `sbatch` (Slurm)
- `Clang` for building vendored tools (lib-hts)

Additional requirements for `souporcellx run`, depending on mode:
- with `--vcf-manifest` and no `--remap`: no extra external tools
- with `--remap` (and `--vcf-manifest`): `minimap2`, `samtools`
- with `--remap` and no `--vcf-manifest` (de novo): all of the above plus `freebayes`, `bcftools`, `bgzip`, `tabix`

`souporcellx validate` only checks manifest content and does not perform external tool PATH checks.

## Install

```bash
cargo install --path .
```

## Quick start

```bash
ml Clang/18.1.8-GCCcore-13.3.0 # need this for lib-hts

cd ~/develop/souporcellx # or whereever you clone this repo
cargo install --path .
souporcellx tools fetch       # clone upstream vartrix & souporcell repos into vendor/
souporcellx tools bootstrap   # build vendored tools from source

# standard run (VCF manifest, no remap): no minimap2/freebayes required
souporcellx run --sample-manifest samples.csv --vcf-manifest vcfs.csv --ref genome.fa --workdir /path/to/run --ks 1,2,3,4

# skip coverage filtering if desired
souporcellx run --sample-manifest samples.csv --vcf-manifest vcfs.csv --ref genome.fa --workdir /path/to/run --ks 1,2,3,4 --skip-coverage-filter

# remap + VCF manifest
ml minimap2/2.29-GCCcore-13.3.0
ml SAMtools/1.21-GCC-13.3.0
souporcellx run --sample-manifest samples.csv --vcf-manifest vcfs.csv --ref genome.fa --workdir /path/to/run --ks 1,2,3,4 --remap

# de novo variant calling (remap + freebayes, no VCF manifest needed)
ml freebayes/1.3.2-GCCcore-8.3.0
# bcftools/bgzip/tabix must also be on PATH in this mode
souporcellx run --sample-manifest samples.csv --ref genome.fa --workdir /path/to/run --ks 1,2,3,4 --remap
```

## Commands

```bash
souporcellx tools fetch       # clone upstream repos into vendor/
souporcellx tools bootstrap   # build vendored tools from source
souporcellx tools update      # pull latest sources and rebuild
souporcellx tools show        # print managed tool paths
souporcellx manifest --cellranger-dirs /path/to/sample1 /path/to/sample2  # generate sample manifest from Cell Ranger outputs

souporcellx stage --cellranger-dirs /path/to/sample1 /path/to/sample2 --dest /local/path  # copy Cell Ranger outputs locally

souporcellx validate --sample-manifest samples.csv --vcf-manifest vcfs.csv
souporcellx run --sample-manifest samples.csv --vcf-manifest vcfs.csv --ref genome.fa --workdir /path/to/run --ks 1,2,3,4
souporcellx filter-vcf --vcf input.vcf --bams s1.bam s2.bam --min-cov 20 --output filtered.vcf  # standalone VCF coverage filter
```

By default `run` performs a dry run, printing every `sbatch` command that would be executed. Add `--submit` to actually submit jobs to Slurm.

## VCF coverage filtering

By default, `run` filters each VCF to retain only variants with sufficient read depth in the BAM files before running vartrix. This mirrors the coverage filtering performed by `souporcell_pipeline.py` (`samtools depth` + `bedtools intersect`), but implemented natively in Rust using `rust-htslib` — no `samtools` or `bedtools` required.

For each (VCF, group) pair, a `covfilt_*` Slurm job runs before vartrix. It opens indexed BAMs for all samples in the group, queries depth at each VCF position, and keeps variants where combined depth >= `min_alt + min_ref` (default: 10 + 10 = 20) and < 100,000.

The same `--min-alt` and `--min-ref` values are also passed to the souporcell binary for cell-level locus filtering during clustering.

| Option | Default | Description |
|---|---|---|
| `--min-alt` | `10` | Minimum alt read depth for coverage filtering and souporcell clustering |
| `--min-ref` | `10` | Minimum ref read depth for coverage filtering and souporcell clustering |
| `--skip-coverage-filter` | off | Skip VCF coverage filtering (pass raw VCFs directly to vartrix) |

The `filter-vcf` subcommand can also be used standalone outside the pipeline:

```bash
souporcellx filter-vcf --vcf common_variants.vcf --bams sample1.bam sample2.bam --min-cov 20 --output filtered.vcf

```

**Note:** BAM files must be indexed (`.bai` files present).

## Generating a sample manifest

> **Requires Cell Ranger >= 5.** Earlier versions use a different output layout and are not supported.

The `manifest` command auto-generates a sample manifest from one or more Cell Ranger output directories. It auto-detects whether the output was produced by `cellranger count` or `cellranger multi` and locates the BAM and barcodes files accordingly.

| Layout | BAM | Barcodes |
|---|---|---|
| `cellranger count` | `outs/possorted_genome_bam.bam` | `outs/filtered_feature_bc_matrix/barcodes.tsv.gz` |
| `cellranger multi` | `outs/per_sample_outs/*/count/sample_alignments.bam` | `outs/per_sample_outs/*/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz` |

For `cellranger count`, the `library_id` is derived from the directory name. For `cellranger multi`, each sub-sample under `per_sample_outs/` becomes its own row with the sub-sample folder name as `library_id`.

```bash
# Print to stdout
souporcellx manifest --cellranger-dirs /path/to/sample1 /path/to/sample2

# Write to file with a custom group ID
souporcellx manifest --cellranger-dirs /path/to/sample1 /path/to/sample2 \
    --group-id my_cohort --output samples.csv

# With custom barcode prefixes (one per sample)
souporcellx manifest --cellranger-dirs /path/to/sample1 /path/to/sample2 \
    --group-id my_cohort --prefixes S1 S2 --output samples.csv

```


| Option | Default | Description |
|---|---|---|
| `--cellranger-dirs` | *(required)* | One or more Cell Ranger output directories |
| `--group-id` | `group1` | Group ID assigned to all rows |
| `--prefixes` | *(uses library_id)* | Barcode prefixes, one per sample. Adds a `prefix` column to the output |
| `--output` | stdout | Write manifest to a file |

## Remapping (`--remap`)

The `--remap` flag enables a per-sample BAM remapping stage before vartrix. This mirrors the remapping step in the original `souporcell_pipeline.py`: cell barcodes and UMIs are extracted from BAM tags, reads are realigned with `minimap2`, and tags are reattached.

The remap stage runs as an independent Slurm job per sample, and all downstream jobs depend on it.

**Additional PATH requirements when using `--remap`:**
- `samtools`
- `minimap2`

```bash
# Remap with de novo variant calling (no VCF manifest needed)
souporcellx run --sample-manifest samples.csv \
    --ref genome.fa --workdir /path/to/run --ks 1,2,3,4 \
    --remap --submit

# Remap can also be combined with a VCF manifest if you already have variant panels
souporcellx run --sample-manifest samples.csv --vcf-manifest vcfs.csv \
    --ref genome.fa --workdir /path/to/run --ks 1,2,3,4 \
    --remap --submit
```

| Option | Default | Description |
|---|---|---|
| `--remap` | off | Enable remapping (renamer -> minimap2 -> retag -> sort/index) |
| `--remap-threads` | `24` | Threads for minimap2 and samtools during remapping |
| `--umi-tag` | `UB` | BAM tag for UMIs |
| `--cell-tag` | `CB` | BAM tag for cell barcodes |
| `--no-umi` | off | Set if BAM files lack UMI tags |

## De novo variant calling (`--remap` without `--vcf-manifest`)

When `--remap` is used **without** a VCF manifest, souporcellx runs de novo variant discovery with `freebayes` instead of using a precomputed VCF panel. This is the fully self-contained mode -- no external VCF is needed.

De novo variant calling **without** remapping is not supported. If `--remap` is not set, you must provide `--vcf-manifest`.

The pipeline in this mode:
1. **Remap** each sample (per-sample Slurm jobs)
2. **Freebayes** per group -- merges the group's remapped BAMs, splits the genome into regions, runs freebayes in parallel, and produces a bgzipped + indexed VCF (`variants.vcf.gz`)
3. **Vartrix** -> **combine** -> **souporcell** -> **troublet** (same as standard mode, using the freebayes-discovered VCF)

**Additional PATH requirements for de novo mode:**
- Everything required by `--remap` (above), plus:
- `freebayes`
- `bcftools`
- `bgzip`
- `tabix`

```bash
# De novo variant calling (no VCF manifest needed)
souporcellx run --sample-manifest samples.csv \
    --ref genome.fa --workdir /path/to/run --ks 1,2,3,4 \
    --remap --submit
```

The freebayes step uses `--min-cov` (derived from `--min-alt + --min-ref`, default 20) and runs with the same thread count as `--remap-threads`.

## Souporcell3 mode

The upstream souporcell binary supports a `souporcell3` mode that iterates multiple times to refine clustering results, with bad cluster detection and reinitialization. This is recommended for datasets with a high number of donors (>16).

```bash
# Enable souporcell3 with k harmonic means clustering
souporcellx run --sample-manifest samples.csv --vcf-manifest vcfs.csv \
    --ref genome.fa --workdir /path/to/run --ks 1,2,3,4 \
    --souporcell3 --clustering-method khm --submit
```

| Option | Default | Description |
|---|---|---|
| `--souporcell3` | off | Enable iterative refinement with bad cluster reinitialization |
| `--clustering-method` | *(unset, uses souporcell default: em)* | Clustering method: `em` (expectation maximization) or `khm` (k harmonic means) |

## Sample manifest

```csv
group_id,library_id,bam,barcodes
merge_R1D2R2D1,AML_MRD_R1_D2_A1,/path/sample1.bam,/path/barcodes1.tsv.gz
merge_R1D2R2D1,AML_MRD_R1_D2_A2,/path/sample2.bam,/path/barcodes2.tsv.gz
merge_R1D2R2D1,AML_MRD_R1_D2_B1,/path/sample3.bam,/path/barcodes3.tsv.gz
```

An optional `prefix` column controls the barcode prefix used when combining matrices. If omitted, `library_id` is used as the prefix.

```csv
group_id,library_id,bam,barcodes,prefix
merge_R1D2R2D1,AML_MRD_R1_D2_A1,/path/sample1.bam,/path/barcodes1.tsv.gz,R1D2A1
merge_R1D2R2D1,AML_MRD_R1_D2_A2,/path/sample2.bam,/path/barcodes2.tsv.gz,R1D2A2
merge_R1D2R2D1,AML_MRD_R1_D2_B1,/path/sample3.bam,/path/barcodes3.tsv.gz,R1D2B1
```

## VCF manifest

```csv
vcf_id,vcf_path
kg1k,/path/filtered_2p_1kgenomes_chr.vcf
common_highconf,/path/common_highconf.vcf
```

### VCF requirements

Vartrix only uses **CHROM**, **POS**, **REF**, and **ALT** from each VCF record. All other fields (ID, QUAL, FILTER, INFO, FORMAT/genotypes) are ignored — the FILTER column does not need to contain `PASS`.

Minimum requirements:
- **Biallelic SNPs/indels only** — multi-allelic records (>2 alleles) are silently skipped
- **Valid nucleotide characters in ALT** — records with non-sequence-resolved ALT alleles are skipped
- **CHROM names must match** the reference FASTA and BAM headers (e.g. `chr1` vs `1`)
- **Format**: plain `.vcf` or bgzipped `.vcf.gz` (bgzipped files require a `.tbi` index)

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

## Staging Cell Ranger outputs (`stage`)

> **Note:** This command has nothing to do with souporcell or genotype demultiplexing. It is a general-purpose utility for copying Cell Ranger output files from a remote or shared filesystem to a local destination for use in downstream R or Python workflows.

The `stage` command copies a curated subset of Cell Ranger output files — web summaries, metrics CSVs, count matrices, and optionally VDJ results — to a local directory. Large files such as BAMs, FASTQ files, and BAM indexes are intentionally excluded. The source directory structure is mirrored under `--dest`.

It auto-detects the Cell Ranger layout:

| Layout | Detected by | Files staged |
|---|---|---|
| `cellranger count` (CR 5+) | presence of `outs/possorted_genome_bam.bam` | `outs/*.html`, `outs/*.csv`, `outs/*matrix.h5` |
| `cellranger multi` (CR 6+) | presence of `outs/per_sample_outs/` | per-sample `*.html`, `*.csv`; `count/*matrix.h5`; optionally VDJ |

```bash
# Stage two Cell Ranger runs to a local directory
souporcellx stage \
    --cellranger-dirs /cluster/runs/sample1 /cluster/runs/sample2 \
    --dest /local/data/project1

# Include VDJ output files
souporcellx stage \
    --cellranger-dirs /cluster/runs/sample1 \
    --dest /local/data/project1 \
    --include-vdj

# Also copy souporcell cluster results from a merge directory
souporcellx stage \
    --cellranger-dirs /cluster/runs/sample1 /cluster/runs/sample2 \
    --souporcell-dirs /cluster/runs/merge_group1 \
    --dest /local/data/project1

# Dry run: print what would be copied without copying
souporcellx stage \
    --cellranger-dirs /cluster/runs/sample1 \
    --dest /local/data/project1 \
    --dry-run
```

| Option | Default | Description |
|---|---|---|
| `--cellranger-dirs` | *(required)* | One or more Cell Ranger output directories |
| `--dest` | *(required)* | Local destination directory |
| `--souporcell-dirs` | *(none)* | Directories containing `souporcell_*/` result subdirs; copies `clusters.tsv` and `log.tsv` |
| `--include-vdj` | off | Also copy VDJ files (`.csv`, `.fasta`, `.tsv`, `.json`) |
| `--dry-run` | off | Print files that would be copied without copying them |
