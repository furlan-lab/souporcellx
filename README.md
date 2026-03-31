# souporcellx

A Rust-first orchestrator for souporcell workflows with grouped multi-BAM support, multiple VCF panels, and pinned local builds for Rust tools.

## Opinionated design

This package manages local builds of:
- `vartrix`
- `souporcell`
- `troublet`

It expects these to already exist on the cluster `PATH`:
- `minimap2`
- `freebayes`
- `sbatch` when using submission mode

## Commands

```bash
cargo run -- tools bootstrap
cargo run -- validate --sample-manifest samples.tsv --vcf-manifest vcfs.tsv
cargo run -- run --sample-manifest samples.tsv --vcf-manifest vcfs.tsv --ref genome.fa --workdir /path/to/run --ks 1,2,3,4
```

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

## Vendoring upstream sources

The expected layout is:

```text
vendor/
  vartrix/
  souporcell/
```

Where:
- `vendor/vartrix` is the upstream VarTrix repository root.
- `vendor/souporcell` is the upstream Souporcell repository root containing `souporcell/` and `troublet/` Rust subdirectories.

`tools bootstrap` will build:
- `vendor/vartrix`
- `vendor/souporcell/souporcell`
- `vendor/souporcell/troublet`

and write a small registry file under the local cache directory.
