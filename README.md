# souporcellx

A Rust-first orchestrator for souporcell workflows on **Slurm clusters**. Supports grouped multi-BAM runs, multiple VCF panels, and pinned local builds for Rust tools.

> **Linux only.** souporcellx is designed exclusively for Linux HPC environments running the Slurm workload manager. It generates and submits jobs via `sbatch`.

## Requirements

The following must be available on the cluster `PATH`:
- `minimap2`
- `freebayes`
- `sbatch` (Slurm)
- `Clang` for building lib-hts

## Install

```bash
cargo install --path .
```

## Quick start

```bash
cargo install --path .
souporcellx tools fetch       # clone upstream vartrix & souporcell repos into vendor/
souporcellx tools bootstrap   # build vendored tools from source

ml freebayes/1.3.2-GCCcore-8.3.0
ml minimap2/2.29-GCCcore-13.3.0
ml Clang/18.1.8-GCCcore-13.3.0

souporcellx run --sample-manifest samples.csv --vcf-manifest vcfs.csv --ref genome.fa --workdir /path/to/run --ks 1,2,3,4

# skip coverage filtering if desired
souporcellx run --sample-manifest samples.csv --vcf-manifest vcfs.csv --ref genome.fa --workdir /path/to/run --ks 1,2,3,4 --skip-coverage-filter
```

## Commands

```bash
souporcellx tools fetch       # clone upstream repos into vendor/
souporcellx tools bootstrap   # build vendored tools from source
souporcellx tools update      # pull latest sources and rebuild
souporcellx tools show        # print managed tool paths
souporcellx manifest --cellranger-dirs /path/to/sample1 /path/to/sample2  # generate sample manifest from Cell Ranger outputs

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

## real test
export id=AML_MRD_DL3
export wd="/hpc/temp/furlan_s/${id}"
mergename="merge_R1D2R2D1_oldwf"
export OUT=${wd}/${mergename}
cd $OUT
VCF=/fh/fast/furlan_s/grp/refs/vcf/GRCh38/filtered_2p_1kgenomes_chr.vcf
BAM=/hpc/temp/furlan_s/AML_MRD_DL3/merge_R1D2R2D1_oldwf/out.sorted.bam
souporcellx filter-vcf --vcf $VCF --bams $BAM --min-cov 20 --output filtered.vcf
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


#REAL EXAMPLE
ROOT=/home/sfurlan/develop/souporcellx/data
GROUPID=toy
mkdir -p $ROOT/$GROUPID
cd $ROOT/$GROUPID

REF=/fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2024-A/fasta/genome.fa


cat > sample_mani.csv << 'EOL'
group_id,library_id,bam,barcodes,prefix
merge_p1p2,AML_MRD_R1_D2_A1,../toy1.bam,../barcodes1.tsv.gz,AML_MRD_R1_D2_A1_
merge_p1p2,AML_MRD_R2_D1_B2,../toy2.bam,../barcodes2.tsv.gz,AML_MRD_R2_D1_B2_
EOL

cat > vcfs.csv << 'EOL'
vcf_id,vcf_path
kg1k,/fh/fast/furlan_s/grp/refs/vcf/GRCh38/filtered_2p_1kgenomes_chr.vcf
EOL

souporcellx validate --sample-manifest sample_mani.csv --vcf-manifest vcfs.csv

souporcellx run --sample-manifest sample_mani.csv \
                --vcf-manifest vcfs.csv \
                --workdir $ROOT/$GROUPID \
                --ref $REF

# souporcellx manifest --cellranger-dirs $ROOT/AML_MRD_R1_D2_A1 $ROOT/AML_MRD_R1_D2_A2 $ROOT/AML_MRD_R1_D2_B1 $ROOT/AML_MRD_R2_D1_B2 --group-id $GROUPID --output sample_mani.csv

# cat > vcfs.csv << 'EOL'
# vcf_id,vcf_path
# kg1k,/fh/fast/furlan_s/grp/refs/vcf/GRCh38/filtered_2p_1kgenomes_chr.vcf
# kg1k_deep,/hpc/temp/furlan_s/AML_MRD_DL3/merged_vcf/common_variants_grch38_1kg30x_maf2pct.vcf.gz
# EOL



```

## REMOVE THIS 

Will run old implementation of souporcell with mergebams to compare results with the above new souporcellx

```sh

### Create toy dataset for souporcellx
screen

grabnode
24
200
7
N

cd ~/develop/souporcellx
ml R/4.4.2-gfbf-2024a
R
```


```R


#devtools::install_github("furlan-lab/mergebamsR")
library(mergebamsR)

bam1 <- "/fh/working/furlan_s/AML_MRD/AML_MRD1/AML_MRD_R1_D2_A1/outs/per_sample_outs/AML_MRD_R1_D2_A1/count/sample_alignments.bam"
bam2 <- "/fh/working/furlan_s/AML_MRD/AML_MRD1/AML_MRD_R2_D1_B2/outs/per_sample_outs/AML_MRD_R2_D1_B2/count/sample_alignments.bam"

cbs <- read.table("data/good_cb.tsv")$x

#subset to 500 cells
cbs1 <- cbs[grep("AML_MRD_R1_D2_A1", cbs)]
cbs2 <- cbs[grep("AML_MRD_R2_D1_B2", cbs)]

cbs1 <- cbs1[sample(1:length(cbs1), 250)]
cbs2 <- cbs2[sample(1:length(cbs2), 250)]

right <- function(x, n) substr(x, nchar(x) - n + 1, nchar(x))
cbs1 <- right(cbs1, 18)
cbs2 <- right(cbs2, 18)


subsetbam(inputbam = bam1, outputbams = "data/toy1.bam", features = list(cbs1), cores=24, split_bam = TRUE)
subsetbam(inputbam = bam2, outputbams = "data/toy2.bam", features = list(cbs2), cores=24, split_bam = TRUE)

write.table(data.frame(cbs1), gzfile("data/barcodes1.tsv.gz"), col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(data.frame(cbs2), gzfile("data/barcodes2.tsv.gz"), col.names = FALSE, row.names = FALSE, quote = FALSE)

q()
N
```

```sh




```

```sh
#####SOUPORCELL after cellranger complete######
ml Singularity/3.5.3
ml Python/3.9.6-GCCcore-11.2.0
ml SAMtools/1.14-GCC-11.2.0

export id=AML_MRD_DL3
export wd="/hpc/temp/furlan_s/${id}"
mergename="merge_R1D2R2D1_oldwf"
#ls -alh $fq2 #optionally check fastq folder                             
samps=(AML_MRD_R1_D2_A1 AML_MRD_R1_D2_A2 AML_MRD_R1_D2_B1 AML_MRD_R2_D1_B2) #edit this
export OUT=${wd}/${mergename}
mkdir $OUT
cd $OUT
export K=1
labels=()
bamvar=()
bcvar=()
for sample in ${samps[@]}; do
bamvar+=($wd/$sample/outs/per_sample_outs/$sample/count/sample_alignments.bam)
bcvar+=($wd/$sample/outs/per_sample_outs/$sample/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz)
labels+=("${sample}_")
done
export bams=$(IFS=, ; echo "${bamvar[*]}")
export bcs=$(IFS=, ; echo "${bcvar[*]}")
export samples=$(IFS=, ; echo "${labels[*]}")
RUNNUM=$(sbatch -n 1 -c 1 -p campus-new --mem-per-cpu=40000MB --wrap='/home/sfurlan/develop/mergebams/target/release/mergebams \
    --inputs $bams \
    --labels $samples \
    --bcs $bcs \
    --out $OUT')
LASTRUNTEMP="${RUNNUM//'Submitted batch job '}"
export LASTRUN="${LASTRUNTEMP//' on cluster gizmo'}"
RUNNUM=$(sbatch -n 1 -c 24 -p campus-new --dependency=afterok:$LASTRUN --mem-per-cpu=16000MB --wrap='samtools sort -@ 24 out_bam.bam -o out.sorted.bam')
#RUNNUM=$(sbatch -n 1 -c 24 -p campus-new -M gizmo --mem-per-cpu=16000MB --wrap='samtools sort -@ 24 out_bam.bam -o out.sorted.bam')
LASTRUNTEMP="${RUNNUM//'Submitted batch job '}"
export LASTRUN="${LASTRUNTEMP//' on cluster gizmo'}"
RUNNUM=$(sbatch -n 1 -c 12 -p campus-new --dependency=afterok:$LASTRUN --mem-per-cpu=16000MB --wrap='samtools index -@ 12 out.sorted.bam')
LASTRUNTEMP="${RUNNUM//'Submitted batch job '}"
export LASTRUN="${LASTRUNTEMP//' on cluster gizmo'}"


export REF=/fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2020-A/fasta/genome.fa
export VCF=/fh/fast/furlan_s/grp/refs/vcf/GRCh38/filtered_2p_1kgenomes_chr.vcf
export sif=/home/sfurlan/sifs/souporcell.sif
export SINGULARITY_BINDPATH="/fh/scratch,/fh/fast,/shared,/hpc"
export OUTDIR=$OUT/souporcell_1
mkdir -p $OUTDIR
# RUNNUM=$(sbatch -n 1 -c 35 -p campus-new --dependency=afterok:$LASTRUN --mem-per-cpu=16000MB --wrap='singularity exec $sif souporcell_pipeline.py \
#                   -i out.sorted.bam -b out_barcodes.tsv.gz -f $REF --common_variants $VCF --skip_remap True \
#                   -t 35 -o $OUTDIR -k $K')

squeue -u $USER -h -o %i | awk '$1 > 50203326' | xargs scancel

RUNNUM=$(sbatch -n 1 -c 35 -p campus-new --mem-per-cpu=16000MB --wrap='singularity exec $sif souporcell_pipeline.py \
                  -i out.sorted.bam -b out_barcodes.tsv.gz -f $REF --common_variants $VCF --skip_remap True \
                  -t 35 -o $OUTDIR -k $K')
LASTRUNTEMP="${RUNNUM//'Submitted batch job '}"
export LASTRUN="${LASTRUNTEMP//' on cluster gizmo'}"

ks=(1 2 3 4 5 6 7 8)
for k in ${ks[@]}; do
export K=$k
export OUTDIR=${OUT}/souporcell_${K}
mkdir -p $OUTDIR
RUNNUM=$(sbatch -n 1 -c 35 -p campus-new --dependency=afterok:$LASTRUN --mem-per-cpu=16000MB --wrap='singularity exec $sif /opt/souporcell/souporcell/target/release/souporcell -a souporcell_1/alt.mtx -r souporcell_1/ref.mtx -b out_barcodes.tsv.gz --min_alt 2 --min_ref 2 -k $K -t 35 > $OUTDIR/clusters_tmp.tsv 2> $OUTDIR/log.tsv')
LASTRUNTEMP="${RUNNUM//'Submitted batch job '}"
export LASTRUN2="${LASTRUNTEMP//' on cluster gizmo'}"
sbatch -n 1 -c 1 -p campus-new --dependency=afterok:$LASTRUN2 --mem-per-cpu=16000MB --wrap='singularity exec $sif troublet -a souporcell_1/alt.mtx -r souporcell_1/ref.mtx --clusters $OUTDIR/clusters_tmp.tsv > $OUTDIR/clusters.tsv'
done


```

| Option | Default | Description |
|---|---|---|
| `--cellranger-dirs` | *(required)* | One or more Cell Ranger output directories |
| `--group-id` | `group1` | Group ID assigned to all rows |
| `--prefixes` | *(uses library_id)* | Barcode prefixes, one per sample. Adds a `prefix` column to the output |
| `--output` | stdout | Write manifest to a file |

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
