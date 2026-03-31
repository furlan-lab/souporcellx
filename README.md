# souporcellx

A Rust-first orchestrator for souporcell workflows on **Slurm clusters**. Supports grouped multi-BAM runs, multiple VCF panels, and pinned local builds for Rust tools.

> **Linux only.** souporcellx is designed exclusively for Linux HPC environments running the Slurm workload manager. It generates and submits jobs via `sbatch`.

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

ml freebayes/1.3.2-GCCcore-8.3.0
ml minimap2/2.29-GCCcore-13.3.0

souporcellx run --sample-manifest samples.csv --vcf-manifest vcfs.csv --ref genome.fa --workdir /path/to/run --ks 1,2,3,4
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
```

By default `run` performs a dry run, printing every `sbatch` command that would be executed. Add `--submit` to actually submit jobs to Slurm.

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
ROOT=/hpc/temp/furlan_s/AML_MRD_DL3
GROUPID=merge_R1D2R2D1
mkdir -p $ROOT/$GROUPID
cd $ROOT/$GROUPID

REF=/fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2024-A/fasta/genome.fa
souporcellx manifest --cellranger-dirs $ROOT/AML_MRD_R1_D2_A1 $ROOT/AML_MRD_R1_D2_A2 $ROOT/AML_MRD_R1_D2_B1 $ROOT/AML_MRD_R2_D1_B2 --group-id $GROUPID --output sample_mani.csv

cat > vcfs.csv << 'EOL'
vcf_id,vcf_path
kg1k,/fh/fast/furlan_s/grp/refs/vcf/GRCh38/filtered_2p_1kgenomes_chr.vcf
kg1k_deep,/hpc/temp/furlan_s/AML_MRD_DL3/merged_vcf/common_variants_grch38_1kg30x_maf2pct.vcf.gz
EOL

cat > vcfs.csv << 'EOL'
vcf_id,vcf_path
kg1k,/fh/fast/furlan_s/grp/refs/vcf/GRCh38/filtered_2p_1kgenomes_chr.vcf
EOL

souporcellx validate --sample-manifest sample_mani.csv --vcf-manifest vcfs.csv

souporcellx run --sample-manifest sample_mani.csv \
                --vcf-manifest vcfs.csv \
                --workdir $ROOT/$GROUPID \
                --ref $REF \
                --submit

```

## REMOVE THIS 

Will run old implementation of souporcell with mergebams to compare results with the above new souporcellx

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
export REF=/shared/biodata/ngs/Reference/10X/refdata-gex-GRCh38-2024-A/fasta/genome.fa
export VCF=/fh/fast/furlan_s/grp/refs/vcf/GRCh38/filtered_2p_1kgenomes_chr.vcf
export sif=/home/sfurlan/sifs/souporcell.sif
export SINGULARITY_BINDPATH="/fh/scratch,/fh/fast,/shared,/hpc"
export OUTDIR=$OUT/souporcell_1
mkdir -p $OUTDIR
RUNNUM=$(sbatch -n 1 -c 35 -p campus-new --dependency=afterok:$LASTRUN --mem-per-cpu=16000MB --wrap='singularity exec $sif souporcell_pipeline.py \
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

Note that vartrix handles vcf.gz with the rust-htslib library; just make sure it is bgzipped and has a .tbi file.

```csv
vcf_id,vcf_path
kg1k,/path/filtered_2p_1kgenomes_chr.vcf
common_highconf,/path/common_highconf.vcf
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
