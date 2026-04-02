## Test remap on toy data

Tests the souporcellx remap pipeline on the toy dataset with a VCF manifest.

```sh
screen
grabnode
24
200
7
N

ml SAMtools
ml freebayes/1.3.2-GCCcore-8.3.0
ml minimap2/2.29-GCCcore-13.3.0
ml Clang/18.1.8-GCCcore-13.3.0

ROOT=/home/sfurlan/develop/souporcellx/data
GROUPID=toy_remap
mkdir -p $ROOT/$GROUPID
cd $ROOT/$GROUPID

REF=/fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2024-A/fasta/genome.fa

cat > sample_mani.csv << 'EOL'
group_id,library_id,bam,barcodes,prefix
merge_p1p2,AML_MRD_R1_D2_A1,../toy1_sorted.bam,../barcodes1.tsv.gz,AML_MRD_R1_D2_A1
merge_p1p2,AML_MRD_R2_D1_B2,../toy2_sorted.bam,../barcodes2.tsv.gz,AML_MRD_R2_D1_B2
EOL

cat > vcfs.csv << 'EOL'
vcf_id,vcf_path
kg1k,/fh/fast/furlan_s/grp/refs/vcf/GRCh38/filtered_2p_1kgenomes_chr.vcf
EOL

souporcellx validate --sample-manifest sample_mani.csv --vcf-manifest vcfs.csv

# Dry run first to inspect sbatch commands
souporcellx run --sample-manifest sample_mani.csv \
                --vcf-manifest vcfs.csv \
                --workdir $ROOT/$GROUPID \
                --ref $REF \
                --ks 1,2,3,4 \
                --remap

# Submit
souporcellx run --sample-manifest sample_mani.csv \
                --vcf-manifest vcfs.csv \
                --workdir $ROOT/$GROUPID \
                --ref $REF \
                --ks 1,2,3,4 \
                --remap \
                --submit

# --- Compare remapped vs non-remapped results ---
# After jobs complete, compare cluster assignments
NOREMAP=$ROOT/toy/kg1k/souporcell_2/merge_p1p2/clusters.tsv
REMAP=$ROOT/$GROUPID/kg1k/souporcell_2/merge_p1p2/clusters.tsv

echo "=== no-remap ==="
head $NOREMAP
echo "=== remapped ==="
head $REMAP

# diff cluster assignments (column 2)
paste <(awk '{print $2}' $NOREMAP) <(awk '{print $2}' $REMAP) | awk '$1 != $2 {diff++} END {print "discordant:", diff+0, "/", NR}'

# compare log-likelihoods
echo "no-remap log_loss_singlet:"
awk '{sum += $4} END {print sum/NR}' $NOREMAP
echo "remapped log_loss_singlet:"
awk '{sum += $4} END {print sum/NR}' $REMAP
```


## Real test of filtering vcf


```sh

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

## Test workflow on toy and real data

```sh


#TOY EXAMPLE
ROOT=/home/sfurlan/develop/souporcellx/data
GROUPID=toy
mkdir -p $ROOT/$GROUPID
cd $ROOT/$GROUPID

REF=/fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2024-A/fasta/genome.fa


cat > sample_mani.csv << 'EOL'
group_id,library_id,bam,barcodes,prefix
merge_p1p2,AML_MRD_R1_D2_A1,../toy1_sorted.bam,../barcodes1.tsv.gz,AML_MRD_R1_D2_A1
merge_p1p2,AML_MRD_R2_D1_B2,../toy2_sorted.bam,../barcodes2.tsv.gz,AML_MRD_R2_D1_B2
EOL

cat > vcfs.csv << 'EOL'
vcf_id,vcf_path
kg1k,/fh/fast/furlan_s/grp/refs/vcf/GRCh38/filtered_2p_1kgenomes_chr.vcf
kg1k_deep,/hpc/temp/furlan_s/AML_MRD_DL3/merged_vcf/common_variants_grch38_1kg30x_maf2pct.vcf.gz
EOL

souporcellx validate --sample-manifest sample_mani.csv --vcf-manifest vcfs.csv

souporcellx run --sample-manifest sample_mani.csv \
                --vcf-manifest vcfs.csv \
                --workdir $ROOT/$GROUPID \
                --ref $REF \
                --submit


DEEP_FILT="kg1k/filtered/merge_p1p2/filtered.vcf"
SHALLOW_FILT="kg1k_deep/filtered/merge_p1p2/filtered.vcf"

grep -v '^#' "$DEEP_FILT" | awk '{print $1"\t"$2}' | sort -u > /tmp/deep_filt_pos.txt
grep -v '^#' "$SHALLOW_FILT" | awk '{print $1"\t"$2}' | sort -u > /tmp/shallow_filt_pos.txt

echo "deep total:    $(wc -l < /tmp/deep_filt_pos.txt)"
echo "shallow total: $(wc -l < /tmp/shallow_filt_pos.txt)"
echo "shared:        $(comm -12 /tmp/deep_filt_pos.txt /tmp/shallow_filt_pos.txt | wc -l)"
echo "deep only:     $(comm -23 /tmp/deep_filt_pos.txt /tmp/shallow_filt_pos.txt | wc -l)"
echo "shallow only:  $(comm -13 /tmp/deep_filt_pos.txt /tmp/shallow_filt_pos.txt | wc -l)"

#REAL EXAMPLE
ROOT=/hpc/temp/furlan_s/AML_MRD_DL3
GROUPID=merge_R1D2R2D1
mkdir -p $ROOT/$GROUPID
cd $ROOT/$GROUPID
souporcellx manifest --cellranger-dirs $ROOT/AML_MRD_R1_D2_A1 $ROOT/AML_MRD_R1_D2_A2 $ROOT/AML_MRD_R1_D2_B1 $ROOT/AML_MRD_R2_D1_B2 --group-id $GROUPID --output sample_mani.csv

cat > vcfs.csv << 'EOL'
vcf_id,vcf_path
kg1k,/fh/fast/furlan_s/grp/refs/vcf/GRCh38/filtered_2p_1kgenomes_chr.vcf
kg1k_deep,/hpc/temp/furlan_s/AML_MRD_DL3/merged_vcf/common_variants_grch38_1kg30x_maf2pct.vcf.gz
EOL

souporcellx validate --sample-manifest sample_mani.csv --vcf-manifest vcfs.csv

REF=/fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2024-A/fasta/genome.fa

souporcellx run --sample-manifest sample_mani.csv \
                --vcf-manifest vcfs.csv \
                --workdir $ROOT/$GROUPID \
                --ref $REF \
                --submit
cd /hpc/temp/furlan_s/AML_MRD_DL3/merge_R1D2R2D1

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
ml  SAMtools
samtools sort data/toy1.bam -o data/toy1_sorted.bam
samtools index data/toy1_sorted.bam
samtools sort data/toy2.bam -o data/toy2_sorted.bam
samtools index data/toy2_sorted.bam

rm data/toy1.bam data/toy2.bam

ROOT=/home/sfurlan/develop/souporcellx/data
GROUPID=toy
mkdir -p $ROOT/$GROUPID
export OUT=$ROOT/$GROUPID

ml Singularity/3.5.3
ml Python/3.9.6-GCCcore-11.2.0
ml SAMtools/1.14-GCC-11.2.0

cd $OUT
export K=1

samps=(AML_MRD_R1_D2_A1 AML_MRD_R2_D1_B2) 
labels=(AML_MRD_R1_D2_A1_ AML_MRD_R2_D1_B2_)
bamvar=(../toy1_sorted.bam ../toy2_sorted.bam)
bcvar=(../barcodes1.tsv.gz ../barcodes1.tsv.gz)

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


export REF=/fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2024-A/fasta/genome.fa
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
