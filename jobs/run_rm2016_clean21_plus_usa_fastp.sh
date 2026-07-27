#!/bin/bash
#BSUB -J RM2016_C21_USA_FASTP[1-11]%2
#BSUB -q single_chassis
#BSUB -n 2
#BSUB -R "rusage[mem=16]"
#BSUB -M 16
#BSUB -W 24:00
#BSUB -o logs/RM2016_C21_USA_FASTP_%I.out
#BSUB -e logs/RM2016_C21_USA_FASTP_%I.err

set -euo pipefail

source /usr/local/apps/miniconda20240526/etc/profile.d/conda.sh
conda activate /rsstu/users/s/sleblan/MismatchRepair/Project/envs/rr_meta

cd /rsstu/users/s/sleblan/MismatchRepair/Project/pipeline/rr-pipeline

SHEET="rr-analysis/hendriksen_resistancemap_2016/metadata/batches/hendriksen_resistancemap_clean21_plus_USA_chicago_11_fastp.tsv"
OUTROOT="rr-analysis/hendriksen_resistancemap_2016/results/fastp_clean21_plus_USA_chicago_11"
CKPT="rr-analysis/hendriksen_resistancemap_2016/checkpoints"

LINE=$((LSB_JOBINDEX + 1))

SAMPLE=$(awk -F'\t' -v line="$LINE" 'NR==line {print $1}' "$SHEET")
R1=$(awk -F'\t' -v line="$LINE" 'NR==line {print $2}' "$SHEET")
R2=$(awk -F'\t' -v line="$LINE" 'NR==line {print $3}' "$SHEET")

mkdir -p "${OUTROOT}/${SAMPLE}" "$CKPT"

if [[ -s "${CKPT}/${SAMPLE}_RM2016_C21_USA_FASTP_COMPLETE.txt" ]]; then
  echo "Already complete: $SAMPLE"
  exit 0
fi

echo "Running FASTP for $SAMPLE"
echo "R1: $R1"
echo "R2: $R2"

fastp \
  -i "$R1" \
  -I "$R2" \
  -o "${OUTROOT}/${SAMPLE}/${SAMPLE}_R1.fastp.fastq.gz" \
  -O "${OUTROOT}/${SAMPLE}/${SAMPLE}_R2.fastp.fastq.gz" \
  -h "${OUTROOT}/${SAMPLE}/${SAMPLE}.fastp.html" \
  -j "${OUTROOT}/${SAMPLE}/${SAMPLE}.fastp.json" \
  --thread 2 \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 50

test -s "${OUTROOT}/${SAMPLE}/${SAMPLE}_R1.fastp.fastq.gz"
test -s "${OUTROOT}/${SAMPLE}/${SAMPLE}_R2.fastp.fastq.gz"
test -s "${OUTROOT}/${SAMPLE}/${SAMPLE}.fastp.html"
test -s "${OUTROOT}/${SAMPLE}/${SAMPLE}.fastp.json"

echo "$SAMPLE FASTP complete" > "${CKPT}/${SAMPLE}_RM2016_C21_USA_FASTP_COMPLETE.txt"

echo "DONE: $SAMPLE"
