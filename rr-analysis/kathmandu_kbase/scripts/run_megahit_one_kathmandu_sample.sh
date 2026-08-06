#!/bin/bash
set -euo pipefail

SAMPLE_ID="${1:?Usage: run_megahit_one_kathmandu_sample.sh SAMPLE_ID}"

PROJECT_ROOT="/rsstu/users/s/sleblan/MismatchRepair/Project/pipeline/rr-pipeline"
SHEET="${PROJECT_ROOT}/rr-analysis/kathmandu_kbase/metadata/kathmandu_megahit_remaining_after_diamond_complete_samples.tsv"

OUT_ROOT="${PROJECT_ROOT}/rr-analysis/kathmandu_kbase/results/megahit_kathmandu_remaining_after_diamond_complete"
CHECKPOINT_ROOT="${PROJECT_ROOT}/rr-analysis/kathmandu_kbase/checkpoints"

mkdir -p "$OUT_ROOT" "$CHECKPOINT_ROOT"

source /usr/local/apps/miniconda20240526/etc/profile.d/conda.sh
conda activate /rsstu/users/s/sleblan/MismatchRepair/Project/envs/rr_mobility

cd "$PROJECT_ROOT"

echo "===== MEGAHIT Kathmandu single-sample run ====="
echo "Sample: $SAMPLE_ID"
echo "Host: $(hostname)"
echo "Start: $(date)"
echo "Project: $PROJECT_ROOT"
echo "Conda env: $CONDA_PREFIX"
echo "MEGAHIT: $(which megahit)"
megahit --version || true

R1=$(awk -F'\t' -v S="$SAMPLE_ID" 'NR>1 && $1==S {print $4}' "$SHEET")
R2=$(awk -F'\t' -v S="$SAMPLE_ID" 'NR>1 && $1==S {print $5}' "$SHEET")

if [[ -z "${R1:-}" || -z "${R2:-}" ]]; then
  echo "ERROR: Sample $SAMPLE_ID not found in $SHEET"
  exit 1
fi

if [[ ! -s "$R1" ]]; then
  echo "ERROR: R1 missing or empty: $R1"
  exit 1
fi

if [[ ! -s "$R2" ]]; then
  echo "ERROR: R2 missing or empty: $R2"
  exit 1
fi

OUT_DIR="${OUT_ROOT}/${SAMPLE_ID}"
CHECKPOINT="${CHECKPOINT_ROOT}/${SAMPLE_ID}_MEGAHIT_COMPLETE.txt"

if [[ -s "$CHECKPOINT" && -s "${OUT_DIR}/final.contigs.fa" ]]; then
  echo "Already complete: $SAMPLE_ID"
  exit 0
fi

if [[ -d "$OUT_DIR" && ! -s "${OUT_DIR}/final.contigs.fa" ]]; then
  BACKUP="${OUT_DIR}_failed_$(date +%Y%m%d_%H%M%S)"
  echo "Existing incomplete output found. Moving to: $BACKUP"
  mv "$OUT_DIR" "$BACKUP"
fi

TMP_DIR="${OUT_ROOT}/${SAMPLE_ID}_tmp"
rm -rf "$TMP_DIR"
mkdir -p "$TMP_DIR"

echo "R1: $R1"
echo "R2: $R2"
echo "OUT_DIR: $OUT_DIR"
echo "TMP_DIR: $TMP_DIR"

megahit \
  -1 "$R1" \
  -2 "$R2" \
  -o "$OUT_DIR" \
  --out-prefix "$SAMPLE_ID" \
  --num-cpu-threads 2 \
  --min-contig-len 1000 \
  --tmp-dir "$TMP_DIR"

if [[ ! -s "${OUT_DIR}/final.contigs.fa" && -s "${OUT_DIR}/${SAMPLE_ID}.contigs.fa" ]]; then
  ln -s "${SAMPLE_ID}.contigs.fa" "${OUT_DIR}/final.contigs.fa"
fi

if [[ ! -s "${OUT_DIR}/final.contigs.fa" ]]; then
  echo "ERROR: MEGAHIT finished but final contigs file not found."
  exit 1
fi

CONTIGS=$(grep -c '^>' "${OUT_DIR}/final.contigs.fa" || true)

{
  echo "${SAMPLE_ID} MEGAHIT COMPLETE"
  echo "date=$(date)"
  echo "host=$(hostname)"
  echo "out_dir=${OUT_DIR}"
  echo "contigs=${CONTIGS}"
} > "$CHECKPOINT"

rm -rf "$TMP_DIR"

echo "Completed sample: $SAMPLE_ID"
echo "Contigs: $CONTIGS"
echo "Checkpoint: $CHECKPOINT"
echo "End: $(date)"
