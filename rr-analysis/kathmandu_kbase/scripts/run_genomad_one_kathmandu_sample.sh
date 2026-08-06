#!/bin/bash
set -euo pipefail

SAMPLE_ID="${1:?Usage: run_genomad_one_kathmandu_sample.sh SAMPLE_ID}"

PROJECT_ROOT="/rsstu/users/s/sleblan/MismatchRepair/Project/pipeline/rr-pipeline"
SHEET="${PROJECT_ROOT}/rr-analysis/kathmandu_kbase/metadata/kathmandu_genomad_13_completed_megahit_samples.tsv"

OUT_ROOT="${PROJECT_ROOT}/rr-analysis/kathmandu_kbase/results/genomad_kathmandu_13_completed_megahit"
CHECKPOINT_ROOT="${PROJECT_ROOT}/rr-analysis/kathmandu_kbase/checkpoints"
DB="/rsstu/users/s/sleblan/MismatchRepair/Project/references/genomad_db/genomad_db"

mkdir -p "$OUT_ROOT" "$CHECKPOINT_ROOT"

source /usr/local/apps/miniconda20240526/etc/profile.d/conda.sh
conda activate /rsstu/users/s/sleblan/MismatchRepair/Project/envs/rr_mobility

cd "$PROJECT_ROOT"

echo "===== geNomAD Kathmandu single-sample run ====="
echo "Sample: $SAMPLE_ID"
echo "Host: $(hostname)"
echo "Start: $(date)"
echo "Project: $PROJECT_ROOT"
echo "Conda env: $CONDA_PREFIX"
echo "geNomAD: $(which genomad)"
genomad --version || true

CONTIGS=$(awk -F'\t' -v S="$SAMPLE_ID" 'NR>1 && $1==S {print $2}' "$SHEET")

if [[ -z "${CONTIGS:-}" ]]; then
  echo "ERROR: Sample $SAMPLE_ID not found in $SHEET"
  exit 1
fi

if [[ ! -s "$CONTIGS" ]]; then
  echo "ERROR: contigs missing or empty: $CONTIGS"
  exit 1
fi

if [[ ! -d "$DB" ]]; then
  echo "ERROR: geNomAD database directory not found: $DB"
  exit 1
fi

OUT_DIR="${OUT_ROOT}/${SAMPLE_ID}"
CHECKPOINT="${CHECKPOINT_ROOT}/${SAMPLE_ID}_GENOMAD_COMPLETE.txt"

if [[ -s "$CHECKPOINT" ]]; then
  echo "Already complete: $SAMPLE_ID"
  exit 0
fi

if [[ -d "$OUT_DIR" ]]; then
  BACKUP="${OUT_DIR}_previous_$(date +%Y%m%d_%H%M%S)"
  echo "Existing geNomAD output found. Moving to: $BACKUP"
  mv "$OUT_DIR" "$BACKUP"
fi

mkdir -p "$OUT_DIR"

echo "CONTIGS: $CONTIGS"
echo "OUT_DIR: $OUT_DIR"
echo "DB: $DB"

genomad end-to-end \
  "$CONTIGS" \
  "$OUT_DIR" \
  "$DB" \
  --threads 4 \
  --cleanup

{
  echo "${SAMPLE_ID} GENOMAD COMPLETE"
  echo "date=$(date)"
  echo "host=$(hostname)"
  echo "contigs=${CONTIGS}"
  echo "out_dir=${OUT_DIR}"
} > "$CHECKPOINT"

echo "Completed sample: $SAMPLE_ID"
echo "Checkpoint: $CHECKPOINT"
echo "End: $(date)"
