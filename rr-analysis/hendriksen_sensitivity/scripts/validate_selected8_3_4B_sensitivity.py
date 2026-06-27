#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

TRUTH = Path("rr-analysis/hendriksen_sensitivity/truth_set/AUS18_Hendriksen_reference_drug_classes.tsv")
DETECTED = Path("rr-analysis/hendriksen_sensitivity/results/rgi_bwt_batch01_selected8_summary/normalized_by_clean_read_depth/selected8_rgi_drug_class_summary_normalized.tsv")

OUT_DIR = Path("rr-analysis/hendriksen_sensitivity/results/stage3_4B_selected8_sensitivity")
OUT_DIR.mkdir(parents=True, exist_ok=True)

truth = pd.read_csv(TRUTH, sep="\t")
det = pd.read_csv(DETECTED, sep="\t")

truth_classes = sorted(
    truth["hendriksen_reference_class"]
    .dropna()
    .astype(str)
    .str.strip()
    .unique()
)

samples = sorted(det["sample_id"].dropna().astype(str).unique())

rows = []
detected_rows = []

for sample in samples:
    sample_det = det[det["sample_id"] == sample].copy()
    detected_classes = set(sample_det["drug_class"].dropna().astype(str).str.strip())

    for cls in sorted(detected_classes):
        detected_rows.append({
            "sample_id": sample,
            "detected_drug_class": cls
        })

    for ref_cls in truth_classes:
        recovered = ref_cls in detected_classes
        rows.append({
            "sample_id": sample,
            "hendriksen_reference_class": ref_cls,
            "recovered_by_CARD_RGI": recovered,
            "status": "recovered" if recovered else "missed",
        })

by_class = pd.DataFrame(rows)

by_sample = (
    by_class.groupby("sample_id", as_index=False)
    .agg(
        expected_classes=("hendriksen_reference_class", "nunique"),
        recovered_classes=("recovered_by_CARD_RGI", "sum"),
    )
)

by_sample["missed_classes"] = by_sample["expected_classes"] - by_sample["recovered_classes"]
by_sample["sensitivity_percent"] = (
    by_sample["recovered_classes"] / by_sample["expected_classes"] * 100
)

overall_expected = len(samples) * len(truth_classes)
overall_recovered = int(by_class["recovered_by_CARD_RGI"].sum())
overall_missed = overall_expected - overall_recovered
overall_sensitivity = overall_recovered / overall_expected * 100 if overall_expected else 0

overall = pd.DataFrame([{
    "n_samples": len(samples),
    "reference_classes_per_sample": len(truth_classes),
    "total_expected_class_calls": overall_expected,
    "total_recovered_class_calls": overall_recovered,
    "total_missed_class_calls": overall_missed,
    "pooled_sensitivity_percent": overall_sensitivity,
}])

detected_df = pd.DataFrame(detected_rows)

by_class.to_csv(OUT_DIR / "selected8_3_4B_class_level_recovery_by_sample.tsv", sep="\t", index=False)
by_sample.to_csv(OUT_DIR / "selected8_3_4B_sensitivity_by_sample.tsv", sep="\t", index=False)
overall.to_csv(OUT_DIR / "selected8_3_4B_sensitivity_overall.tsv", sep="\t", index=False)
detected_df.to_csv(OUT_DIR / "selected8_3_4B_detected_drug_classes.tsv", sep="\t", index=False)

missed = by_class[by_class["status"] == "missed"].copy()
missed.to_csv(OUT_DIR / "selected8_3_4B_missed_reference_classes.tsv", sep="\t", index=False)

print("Done.")
print(f"Samples evaluated: {len(samples)}")
print(f"Reference classes per sample: {len(truth_classes)}")
print(f"Total expected class calls: {overall_expected}")
print(f"Recovered class calls: {overall_recovered}")
print(f"Missed class calls: {overall_missed}")
print(f"Pooled sensitivity: {overall_sensitivity:.2f}%")
print(f"Output directory: {OUT_DIR}")
