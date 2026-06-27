#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

IN_DIR = Path("rr-analysis/hendriksen_sensitivity/rgi_bwt_batch01_selected8")
OUT_DIR = Path("rr-analysis/hendriksen_sensitivity/results/rgi_bwt_batch01_selected8_summary")
OUT_DIR.mkdir(parents=True, exist_ok=True)

files = sorted(IN_DIR.glob("*_rgi_bwt.gene_mapping_data.txt"))

if not files:
    raise SystemExit(f"No RGI gene_mapping_data files found in {IN_DIR}")

rows = []

for f in files:
    sample_id = f.name.replace("_rgi_bwt.gene_mapping_data.txt", "")
    df = pd.read_csv(f, sep="\t")

    needed = [
        "ARO Term",
        "ARO Accession",
        "AMR Gene Family",
        "Drug Class",
        "Resistance Mechanism",
        "All Mapped Reads",
        "Average Percent Coverage",
    ]

    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise SystemExit(f"{sample_id}: missing columns: {missing}")

    keep = df[needed].copy()
    keep.insert(0, "sample_id", sample_id)

    keep["All Mapped Reads"] = pd.to_numeric(keep["All Mapped Reads"], errors="coerce").fillna(0)
    keep["Average Percent Coverage"] = pd.to_numeric(
        keep["Average Percent Coverage"], errors="coerce"
    ).fillna(0)

    rows.append(keep)

all_hits = pd.concat(rows, ignore_index=True)

# Gene-level table
all_hits.to_csv(OUT_DIR / "selected8_rgi_bwt_gene_level_hits.tsv", sep="\t", index=False)

# Drug-class summary
drug_rows = []
for _, r in all_hits.iterrows():
    sample_id = r["sample_id"]
    mapped_reads = r["All Mapped Reads"]
    coverage = r["Average Percent Coverage"]
    aro = r["ARO Term"]
    aro_acc = r["ARO Accession"]

    drug_classes = str(r["Drug Class"]).split(";")
    for dc in drug_classes:
        dc = dc.strip()
        if dc and dc.lower() != "nan":
            drug_rows.append({
                "sample_id": sample_id,
                "drug_class": dc,
                "ARO Term": aro,
                "ARO Accession": aro_acc,
                "All Mapped Reads": mapped_reads,
                "Average Percent Coverage": coverage,
            })

drug_df = pd.DataFrame(drug_rows)

drug_summary = (
    drug_df.groupby(["sample_id", "drug_class"], as_index=False)
    .agg(
        n_ARO_terms=("ARO Term", "nunique"),
        total_mapped_reads=("All Mapped Reads", "sum"),
        mean_percent_coverage=("Average Percent Coverage", "mean"),
    )
    .sort_values(["sample_id", "total_mapped_reads"], ascending=[True, False])
)

drug_summary.to_csv(OUT_DIR / "selected8_rgi_bwt_drug_class_summary.tsv", sep="\t", index=False)

# Resistance mechanism summary
mech_rows = []
for _, r in all_hits.iterrows():
    sample_id = r["sample_id"]
    mapped_reads = r["All Mapped Reads"]
    coverage = r["Average Percent Coverage"]
    aro = r["ARO Term"]
    aro_acc = r["ARO Accession"]

    mechanisms = str(r["Resistance Mechanism"]).split(";")
    for mech in mechanisms:
        mech = mech.strip()
        if mech and mech.lower() != "nan":
            mech_rows.append({
                "sample_id": sample_id,
                "resistance_mechanism": mech,
                "ARO Term": aro,
                "ARO Accession": aro_acc,
                "All Mapped Reads": mapped_reads,
                "Average Percent Coverage": coverage,
            })

mech_df = pd.DataFrame(mech_rows)

mech_summary = (
    mech_df.groupby(["sample_id", "resistance_mechanism"], as_index=False)
    .agg(
        n_ARO_terms=("ARO Term", "nunique"),
        total_mapped_reads=("All Mapped Reads", "sum"),
        mean_percent_coverage=("Average Percent Coverage", "mean"),
    )
    .sort_values(["sample_id", "total_mapped_reads"], ascending=[True, False])
)

mech_summary.to_csv(OUT_DIR / "selected8_rgi_bwt_resistance_mechanism_summary.tsv", sep="\t", index=False)

# Per-sample compact summary
sample_summary = (
    all_hits.groupby("sample_id", as_index=False)
    .agg(
        n_ARO_terms=("ARO Term", "nunique"),
        n_AMR_gene_families=("AMR Gene Family", "nunique"),
        total_mapped_reads=("All Mapped Reads", "sum"),
        mean_percent_coverage=("Average Percent Coverage", "mean"),
    )
    .sort_values("sample_id")
)

sample_summary.to_csv(OUT_DIR / "selected8_rgi_bwt_sample_summary.tsv", sep="\t", index=False)

print("Done.")
print(f"Input files: {len(files)}")
print(f"Output directory: {OUT_DIR}")
print(f"Gene-level rows: {len(all_hits)}")
print(f"Drug-class summary rows: {len(drug_summary)}")
print(f"Mechanism summary rows: {len(mech_summary)}")
