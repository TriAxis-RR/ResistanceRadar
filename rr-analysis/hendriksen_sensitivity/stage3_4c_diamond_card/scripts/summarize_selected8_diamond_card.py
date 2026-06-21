#!/usr/bin/env python3

from pathlib import Path
import re
import pandas as pd

ROOT = Path("rr-analysis/hendriksen_sensitivity")

DIAMOND_DIR = ROOT / "stage3_4c_diamond_card" / "results_cleanfasta_selected8"
RGI_GENE_FILE = ROOT / "results" / "rgi_bwt_batch01_selected8_summary" / "normalized_by_clean_read_depth" / "selected8_rgi_gene_level_hits_normalized.tsv"
OUTDIR = ROOT / "stage3_4c_diamond_card" / "results_summary_selected8"
OUTDIR.mkdir(parents=True, exist_ok=True)

diamond_cols = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]

def parse_sseqid(s):
    s = str(s)
    aro_match = re.search(r"ARO:(\d+)", s)
    aro_id = f"ARO:{aro_match.group(1)}" if aro_match else ""
    parts = s.split("|")
    card_hit = parts[-1] if parts else s
    return aro_id, card_hit

all_hits = []
sample_rows = []

for f in sorted(DIAMOND_DIR.glob("*/*_diamond_card_merged.tsv")):
    sample = f.parent.name
    df = pd.read_csv(f, sep="\t", names=diamond_cols, header=None)

    df["sample_id"] = sample
    df[["aro_id", "card_hit_name"]] = df["sseqid"].apply(
        lambda x: pd.Series(parse_sseqid(x))
    )

    total_hits = len(df)
    unique_reads = df["qseqid"].nunique()
    unique_aro = df["aro_id"].replace("", pd.NA).dropna().nunique()
    unique_card_hits = df["card_hit_name"].nunique()
    mean_pident = df["pident"].mean()
    median_pident = df["pident"].median()
    mean_bitscore = df["bitscore"].mean()
    median_bitscore = df["bitscore"].median()

    sample_rows.append({
        "sample_id": sample,
        "diamond_total_hits": total_hits,
        "diamond_unique_reads": unique_reads,
        "diamond_unique_aro_ids": unique_aro,
        "diamond_unique_card_hit_names": unique_card_hits,
        "mean_pident": mean_pident,
        "median_pident": median_pident,
        "mean_bitscore": mean_bitscore,
        "median_bitscore": median_bitscore,
    })

    hit_summary = (
        df.groupby(["sample_id", "aro_id", "card_hit_name"], dropna=False)
          .agg(
              diamond_read_hits=("qseqid", "count"),
              unique_reads=("qseqid", "nunique"),
              mean_pident=("pident", "mean"),
              median_pident=("pident", "median"),
              mean_bitscore=("bitscore", "mean"),
              median_bitscore=("bitscore", "median"),
              max_bitscore=("bitscore", "max"),
          )
          .reset_index()
    )
    all_hits.append(hit_summary)

diamond_hit_summary = pd.concat(all_hits, ignore_index=True)
diamond_sample_summary = pd.DataFrame(sample_rows)

diamond_hit_summary.to_csv(
    OUTDIR / "selected8_diamond_card_hit_summary.tsv",
    sep="\t",
    index=False
)

diamond_sample_summary.to_csv(
    OUTDIR / "selected8_diamond_card_sample_summary.tsv",
    sep="\t",
    index=False
)

# RGI comparison by ARO if normalized RGI gene-level table exists.
if RGI_GENE_FILE.exists():
    rgi = pd.read_csv(RGI_GENE_FILE, sep="\t")

    sample_col = "sample_id"
    rgi_aro_col = "ARO Accession"

    rgi_simple = (
        rgi[[sample_col, rgi_aro_col]]
        .dropna()
        .rename(columns={rgi_aro_col: "aro_id"})
    )

    rgi_simple["aro_id"] = rgi_simple["aro_id"].astype(str)
    rgi_simple["aro_id"] = rgi_simple["aro_id"].str.replace(r"\.0$", "", regex=True)
    rgi_simple["aro_id"] = rgi_simple["aro_id"].apply(
        lambda x: x if x.startswith("ARO:") else f"ARO:{x}"
    )

    diamond_aro = diamond_hit_summary[["sample_id", "aro_id"]].drop_duplicates()
    diamond_aro = diamond_aro[
        diamond_aro["aro_id"].astype(str).str.startswith("ARO:")
    ]

    rows = []
    for sample in sorted(set(diamond_aro["sample_id"]) | set(rgi_simple["sample_id"])):
        dset = set(diamond_aro.loc[diamond_aro["sample_id"] == sample, "aro_id"].dropna())
        rset = set(rgi_simple.loc[rgi_simple["sample_id"] == sample, "aro_id"].dropna())

        overlap = dset & rset
        diamond_only = dset - rset
        rgi_only = rset - dset
        union = dset | rset

        rows.append({
            "sample_id": sample,
            "diamond_unique_aro": len(dset),
            "rgi_unique_aro": len(rset),
            "overlap_aro": len(overlap),
            "diamond_only_aro": len(diamond_only),
            "rgi_only_aro": len(rgi_only),
            "jaccard_concordance": len(overlap) / len(union) if union else 0,
            "overlap_aro_ids": ";".join(sorted(overlap)),
            "diamond_only_aro_ids": ";".join(sorted(diamond_only)),
            "rgi_only_aro_ids": ";".join(sorted(rgi_only)),
        })

    pd.DataFrame(rows).to_csv(
        OUTDIR / "selected8_diamond_card_vs_rgi_aro_concordance.tsv",
        sep="\t",
        index=False
    )
else:
    print(f"WARNING: RGI gene-level file not found: {RGI_GENE_FILE}")



checkpoint = OUTDIR / "SELECTED8_STAGE3_4C_DIAMOND_CARD_SUMMARY_COMPLETE.txt"
checkpoint.write_text(
    "Selected8 Stage 3.4C DIAMOND/CARD summary complete.\n"
    f"Input directory: {DIAMOND_DIR}\n"
    f"Output directory: {OUTDIR}\n"
    "Generated: selected8_diamond_card_hit_summary.tsv\n"
    "Generated: selected8_diamond_card_sample_summary.tsv\n"
    "Generated if possible: selected8_diamond_card_vs_rgi_aro_concordance.tsv\n"
)

print("DONE")
print(f"Wrote outputs to: {OUTDIR}")
