#!/usr/bin/env python3

import re
from pathlib import Path
import pandas as pd

BOWTIE_DIR = Path("rr-analysis/hendriksen_sensitivity/results/bowtie2_batch01_selected8")
SUMMARY_DIR = Path("rr-analysis/hendriksen_sensitivity/results/rgi_bwt_batch01_selected8_summary")

SAMPLE_SUMMARY = SUMMARY_DIR / "selected8_rgi_bwt_sample_summary.tsv"
DRUG_SUMMARY = SUMMARY_DIR / "selected8_rgi_bwt_drug_class_summary.tsv"
GENE_HITS = SUMMARY_DIR / "selected8_rgi_bwt_gene_level_hits.tsv"

OUT_DIR = SUMMARY_DIR / "normalized_by_clean_read_depth"
OUT_DIR.mkdir(parents=True, exist_ok=True)

def parse_bowtie2_summary(path):
    text = path.read_text()

    input_pairs = None
    clean_pairs = None
    overall_alignment_rate = None

    m = re.search(r"^(\d+)\s+reads; of these:", text, flags=re.MULTILINE)
    if m:
        input_pairs = int(m.group(1))

    m = re.search(r"^\s*(\d+)\s+\([^)]+\)\s+aligned concordantly 0 times", text, flags=re.MULTILINE)
    if m:
        clean_pairs = int(m.group(1))

    m = re.search(r"([0-9.]+)% overall alignment rate", text)
    if m:
        overall_alignment_rate = float(m.group(1))

    if input_pairs is None or clean_pairs is None:
        raise ValueError(f"Could not parse input/clean pairs from {path}")

    removed_pairs = input_pairs - clean_pairs
    removed_percent = 100.0 * removed_pairs / input_pairs if input_pairs else 0.0

    return {
        "sample_id": path.name.replace(".bowtie2.summary.txt", ""),
        "input_read_pairs": input_pairs,
        "clean_read_pairs": clean_pairs,
        "bowtie2_removed_pairs": removed_pairs,
        "bowtie2_removed_percent_from_concordant_mapping": removed_percent,
        "bowtie2_overall_alignment_rate_percent": overall_alignment_rate,
    }

# Parse Bowtie2 summaries
bt_rows = []
for f in sorted(BOWTIE_DIR.glob("*/*.bowtie2.summary.txt")):
    bt_rows.append(parse_bowtie2_summary(f))

bt = pd.DataFrame(bt_rows).sort_values("sample_id")
bt.to_csv(OUT_DIR / "selected8_bowtie2_clean_read_depth.tsv", sep="\t", index=False)

# Sample-level normalization
sample = pd.read_csv(SAMPLE_SUMMARY, sep="\t")
sample_norm = sample.merge(bt, on="sample_id", how="left")

sample_norm["total_ARG_mapped_reads_per_million_clean_pairs"] = (
    sample_norm["total_mapped_reads"] / sample_norm["clean_read_pairs"] * 1_000_000
)

sample_norm.to_csv(
    OUT_DIR / "selected8_rgi_sample_summary_normalized.tsv",
    sep="\t",
    index=False
)

# Drug-class normalization
drug = pd.read_csv(DRUG_SUMMARY, sep="\t")
drug_norm = drug.merge(bt[["sample_id", "clean_read_pairs"]], on="sample_id", how="left")

drug_norm["drug_class_mapped_reads_per_million_clean_pairs"] = (
    drug_norm["total_mapped_reads"] / drug_norm["clean_read_pairs"] * 1_000_000
)

drug_norm.to_csv(
    OUT_DIR / "selected8_rgi_drug_class_summary_normalized.tsv",
    sep="\t",
    index=False
)

# Gene-level normalization
gene = pd.read_csv(GENE_HITS, sep="\t")
gene_norm = gene.merge(bt[["sample_id", "clean_read_pairs"]], on="sample_id", how="left")

gene_norm["ARO_mapped_reads_per_million_clean_pairs"] = (
    gene_norm["All Mapped Reads"] / gene_norm["clean_read_pairs"] * 1_000_000
)

gene_norm.to_csv(
    OUT_DIR / "selected8_rgi_gene_level_hits_normalized.tsv",
    sep="\t",
    index=False
)

print("Done.")
print(f"Bowtie2 summaries parsed: {len(bt)}")
print(f"Output directory: {OUT_DIR}")
print("Wrote:")
print("  selected8_bowtie2_clean_read_depth.tsv")
print("  selected8_rgi_sample_summary_normalized.tsv")
print("  selected8_rgi_drug_class_summary_normalized.tsv")
print("  selected8_rgi_gene_level_hits_normalized.tsv")
