#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPO = Path("/rsstu/users/s/sleblan/MismatchRepair/Project/pipeline/rr-pipeline")
BASE = REPO / "rr-analysis/kathmandu_kbase"
OUTDIR = BASE / "plots/arg_types_rgi_bwt"
OUTDIR.mkdir(parents=True, exist_ok=True)

SAMPLES = [
    ("E_1-2", "effluent"),
    ("E_3-4", "effluent"),
    ("E_5-6", "effluent"),
    ("E_7-8", "effluent"),
    ("E_9-10", "effluent"),
    ("E_11-12", "effluent"),
    ("I_1-2", "influent"),
    ("I_3-4", "influent"),
    ("I_5-6", "influent"),
    ("I_7-8", "influent"),
    ("I_9-10", "influent"),
    ("I_11-12", "influent"),
    ("RW_1", "river"),
    ("RW_2", "river"),
    ("RW_3", "river"),
    ("RW_4", "river"),
    ("RW_5", "river"),
    ("RW_6", "river"),
]


def first_existing(paths):
    for p in paths:
        if p.exists():
            return p
    return None


def rgi_gene_path(sample):
    candidates = [
        BASE / f"results/rgi_bwt_batch01_2sample/{sample}_rgi_bwt.gene_mapping_data.txt",
        BASE / f"results/rgi_bwt_batch02_ready15/{sample}_rgi_bwt.gene_mapping_data.txt",
        BASE / f"results/rgi_bwt_I_5-6/{sample}_rgi_bwt.gene_mapping_data.txt",
    ]
    return first_existing(candidates)


def safe_float(x):
    try:
        if x is None:
            return 0.0
        x = str(x).replace(",", "").strip()
        if x in {"", "N/A", "NA", "no data"}:
            return 0.0
        return float(x)
    except Exception:
        return 0.0


def split_terms(x):
    if x is None:
        return ["Unknown"]
    x = str(x).strip()
    if not x or x in {"N/A", "NA", "no data"}:
        return ["Unknown"]

    # CARD often separates multiple classes/mechanisms with semicolon.
    parts = []
    for block in x.split(";"):
        block = block.strip()
        if block:
            parts.append(block)

    return parts if parts else ["Unknown"]


records = []

for sample, matrix in SAMPLES:
    path = rgi_gene_path(sample)
    if path is None:
        print(f"MISSING RGI gene mapping: {sample}")
        continue

    with open(path, "r", errors="ignore") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            aro = row.get("ARO Term", "").strip()
            family = row.get("AMR Gene Family", "").strip()
            drug_class = row.get("Drug Class", "").strip()
            mechanism = row.get("Resistance Mechanism", "").strip()
            mapped_reads = safe_float(row.get("All Mapped Reads"))

            if not aro:
                continue

            records.append({
                "sample_id": sample,
                "matrix": matrix,
                "aro_term": aro,
                "amr_gene_family": family if family else "Unknown",
                "drug_class": drug_class if drug_class else "Unknown",
                "resistance_mechanism": mechanism if mechanism else "Unknown",
                "all_mapped_reads": mapped_reads,
                "source_file": str(path.relative_to(REPO)),
            })


print(f"Loaded RGI gene records: {len(records)}")


# Expanded records: split multi-class fields.
expanded = []
for r in records:
    classes = split_terms(r["drug_class"])
    mechanisms = split_terms(r["resistance_mechanism"])
    families = split_terms(r["amr_gene_family"])

    for c in classes:
        expanded.append({
            **r,
            "category_type": "Drug Class",
            "category": c,
        })

    for m in mechanisms:
        expanded.append({
            **r,
            "category_type": "Resistance Mechanism",
            "category": m,
        })

    for fam in families:
        expanded.append({
            **r,
            "category_type": "AMR Gene Family",
            "category": fam,
        })


# Write long table.
long_tsv = OUTDIR / "kathmandu_rgi_bwt_arg_types_long.tsv"
with open(long_tsv, "w", newline="") as f:
    fieldnames = [
        "sample_id",
        "matrix",
        "category_type",
        "category",
        "aro_term",
        "amr_gene_family",
        "drug_class",
        "resistance_mechanism",
        "all_mapped_reads",
        "source_file",
    ]
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for r in expanded:
        writer.writerow({k: r.get(k, "") for k in fieldnames})

print(f"Wrote: {long_tsv}")


def summarize(category_type):
    by_category_reads = defaultdict(float)
    by_category_genes = defaultdict(int)
    by_sample_category_reads = defaultdict(float)
    by_matrix_category_reads = defaultdict(float)

    for r in expanded:
        if r["category_type"] != category_type:
            continue

        cat = r["category"]
        reads = float(r["all_mapped_reads"])
        sample = r["sample_id"]
        matrix = r["matrix"]

        by_category_reads[cat] += reads
        by_category_genes[cat] += 1
        by_sample_category_reads[(sample, cat)] += reads
        by_matrix_category_reads[(matrix, cat)] += reads

    return by_category_reads, by_category_genes, by_sample_category_reads, by_matrix_category_reads


def write_summary(category_type, by_category_reads, by_category_genes):
    out = OUTDIR / f"kathmandu_rgi_bwt_{category_type.lower().replace(' ', '_')}_summary.tsv"
    with open(out, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["category_type", "category", "gene_rows", "all_mapped_reads"])
        for cat, reads in sorted(by_category_reads.items(), key=lambda x: x[1], reverse=True):
            writer.writerow([category_type, cat, by_category_genes.get(cat, 0), reads])
    print(f"Wrote: {out}")


def plot_top_bar(category_type, by_category_reads, filename, top_n=15):
    items = sorted(by_category_reads.items(), key=lambda x: x[1], reverse=True)[:top_n]
    if not items:
        print(f"No data for {category_type}")
        return

    labels = [x[0] for x in items][::-1]
    values = [x[1] for x in items][::-1]

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.barh(labels, values)
    ax.set_title(f"Top {category_type.lower()} signals in Kathmandu samples")
    ax.set_xlabel("Sum of RGI BWT all mapped reads")
    ax.set_ylabel(category_type)
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()

    png = OUTDIR / f"{filename}.png"
    svg = OUTDIR / f"{filename}.svg"
    fig.savefig(png, dpi=300)
    fig.savefig(svg)
    plt.close(fig)

    print(f"Wrote: {png}")
    print(f"Wrote: {svg}")


def plot_sample_stacked(category_type, by_sample_category_reads, by_category_reads, filename, top_n=8):
    top_categories = [x[0] for x in sorted(by_category_reads.items(), key=lambda x: x[1], reverse=True)[:top_n]]
    sample_order = [s for s, _ in SAMPLES]

    data = []
    for sample in sample_order:
        row = []
        other = 0.0
        for cat in top_categories:
            row.append(by_sample_category_reads.get((sample, cat), 0.0))

        # Add all non-top as Other.
        for (s, cat), reads in by_sample_category_reads.items():
            if s == sample and cat not in top_categories:
                other += reads
        row.append(other)
        data.append(row)

    labels = top_categories + ["Other"]

    fig, ax = plt.subplots(figsize=(12, 6))
    bottoms = [0.0] * len(sample_order)

    for j, label in enumerate(labels):
        vals = [row[j] for row in data]
        ax.bar(sample_order, vals, bottom=bottoms, label=label)
        bottoms = [b + v for b, v in zip(bottoms, vals)]

    ax.set_title(f"{category_type} composition by sample")
    ax.set_ylabel("Sum of RGI BWT all mapped reads")
    ax.set_xlabel("Sample")
    ax.tick_params(axis="x", rotation=60)
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()

    png = OUTDIR / f"{filename}.png"
    svg = OUTDIR / f"{filename}.svg"
    fig.savefig(png, dpi=300)
    fig.savefig(svg)
    plt.close(fig)

    print(f"Wrote: {png}")
    print(f"Wrote: {svg}")


def plot_matrix_stacked(category_type, by_matrix_category_reads, by_category_reads, filename, top_n=8):
    matrices = ["effluent", "influent", "river"]
    top_categories = [x[0] for x in sorted(by_category_reads.items(), key=lambda x: x[1], reverse=True)[:top_n]]

    data = []
    for matrix in matrices:
        row = []
        other = 0.0
        for cat in top_categories:
            row.append(by_matrix_category_reads.get((matrix, cat), 0.0))

        for (m, cat), reads in by_matrix_category_reads.items():
            if m == matrix and cat not in top_categories:
                other += reads
        row.append(other)
        data.append(row)

    labels = top_categories + ["Other"]

    fig, ax = plt.subplots(figsize=(9, 5))
    bottoms = [0.0] * len(matrices)

    for j, label in enumerate(labels):
        vals = [row[j] for row in data]
        ax.bar(matrices, vals, bottom=bottoms, label=label)
        bottoms = [b + v for b, v in zip(bottoms, vals)]

    ax.set_title(f"{category_type} composition by sample type")
    ax.set_ylabel("Sum of RGI BWT all mapped reads")
    ax.set_xlabel("Sample type")
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()

    png = OUTDIR / f"{filename}.png"
    svg = OUTDIR / f"{filename}.svg"
    fig.savefig(png, dpi=300)
    fig.savefig(svg)
    plt.close(fig)

    print(f"Wrote: {png}")
    print(f"Wrote: {svg}")


for category_type in ["Drug Class", "Resistance Mechanism", "AMR Gene Family"]:
    by_cat_reads, by_cat_genes, by_sample_cat_reads, by_matrix_cat_reads = summarize(category_type)

    stem = category_type.lower().replace(" ", "_")

    write_summary(category_type, by_cat_reads, by_cat_genes)

    plot_top_bar(
        category_type,
        by_cat_reads,
        f"01_top_{stem}_rgi_bwt_reads",
        top_n=15,
    )

    plot_sample_stacked(
        category_type,
        by_sample_cat_reads,
        by_cat_reads,
        f"02_sample_composition_{stem}_rgi_bwt_reads",
        top_n=8,
    )

    plot_matrix_stacked(
        category_type,
        by_matrix_cat_reads,
        by_cat_reads,
        f"03_matrix_composition_{stem}_rgi_bwt_reads",
        top_n=8,
    )


print("\nDone.")
print(f"Output folder: {OUTDIR}")
