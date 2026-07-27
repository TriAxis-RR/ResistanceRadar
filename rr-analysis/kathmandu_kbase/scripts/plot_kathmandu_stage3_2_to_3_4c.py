#!/usr/bin/env python3

import json
import re
import csv
import gzip
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPO = Path("/rsstu/users/s/sleblan/MismatchRepair/Project/pipeline/rr-pipeline")
BASE = REPO / "rr-analysis/kathmandu_kbase"

OUTDIR = BASE / "plots/stage3_2_to_3_4c"
OUTDIR.mkdir(parents=True, exist_ok=True)

SUMMARY_TSV = OUTDIR / "kathmandu_stage3_2_to_3_4c_summary.tsv"

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


def find_fastp_json(sample):
    roots = [
        BASE / "results/fastp_batch01_2sample",
        BASE / "results/fastp_batch02_7sample",
        BASE / "results/fastp_batch02_interleaved_10sample",
    ]
    hits = []
    for root in roots:
        if root.exists():
            hits.extend(root.glob(f"{sample}/**/*fastp*.json"))
            hits.extend(root.glob(f"{sample}/**/*.json"))
    return hits[0] if hits else None


def parse_fastp_json(path):
    if path is None or not path.exists():
        return None, None, None, None

    try:
        with open(path, "r") as f:
            data = json.load(f)
    except Exception:
        return None, None, None, None

    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})

    before_reads = before.get("total_reads")
    after_reads = after.get("total_reads")
    q30_rate = after.get("q30_rate")

    retention_pct = None
    q30_pct = None
    clean_pairs = None

    if before_reads and after_reads:
        retention_pct = 100.0 * float(after_reads) / float(before_reads)

    if q30_rate is not None:
        q30_pct = 100.0 * float(q30_rate)

    if after_reads is not None:
        clean_pairs = int(after_reads) / 2.0

    return before_reads, after_reads, retention_pct, q30_pct


def bowtie_summary_path(sample):
    candidates = [
        BASE / f"results/bowtie2_batch01_2sample/{sample}/{sample}.bowtie2.summary.txt",
        BASE / f"results/bowtie2_kathmandu_batch02_ready15/{sample}/{sample}.bowtie2.summary.txt",
        BASE / f"results/bowtie2_kathmandu_I_5-6/{sample}/{sample}.bowtie2.summary.txt",
    ]
    return first_existing(candidates)


def parse_bowtie2(path):
    if path is None or not path.exists():
        return None

    txt = path.read_text(errors="ignore")

    # Bowtie2 final line usually contains: "x.xx% overall alignment rate"
    m = re.search(r"([0-9.]+)%\s+overall alignment rate", txt)
    if m:
        return float(m.group(1))

    # fallback for custom summary wording
    m = re.search(r"([0-9.]+)%", txt)
    if m:
        return float(m.group(1))

    return None


def rgi_paths(sample):
    candidates_stats = [
        BASE / f"results/rgi_bwt_batch01_2sample/{sample}_rgi_bwt.overall_mapping_stats.txt",
        BASE / f"results/rgi_bwt_batch02_ready15/{sample}_rgi_bwt.overall_mapping_stats.txt",
        BASE / f"results/rgi_bwt_I_5-6/{sample}_rgi_bwt.overall_mapping_stats.txt",
    ]
    candidates_gene = [
        BASE / f"results/rgi_bwt_batch01_2sample/{sample}_rgi_bwt.gene_mapping_data.txt",
        BASE / f"results/rgi_bwt_batch02_ready15/{sample}_rgi_bwt.gene_mapping_data.txt",
        BASE / f"results/rgi_bwt_I_5-6/{sample}_rgi_bwt.gene_mapping_data.txt",
    ]
    return first_existing(candidates_stats), first_existing(candidates_gene)


def parse_number_from_text(patterns, txt):
    for pat in patterns:
        m = re.search(pat, txt, flags=re.I)
        if m:
            try:
                return float(m.group(1).replace(",", ""))
            except Exception:
                pass
    return None


def parse_rgi_stats(stats_path, gene_path, clean_pairs):
    total_arg_reads = None
    detected_genes = None
    arg_per_million_clean_pairs = None

    if stats_path and stats_path.exists():
        txt = stats_path.read_text(errors="ignore")
        total_arg_reads = parse_number_from_text([
            r"total[_\s-]*ARG[_\s-]*mapped[_\s-]*reads[^0-9]*([0-9,.]+)",
            r"mapped[_\s-]*reads[^0-9]*([0-9,.]+)",
            r"total[^0-9]+([0-9,.]+)",
        ], txt)

    if gene_path and gene_path.exists():
        try:
            with open(gene_path, "r", errors="ignore") as f:
                rows = [line for line in f if line.strip()]
            detected_genes = max(0, len(rows) - 1)
        except Exception:
            detected_genes = None

    if total_arg_reads is not None and clean_pairs and clean_pairs > 0:
        arg_per_million_clean_pairs = total_arg_reads / clean_pairs * 1_000_000.0

    return total_arg_reads, detected_genes, arg_per_million_clean_pairs


def diamond_path(sample):
    candidates = [
        BASE / f"stage3_4c_diamond_card/results_cleanfasta_17sample/{sample}/{sample}_diamond_card_blastx.tsv",
        BASE / f"stage3_4c_diamond_card/results_cleanfasta_I_5-6/{sample}/{sample}_diamond_card_blastx.tsv",
    ]
    return first_existing(candidates)


def count_lines(path):
    if path is None or not path.exists():
        return None
    try:
        if str(path).endswith(".gz"):
            with gzip.open(path, "rt", errors="ignore") as f:
                return sum(1 for _ in f)
        with open(path, "r", errors="ignore") as f:
            return sum(1 for _ in f)
    except Exception:
        return None


def sample_sort_key(row):
    order = {"effluent": 0, "influent": 1, "river": 2}
    return (order.get(row["matrix"], 9), row["sample_id"])


rows = []

for sample, matrix in SAMPLES:
    fastp_json = find_fastp_json(sample)
    before_reads, after_reads, fastp_retention_pct, q30_pct = parse_fastp_json(fastp_json)

    clean_pairs = after_reads / 2.0 if after_reads else None

    bt2_path = bowtie_summary_path(sample)
    bowtie2_alignment_pct = parse_bowtie2(bt2_path)
    bowtie2_removed_pct = bowtie2_alignment_pct

    rgi_stats_path, rgi_gene_path = rgi_paths(sample)
    total_arg_reads, detected_genes, arg_per_million = parse_rgi_stats(
        rgi_stats_path, rgi_gene_path, clean_pairs
    )

    dpath = diamond_path(sample)
    diamond_hits = count_lines(dpath)

    rows.append({
        "sample_id": sample,
        "matrix": matrix,
        "fastp_json": str(fastp_json.relative_to(REPO)) if fastp_json else "",
        "fastp_before_reads": before_reads,
        "fastp_after_reads": after_reads,
        "fastp_retention_pct": fastp_retention_pct,
        "fastp_q30_pct": q30_pct,
        "bowtie2_summary": str(bt2_path.relative_to(REPO)) if bt2_path else "",
        "bowtie2_alignment_or_removed_pct": bowtie2_alignment_pct,
        "rgi_overall_stats": str(rgi_stats_path.relative_to(REPO)) if rgi_stats_path else "",
        "rgi_gene_mapping": str(rgi_gene_path.relative_to(REPO)) if rgi_gene_path else "",
        "rgi_total_arg_mapped_reads": total_arg_reads,
        "rgi_detected_gene_rows": detected_genes,
        "rgi_arg_reads_per_million_clean_pairs": arg_per_million,
        "diamond_output": str(dpath.relative_to(REPO)) if dpath else "",
        "diamond_hit_rows": diamond_hits,
        "fastp_done": 1 if fastp_retention_pct is not None else 0,
        "bowtie2_done": 1 if bowtie2_alignment_pct is not None else 0,
        "rgi_bwt_done": 1 if rgi_gene_path else 0,
        "diamond_done": 1 if diamond_hits is not None else 0,
    })

rows = sorted(rows, key=sample_sort_key)

fieldnames = list(rows[0].keys())
with open(SUMMARY_TSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote: {SUMMARY_TSV}")


def valid_xy(ykey):
    xs, ys, labels = [], [], []
    for r in rows:
        y = r.get(ykey)
        if y is not None:
            labels.append(r["sample_id"])
            xs.append(r["sample_id"])
            ys.append(float(y))
    return xs, ys, labels


def save_bar(ykey, title, ylabel, filename, rotate=True):
    xs, ys, labels = valid_xy(ykey)
    if not xs:
        print(f"Skipping {filename}: no data")
        return

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.bar(xs, ys)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Sample")
    if rotate:
        ax.tick_params(axis="x", rotation=60)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()

    png = OUTDIR / f"{filename}.png"
    svg = OUTDIR / f"{filename}.svg"
    fig.savefig(png, dpi=300)
    fig.savefig(svg)
    plt.close(fig)
    print(f"Wrote: {png}")
    print(f"Wrote: {svg}")


save_bar(
    "fastp_retention_pct",
    "Kathmandu FASTP read retention",
    "Retention after filtering (%)",
    "01_fastp_retention_pct"
)

save_bar(
    "fastp_q30_pct",
    "Kathmandu FASTP Q30 after filtering",
    "Q30 bases after filtering (%)",
    "02_fastp_q30_pct"
)

save_bar(
    "bowtie2_alignment_or_removed_pct",
    "Kathmandu Bowtie2 decontamination alignment rate",
    "Reads aligned to decontamination index (%)",
    "03_bowtie2_alignment_pct"
)

save_bar(
    "rgi_arg_reads_per_million_clean_pairs",
    "Kathmandu CARD/RGI BWT normalized ARG signal",
    "ARG mapped reads per million clean pairs",
    "04_rgi_arg_reads_per_million_clean_pairs"
)

save_bar(
    "diamond_hit_rows",
    "Kathmandu DIAMOND/CARD hit counts",
    "DIAMOND blastx hit rows",
    "05_diamond_hit_rows"
)


# Pipeline completion matrix
stages = ["fastp_done", "bowtie2_done", "rgi_bwt_done", "diamond_done"]
stage_labels = ["FASTP", "Bowtie2", "RGI BWT", "DIAMOND"]
matrix = [[r[s] for s in stages] for r in rows]
sample_labels = [r["sample_id"] for r in rows]

fig, ax = plt.subplots(figsize=(8, 7))
ax.imshow(matrix, aspect="auto")
ax.set_yticks(range(len(sample_labels)))
ax.set_yticklabels(sample_labels)
ax.set_xticks(range(len(stage_labels)))
ax.set_xticklabels(stage_labels)
ax.set_title("Kathmandu pipeline completion status")
for i, r in enumerate(matrix):
    for j, val in enumerate(r):
        ax.text(j, i, "OK" if val else "NA", ha="center", va="center")
fig.tight_layout()
png = OUTDIR / "06_pipeline_completion_matrix.png"
svg = OUTDIR / "06_pipeline_completion_matrix.svg"
fig.savefig(png, dpi=300)
fig.savefig(svg)
plt.close(fig)
print(f"Wrote: {png}")
print(f"Wrote: {svg}")


# Scatter: RGI normalized signal vs DIAMOND hits
scatter_rows = [
    r for r in rows
    if r["rgi_arg_reads_per_million_clean_pairs"] is not None
    and r["diamond_hit_rows"] is not None
]

if scatter_rows:
    fig, ax = plt.subplots(figsize=(6, 5))
    x = [float(r["rgi_arg_reads_per_million_clean_pairs"]) for r in scatter_rows]
    y = [float(r["diamond_hit_rows"]) for r in scatter_rows]
    ax.scatter(x, y)

    for r, xi, yi in zip(scatter_rows, x, y):
        ax.annotate(r["sample_id"], (xi, yi), fontsize=8, xytext=(3, 3), textcoords="offset points")

    ax.set_title("RGI BWT signal vs DIAMOND/CARD hits")
    ax.set_xlabel("RGI ARG reads per million clean pairs")
    ax.set_ylabel("DIAMOND hit rows")
    ax.grid(alpha=0.3)
    fig.tight_layout()

    png = OUTDIR / "07_rgi_vs_diamond_scatter.png"
    svg = OUTDIR / "07_rgi_vs_diamond_scatter.svg"
    fig.savefig(png, dpi=300)
    fig.savefig(svg)
    plt.close(fig)
    print(f"Wrote: {png}")
    print(f"Wrote: {svg}")
else:
    print("Skipping RGI vs DIAMOND scatter: incomplete paired data")


print("\nQuick completion counts:")
for stage in stages:
    print(stage, sum(r[stage] for r in rows), "/", len(rows))

print(f"\nSummary table:\n{SUMMARY_TSV}")
