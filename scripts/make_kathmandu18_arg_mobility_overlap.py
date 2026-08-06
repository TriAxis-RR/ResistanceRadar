#!/usr/bin/env python3

from pathlib import Path
import gzip
import pandas as pd

ROOT = Path("rr-analysis/kathmandu_kbase")

RGI_DIRS = [
    ROOT / "results/rgi_main_contigs_genomad_completed13",
    ROOT / "results/rgi_main_contigs_final5_after_genomad_complete",
]

GENOMAD_DIRS = [
    ROOT / "results/genomad_kathmandu_13_completed_megahit",
    ROOT / "results/genomad_kathmandu_final5_after_megahit_complete",
]

OUTDIR = ROOT / "results/arg_mobility_overlap_kathmandu18"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_ALL = OUTDIR / "kathmandu18_ARG_mobility_overlap.tsv"
OUT_ALL_GZ = OUTDIR / "kathmandu18_ARG_mobility_overlap.tsv.gz"
OUT_STRICT = OUTDIR / "kathmandu18_ARG_mobility_overlap_strict_perfect.tsv"
OUT_SUMMARY = OUTDIR / "kathmandu18_ARG_mobility_sample_summary.tsv"
OUT_SUMMARY_STRICT = OUTDIR / "kathmandu18_ARG_mobility_sample_summary_strict_perfect.tsv"
OUT_DRUG = OUTDIR / "kathmandu18_ARG_mobility_drug_class_summary.tsv"
OUT_CKPT = ROOT / "checkpoints/KATHMANDU18_ARG_MOBILITY_OVERLAP_COMPLETE.txt"

SAMPLES_OLD13 = [
    "E_11-12", "E_1-2", "E_3-4", "E_7-8", "E_9-10",
    "I_3-4", "I_5-6", "I_7-8", "I_9-10",
    "RW_2", "RW_4", "RW_5", "RW_6",
]

SAMPLES_FINAL5 = ["E_5-6", "I_1-2", "I_11-12", "RW_1", "RW_3"]
SAMPLES = SAMPLES_OLD13 + SAMPLES_FINAL5


def find_rgi_txt(sample: str) -> Path:
    candidates = []
    for d in RGI_DIRS:
        candidates.extend(d.glob(f"{sample}/{sample}_rgi_main.txt"))
        candidates.extend(d.glob(f"{sample}/{sample}_rgi_main_contigs.txt"))
    candidates = [p for p in candidates if p.is_file() and p.stat().st_size > 0]
    if not candidates:
        raise FileNotFoundError(f"No RGI-main txt found for {sample}")
    return candidates[0]


def find_genomad_summary(sample: str, kind: str) -> Path | None:
    # kind = plasmid or virus
    for d in GENOMAD_DIRS:
        sample_dir = d / sample
        candidates = list(sample_dir.glob(f"{sample}.contigs_summary/{sample}.contigs_{kind}_summary.tsv"))
        candidates += list(sample_dir.glob(f"*.contigs_summary/*.contigs_{kind}_summary.tsv"))
        candidates = [p for p in candidates if p.is_file() and p.stat().st_size > 0]
        if candidates:
            return candidates[0]
    return None


def clean_col(df: pd.DataFrame, name: str, fallback: str = "") -> pd.Series:
    if name in df.columns:
        return df[name].fillna(fallback).astype(str)
    return pd.Series([fallback] * len(df), index=df.index, dtype=str)


def read_mobile_tables(sample: str):
    plasmid_file = find_genomad_summary(sample, "plasmid")
    virus_file = find_genomad_summary(sample, "virus")

    plasmids = {}
    viruses = {}

    if plasmid_file is not None:
        p = pd.read_csv(plasmid_file, sep="\t")
        if "seq_name" in p.columns:
            score_col = "plasmid_score" if "plasmid_score" in p.columns else None
            for _, row in p.iterrows():
                seq = str(row["seq_name"])
                score = row[score_col] if score_col else ""
                plasmids[seq] = score

    if virus_file is not None:
        v = pd.read_csv(virus_file, sep="\t")
        if "seq_name" in v.columns:
            score_col = "virus_score" if "virus_score" in v.columns else None
            for _, row in v.iterrows():
                seq = str(row["seq_name"])
                score = row[score_col] if score_col else ""
                viruses[seq] = score

    return plasmids, viruses, plasmid_file, virus_file


rows = []
missing = []

for sample in SAMPLES:
    rgi_file = find_rgi_txt(sample)
    rgi = pd.read_csv(rgi_file, sep="\t")

    plasmids, viruses, plasmid_file, virus_file = read_mobile_tables(sample)

    if not plasmids and not viruses:
        missing.append((sample, "no nonempty geNomAD plasmid/virus summary"))

    contigs = clean_col(rgi, "Contig")
    rgi = rgi.copy()

    for idx, row in rgi.iterrows():
        contig = str(row.get("Contig", ""))
        in_plasmid = contig in plasmids
        in_virus = contig in viruses

        if in_plasmid and in_virus:
            context = "plasmid_like;virus_like"
        elif in_plasmid:
            context = "plasmid_like"
        elif in_virus:
            context = "virus_like"
        else:
            context = "chromosome_like"

        rows.append({
            "set": "Kathmandu18",
            "sample_id": sample,
            "contig_id": contig,
            "ORF_ID": row.get("ORF_ID", ""),
            "Best_Hit_ARO": row.get("Best_Hit_ARO", ""),
            "ARO": row.get("ARO", ""),
            "Cut_Off": row.get("Cut_Off", ""),
            "Best_Identities": row.get("Best_Identities", ""),
            "Percentage_Length_of_Reference_Sequence": row.get("Percentage Length of Reference Sequence", ""),
            "Drug_Class": row.get("Drug Class", ""),
            "Resistance_Mechanism": row.get("Resistance Mechanism", ""),
            "AMR_Gene_Family": row.get("AMR Gene Family", ""),
            "Antibiotic": row.get("Antibiotic", ""),
            "matched_to_genomad": bool(in_plasmid or in_virus),
            "genomad_context": context,
            "chromosome_score": "" if (in_plasmid or in_virus) else "1",
            "plasmid_score": plasmids.get(contig, ""),
            "virus_score": viruses.get(contig, ""),
        })

overlap = pd.DataFrame(rows)

# Write full overlap
overlap.to_csv(OUT_ALL, sep="\t", index=False)
with gzip.open(OUT_ALL_GZ, "wt") as gz:
    overlap.to_csv(gz, sep="\t", index=False)

# Strict + Perfect only
strict = overlap[overlap["Cut_Off"].isin(["Strict", "Perfect"])].copy()
strict.to_csv(OUT_STRICT, sep="\t", index=False)


def summarize(df: pd.DataFrame, strict_only: bool = False) -> pd.DataFrame:
    out = []
    for sample, g in df.groupby("sample_id", sort=False):
        total = len(g)
        plasmid = g["genomad_context"].str.contains("plasmid_like", na=False).sum()
        virus = g["genomad_context"].str.contains("virus_like", na=False).sum()
        mobile = ((g["genomad_context"].str.contains("plasmid_like", na=False)) |
                  (g["genomad_context"].str.contains("virus_like", na=False))).sum()
        chrom = (g["genomad_context"] == "chromosome_like").sum()
        matched = g["matched_to_genomad"].sum()
        strict_n = (g["Cut_Off"] == "Strict").sum()
        perfect_n = (g["Cut_Off"] == "Perfect").sum()
        loose_n = (g["Cut_Off"] == "Loose").sum()

        rec = {
            "set": "Kathmandu18",
            "sample_id": sample,
        }

        if strict_only:
            rec["strict_perfect_ARG_hits"] = total
        else:
            rec["total_ARG_hits"] = total
            rec["matched_to_genomad"] = matched

        rec.update({
            "chromosome_like_ARG_hits": chrom,
            "plasmid_like_ARG_hits": plasmid,
            "virus_like_ARG_hits": virus,
            "mobile_context_ARG_hits": mobile,
            "unclassified_ARG_hits": 0,
            "percent_ARG_hits_plasmid_like": (100 * plasmid / total) if total else 0,
            "percent_ARG_hits_mobile_context": (100 * mobile / total) if total else 0,
        })

        if strict_only:
            rec["strict_ARG_hits"] = strict_n
            rec["perfect_ARG_hits"] = perfect_n
        else:
            rec["loose_ARG_hits"] = loose_n
            rec["strict_ARG_hits"] = strict_n
            rec["perfect_ARG_hits"] = perfect_n

        out.append(rec)

    return pd.DataFrame(out)


summary = summarize(overlap, strict_only=False)
summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

summary_strict = summarize(strict, strict_only=True)
summary_strict.to_csv(OUT_SUMMARY_STRICT, sep="\t", index=False)

# Drug-class summary, preserving semicolon-delimited classes as reported
drug_rows = []
for sample, g in overlap.groupby("sample_id", sort=False):
    for _, row in g.iterrows():
        classes = str(row.get("Drug_Class", "")).split(";")
        for cls in classes:
            cls = cls.strip()
            if not cls or cls.lower() == "nan":
                continue
            drug_rows.append({
                "set": "Kathmandu18",
                "sample_id": sample,
                "Drug_Class": cls,
                "ARG_hits": 1,
                "mobile_context": bool(row["matched_to_genomad"]),
                "Cut_Off": row.get("Cut_Off", ""),
            })

drug_df = pd.DataFrame(drug_rows)
if not drug_df.empty:
    drug_summary = (
        drug_df.groupby(["set", "sample_id", "Drug_Class"], as_index=False)
        .agg(
            ARG_hits=("ARG_hits", "sum"),
            mobile_context_ARG_hits=("mobile_context", "sum"),
        )
    )
else:
    drug_summary = pd.DataFrame(columns=["set", "sample_id", "Drug_Class", "ARG_hits", "mobile_context_ARG_hits"])

drug_summary.to_csv(OUT_DRUG, sep="\t", index=False)

with open(OUT_CKPT, "w") as f:
    f.write("Kathmandu18 ARG x mobility overlap complete.\n")
    f.write(f"Samples: {len(SAMPLES)}\n")
    f.write(f"Full overlap rows: {len(overlap)}\n")
    f.write(f"Strict+Perfect rows: {len(strict)}\n")
    f.write(f"Output directory: {OUTDIR}\n")
    if missing:
        f.write("Warnings:\n")
        for sample, msg in missing:
            f.write(f"- {sample}: {msg}\n")

print("DONE")
print(f"Samples: {len(SAMPLES)}")
print(f"Full overlap rows: {len(overlap)}")
print(f"Strict+Perfect rows: {len(strict)}")
print(f"Output directory: {OUTDIR}")
