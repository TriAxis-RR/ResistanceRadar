Hendriksen Global Sewage Selected11 Standardized Analysis Release
====================================================================

Generated: 2026-08-03T23:11:24Z
Samples: 11
Portable stage files: 143
Symbolic links: 0

Samples
-------
ERR1713348
ERR1713358
ERR1713367
ERR1713372
ERR1713399
ERR1713409
ERR1713410
ERR1725966
ERR1725973
ERR1725975
ERR1725976

Release contents
----------------
FASTP
  FASTP JSON quality-control reports.

Bowtie2
  Bowtie2 host/decontamination alignment logs.

CARD_RGI_BWT
  RGI BWT read-level gene-mapping and overall-mapping statistics.

CARD_RGI_Main
  Compact contig-level RGI Main tables. All result rows are retained.
  Predicted_DNA, Predicted_Protein, and CARD_Protein_Sequence were
  excluded to create Git-portable files.

DIAMOND_CARD
  Deterministically compressed merged DIAMOND/CARD result tables.

Assembly
  MEGAHIT assembly metrics. Full contig FASTA files are not included.

geNomAD
  Per-sample plasmid and virus summary tables.

ARG_Mobility_Overlap
  Per-sample strict/perfect overlap, strict/perfect sample summary,
  overall sample summary, and drug-class summary tables.

Metadata
--------
sample_manifest.tsv
dataset_inputs_index.tsv
file_checksums.sha256
rgi_main_compact_provenance.tsv
diamond_card_compression_provenance.tsv

Stage file counts
-----------------
FASTP: 11
Bowtie2: 11
CARD_RGI_BWT: 22
CARD_RGI_Main: 11
DIAMOND_CARD: 11
Assembly: 11
geNomAD: 22
ARG_Mobility_Overlap: 44
