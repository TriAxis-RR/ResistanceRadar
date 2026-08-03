Hendriksen Global Sewage Pilot8 Standardized Analysis Release
================================================================

Generated: 2026-08-03T22:44:14Z
Samples: 8
Portable stage files: 104
Symbolic links: 0

Samples
-------
ERR1713331
ERR1713332
ERR1713334
ERR1713337
ERR1713339
ERR1713344
ERR1725947
ERR1725962

Release contents
----------------
FASTP
  FASTP JSON quality-control reports.

Bowtie2
  Bowtie2 host/decontamination alignment logs.

CARD_RGI_BWT
  RGI BWT read-level gene-mapping and overall-mapping statistics.

CARD_RGI_Main
  Compact contig-level RGI Main result tables. All result rows are retained.
  The sequence-heavy columns Predicted_DNA, Predicted_Protein, and
  CARD_Protein_Sequence were excluded to create Git-portable tables.
  See metadata/rgi_main_compact_provenance.tsv.

DIAMOND_CARD
  DIAMOND/CARD merged result tables compressed using deterministic gzip.
  See metadata/diamond_card_compression_provenance.tsv.

Assembly
  MEGAHIT assembly summary metrics. Full contig FASTA files are not included
  because they are large pipeline intermediates.

geNomAD
  Per-sample plasmid and virus summary tables.

ARG_Mobility_Overlap
  Per-sample ARG-mobility overlap, strict/perfect overlap summary,
  sample summary, and drug-class summary tables.

Metadata
--------
sample_manifest.tsv
  Pilot8 sample manifest.

dataset_inputs_index.tsv
  One record for each of the 104 delivered stage files, including size,
  SHA-256 checksum, and validation status.

file_checksums.sha256
  SHA-256 checksums for all delivered stage files. Verify from this release
  directory using:

      sha256sum -c metadata/file_checksums.sha256

rgi_main_compact_provenance.tsv
  Source and delivered-file provenance for compact RGI Main tables.

diamond_card_compression_provenance.tsv
  Source and delivered-file provenance for compressed DIAMOND/CARD tables.

Stage file counts
-----------------
FASTP: 8
Bowtie2: 8
CARD_RGI_BWT: 16
CARD_RGI_Main: 8
DIAMOND_CARD: 8
Assembly: 8
geNomAD: 16
ARG_Mobility_Overlap: 32
