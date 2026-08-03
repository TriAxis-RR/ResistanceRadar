Hendriksen Global Sewage selected8 standardized Git release
================================================================

Dataset
-------
Hendriksen Global Sewage selected8 subset used for AMR read mapping,
DIAMOND/CARD comparison, metagenomic assembly, mobile-element
classification, and ARG–mobility benchmarking.

Samples
-------
ERR1713350
ERR1713354
ERR1713355
ERR1713359
ERR1713361
ERR1713362
ERR1713365
ERR1713369

Directory structure
-------------------
FASTP/<sample_id>/
    Per-sample FASTP JSON quality-control report.

Bowtie2/<sample_id>/
    Per-sample Bowtie2 decontamination summary.

CARD_RGI_BWT/<sample_id>/
    RGI BWT gene_mapping_data and overall_mapping_stats files.

CARD_RGI_Main/<sample_id>/
    Compact RGI Main contig-level result table.

    The compact table retains every result row and 26 of the original
    29 analytical columns. The following sequence-heavy columns were
    omitted to keep the Git release manageable:

    - Predicted_DNA
    - Predicted_Protein
    - CARD_Protein_Sequence

    Source and delivered checksums, row counts, and column counts are
    recorded in metadata/rgi_main_compact_provenance.tsv.

DIAMOND_CARD/<sample_id>/
    Gzip-compressed merged DIAMOND/CARD result table.

    Files were compressed deterministically with gzip. Integrity and
    source checksums are recorded in
    metadata/diamond_card_compression_provenance.tsv.

    Example decompression:
        gzip -dc <sample_id>_diamond_card_merged.tsv.gz > results.tsv

Assembly/<sample_id>/
    Git-safe MEGAHIT assembly metrics, including contig count, total
    assembly length, N50/L50, N90/L90, GC percentage, and length
    summaries.

    Full contig FASTA files are not included in Git.

geNomAD/<sample_id>/
    Per-sample plasmid and virus summary tables.

ARG_Mobility_Overlap/<sample_id>/
    Selected8-specific strict/perfect overlap, sample-summary, and
    drug-class-summary tables extracted from the completed 16-sample
    analysis.

metadata/
    sample_manifest.tsv
        Selected8 sample identifiers.

    dataset_inputs_index.tsv
        File-level inventory, sizes, checksums, and validation status.

    rgi_main_compact_provenance.tsv
        Provenance for compact RGI Main tables.

    diamond_card_compression_provenance.tsv
        Provenance and integrity validation for compressed
        DIAMOND/CARD tables.

    checksums.sha256
        SHA-256 checksums for the release contents.

Release principles
------------------
- Original pipeline outputs remain untouched in their analysis folders.
- This release contains real files rather than Hazel-specific symlinks.
- Files are organized by analytical stage and then by sample.
- No FASTQ files, full contig FASTA files, software databases, or other
  bulk runtime intermediates are included.
- Large RGI Main sequence fields were excluded, but all result rows and
  non-sequence analytical columns were preserved.
