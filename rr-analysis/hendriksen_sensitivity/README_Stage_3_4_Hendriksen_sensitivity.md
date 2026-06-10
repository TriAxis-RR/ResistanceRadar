# ResistanceRadar Stage 3.4: Hendriksen CARD/RGI ARG Detection Sensitivity

## Purpose

This folder contains a preliminary one-sample research analysis evaluating CARD/RGI ARG detection sensitivity on the Hendriksen global wastewater metagenomics dataset for the NSF SBIR pitch.

This is not a validated full ResistanceRadar prototype benchmark.

## Sample analyzed

- ENA project: ERP015409
- ENA run accession: ERR1713332
- Hendriksen sample/site ID: AUS.18
- Mapping basis: Australia, Canberra, collection date 02-04-2016

## Input reads

The input reads were the Stage 3.3 Bowtie2-cleaned paired-end FASTQ files generated after host/contamination removal using the production decontamination reference:

- GRCh38.p14
- PhiX174
- pUC19
- pBR322
- UniVec

Cleaned reads retained:

- 13,881,133 paired-end read pairs

## CARD/RGI analysis

- Tool: RGI BWT
- Aligner: KMA
- CARD database version: 4.0.1
- RGI gene-level CARD rows detected: 467
- Unique CARD/RGI drug classes detected: 31

## Hendriksen truth set

- File: 41467_2019_8853_MOESM7_ESM.xlsx
- Sheet: ResFind.Gene.count
- Truth-set rule: ResFinder gene present in AUS.18 if count > 0
- Reference genes present: 150
- Reference ARG drug classes: 11

## Scoring approach

Primary scoring was performed at the ARG drug-class level because class-level comparison is more defensible across ResFinder and CARD/RGI nomenclatures.

CARD/RGI beta-lactam subclasses were normalized to beta-lactam antibiotic for comparison with the Hendriksen reference class set.

## Sensitivity result

Sensitivity was defined as:

Sensitivity = reference ARG drug classes recovered by CARD/RGI / total Hendriksen reference ARG drug classes present.

For ERR1713332 / AUS.18:

- Reference ARG drug classes: 11
- Recovered by CARD/RGI: 11
- Missed: 0
- Preliminary class-level sensitivity: 100.00%

## Main deliverables

- tables/Table_3_4_Hendriksen_CARD_RGI_sensitivity.csv
- figures/Figure_3_4_Hendriksen_CARD_RGI_sensitivity.pdf
- methods/Section_3_4_methods_note.txt
