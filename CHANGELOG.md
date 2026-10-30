# Changelog

Data releases (`v1.x`). Each section's first paragraph is shown on the site's Downloads page
(`page_tables.py`). A release adds datasets without recalculating earlier ones; a new reference
or pipeline version would mean reprocessing everything and a `v2.0`. The website changes
separately and is tagged by date (`site-YYYY-MM-DD`).

## v1.1 · 2026-10-31

Adds two tissue datasets: Moyo et al. 2020 (15 patients, every leiomyoma MED12 mutant) and Leistico et al. 2021 (20 patients, MED12-mutant and other leiomyomas, with a comparison limited to patients without hormone treatment). The v1.0 numbers are unchanged.

- Samples: GSE128229 (30) and GSE142329 (40), tissue only; the cell experiments in those
  studies (AP-1 and HOXA13 knock-down and over-expression) are not included.
- Processed as their own batch (`results_v1.1/`) with nf-core/rnaseq 3.27.0, GENCODE v44 and
  the v1.0 Salmon index; the batches' gene tables are placed side by side (`merge_salmon.R`).
- Mutation and hormone treatment per patient from the papers' patient tables (Moyo et al.
  Supplementary Table 3; Leistico et al. Table S1); `samples.tsv` gains `patient`, `batch` and
  `hormone_treatment` columns, and a sample sequenced in several runs is one row.
- Datasets and comparisons are listed once, in `datasets.json`.

## v1.0 · 2026-10-30

First release: gene TPM and differential expression for Zuberi et al. 2023 (10 patients) and Buyukcelebi et al. 2023 (MED12-edited cells in 2D and 3D culture).

- GSE204969 (20), GSE226014 (12), GSE226016 (8).
- GEO lists the MED12 genotypes of patients 1007 and 1008 swapped; corrected following the
  study authors (1008 LM is G44D, 1007 LM wild-type).
