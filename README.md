# MED12 RNA-seq explorer

Five public MED12 RNA-seq datasets reprocessed from raw reads the same way, with an
interactive site: **https://brandonlukas.github.io/med12-tpm/**

- **Overview**: a gene search, an illustration of where the samples come from, and why the datasets belong together
- **Genes**: TPM of any gene (or several side by side) in every sample, with DESeq2 results
- **Comparisons**: volcano plot, ranked table and gene lookup for every DESeq2 comparison,
  and agreement between any two comparisons
- **Downloads**: two Excel workbooks (gene TPM; DESeq2 results, one sheet per comparison) and the individual CSV tables
- **Methods**: how everything was processed

| Dataset | Paper | GEO | Samples | Reads |
|---|---|---|---|---|
| Tissue · Zuberi 2023 | Zuberi et al. 2023 JCI Insight | GSE204969 | 10 patients × MM/LM = 20 | SE 75 bp |
| Tissue · Moyo 2020 | Moyo et al. 2020 Nat Commun | GSE128229 | 15 patients × MM/LM = 30 | PE 2×42 |
| Tissue · Leistico 2021 | Leistico et al. 2021 Cell Rep | GSE142329 | 20 patients × MM/LM = 40 | PE 2×42 |
| 2D cells · Buyukcelebi 2023 | Buyukcelebi et al. 2023 Nat Commun | GSE226014 | WT, KO, Mut1, Mut2 × 3 = 12 | PE 2×75 |
| 3D spheroids · Buyukcelebi 2023 | same | GSE226016 | WT ×2, Mut1 ×3, Mut2 ×3 = 8 | PE 2×151 |

**Genotype correction:** GEO has the MED12 genotypes of patients 1007 and 1008 swapped.
1008 LM is the G44D mutant and 1007 LM is wild-type. `samples.tsv` uses the corrected
genotype and keeps GEO's label in `genotype_as_listed_on_GEO_(incorrect)`.

Gene TPM and DESeq2 tables are attached to the
[releases](https://github.com/brandonlukas/med12-tpm/releases); `CHANGELOG.md` lists what each added.
The [methods page](https://brandonlukas.github.io/med12-tpm/methods.html) documents every step.

## Files

- `datasets.json`: the one list of datasets and comparisons, their labels and papers, the release, and each batch's output folder
- `samples.tsv`: one row per sample: dataset, batch, group, patient, mutation, hormone treatment, GSM/SRR IDs, expected read counts
- `CHANGELOG.md`: what each data release added
- `samplesheet.csv`: nf-core/rnaseq samplesheet (strandedness auto-detected)
- `params.yaml`: pipeline parameters (Salmon-only)
- `home.config`: CPU/RAM limits; edit for your machine
- `fetch_reference.sh`: GENCODE v44 fasta + GTF, md5-checked
- `download_fastqs.sh`: runs from NCBI SRA's AWS mirror, read counts checked; optionally only some datasets
- `run_rnaseq.sh`: launches nf-core/rnaseq 3.27.0, optionally on a subset of samples into its own folder; `share.config` caps it to share the machine
- `merge_salmon.R`: places the batches' gene tables side by side (`results_merged/`)
- `split_tpm.py`: per-dataset TPM tables and an Excel workbook
- `deseq2.R`: DESeq2 within each dataset
- `excel_deseq2.py`: one Excel workbook of all DESeq2 comparisons, with a README sheet
- `pca.R`, `pca_svg.py`: per-dataset PCA (DESeq2 vst, top 500 variable genes) → `docs/data/pca.json` → the panels on the Methods page
- `build_site_data.py`: data for the site (`docs/data/`)
- `page_tables.py`: the dataset tables on the Methods and Downloads pages, the version list and the footers
- `illustration/vector.py`: generates the overview's illustration (smooth vector SVG, gradient shading) into `docs/index.html`
- `stamp.py`: stamps the pages' stylesheet and script links with a content hash; run it before committing a change to `docs/`
- `tests/site_test.py`: scripted browser checks of the site (Playwright); `tests/screenshots.py` renders pages to PNG
- `docs/`: the site, served by GitHub Pages: `index.html` + `overview.js` (Overview), `genes.html` + `app.js`
  (Genes), `de.html` + `de.js` (Comparisons), `methods.html`, `checks.html` (evidence for the sample checks, linked from Methods), `common.js`, `style.css`; Plotly 2.35.2 is vendored in `docs/vendor/`
  so the site has no third-party dependency at load time

## Reproducing (Linux / WSL2)

Needs Nextflow ≥ 25.10.4, Docker, sra-tools, pigz, R with DESeq2, apeglm, data.table and
jsonlite, and Python with pandas and openpyxl; about 400 GB of free disk.

```bash
./fetch_reference.sh            # ~1 GB
./download_fastqs.sh 6          # ~200 GB
./run_rnaseq.sh                                      # batch v1.0 into results/, builds the Salmon index
ONLY='^(moyo|leistico)_' OUTDIR=results_v1.1 ./run_rnaseq.sh \
    --salmon_index results/genome/index/salmon \
    --transcript_fasta results/genome/genome.transcripts.fa   # batch v1.1, same index
Rscript merge_salmon.R          # results_merged/
python3 split_tpm.py            # deliverables/*_gene_TPM.csv, .xlsx
Rscript deseq2.R                # deliverables/de/*.csv
python3 excel_deseq2.py         # deliverables/MED12_RNAseq_DESeq2.xlsx
Rscript pca.R && python3 pca_svg.py   # docs/data/pca.json and the Methods PCA panels
python3 build_site_data.py      # docs/data/
python3 page_tables.py          # dataset tables, versions and footers in docs/
```

## Adding a dataset

Earlier datasets are never re-run: TPM is per sample, and DESeq2, PCA and percentiles are
computed within a dataset.

1. Add its samples to `samples.tsv` with a new `batch` (e.g. `v1.2`), its entry to
   `datasets.json` (groups, comparisons, labels), the batch's folder to `datasets.json`'s
   `batches`, its rows to `samplesheet.csv`, and a section to `CHANGELOG.md`; set `release`.
2. `./download_fastqs.sh 6 <dataset_id>`
3. `ONLY='<regex on its sample names>' OUTDIR=results_v1.2 ./run_rnaseq.sh --salmon_index … --transcript_fasta …`
   with the same index, so every batch shares one gene list (`merge_salmon.R` checks this).
4. `Rscript merge_salmon.R`, then the rest of the commands above; `Rscript deseq2.R
   results_merged/salmon.merged.gene_counts_length_scaled.tsv <dataset_id>` runs only the new one.
5. A new group that needs a color goes in `COLORS` in `docs/common.js` and `pca_svg.py`.

Changing the reference or pipeline version means reprocessing every batch, and a `v2.0`.

TPM compares genes within a sample and samples within a dataset. It is not comparable
across datasets (tissue vs cultured cells, different library preparations).
