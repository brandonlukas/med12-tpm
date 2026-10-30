# MED12 RNA-seq explorer

Three public MED12 RNA-seq datasets reprocessed from raw reads the same way, with an
interactive site: **https://brandonlukas.github.io/med12-tpm/**

- **Overview**: a gene search, an illustration of where the samples come from, and why the studies belong together
- **Genes**: TPM of any gene (or several side by side) in every sample, with DESeq2 results
- **Comparisons**: volcano plot, ranked table and gene lookup for every DESeq2 comparison,
  and agreement between any two comparisons
- **Methods**: how everything was processed, with downloads

| Dataset | Paper | GEO | Samples | Reads |
|---|---|---|---|---|
| Tissue | Zuberi et al. 2023 JCI Insight (PMID 37607000) | GSE204969 | 10 patients × MM/LM = 20 | SE 75 bp |
| 2D cells | Buyukcelebi et al. 2023 Nat Commun | GSE226014 | WT, KO, Mut1, Mut2 × 3 = 12 | PE 2×75 |
| 3D spheroids | same | GSE226016 | WT ×2, Mut1 ×3, Mut2 ×3 = 8 | PE 2×151 |

**Genotype correction:** GEO has the MED12 genotypes of patients 1007 and 1008 swapped.
1008 LM is the G44D mutant and 1007 LM is wild-type. `samples.tsv` uses the corrected
genotype and keeps GEO's label in `genotype_as_listed_on_GEO_(incorrect)`.

Gene TPM and DESeq2 tables are attached to the
[v1.0 release](https://github.com/brandonlukas/med12-tpm/releases/tag/v1.0).
The [methods page](https://brandonlukas.github.io/med12-tpm/methods.html) documents every step.

## Files

- `samples.tsv`: sample metadata (corrected genotype, GSM/SRR IDs, expected read counts)
- `samplesheet.csv`: nf-core/rnaseq samplesheet (strandedness auto-detected)
- `params.yaml`: pipeline parameters (Salmon-only)
- `home.config`: CPU/RAM limits; edit for your machine
- `fetch_reference.sh`: GENCODE v44 fasta + GTF, md5-checked
- `download_fastqs.sh`: 40 runs from NCBI SRA's AWS mirror, read counts checked
- `run_rnaseq.sh`: launches nf-core/rnaseq 3.27.0
- `split_tpm.py`: per-dataset TPM tables and an Excel workbook
- `deseq2.R`: DESeq2 within each dataset
- `build_site_data.py`: data for the site (`docs/data/`)
- `illustration/vector.py`: generates the overview's illustration (smooth vector SVG, gradient shading) into `docs/index.html`
- `stamp.py`: stamps the pages' stylesheet and script links with a content hash; run it before committing a change to `docs/`
- `tests/site_test.py`: scripted browser checks of the site (Playwright); `tests/screenshots.py` renders pages to PNG
- `docs/`: the site, served by GitHub Pages: `index.html` + `overview.js` (Overview), `genes.html` + `app.js`
  (Genes), `de.html` + `de.js` (Comparisons), `methods.html`, `common.js`, `style.css`; Plotly 2.35.2 is vendored in `docs/vendor/`
  so the site has no third-party dependency at load time

## Reproducing (Linux / WSL2)

Needs Nextflow ≥ 25.10.4, Docker, sra-tools, pigz, R with DESeq2 and data.table, and
Python with pandas and openpyxl; about 250 GB of free disk.

```bash
./fetch_reference.sh            # ~1 GB
./download_fastqs.sh 4          # ~90 GB, about 1 h
./run_rnaseq.sh                 # about 2 h on 44 cores, including the Salmon index
python3 split_tpm.py            # deliverables/*_gene_TPM.csv, .xlsx
Rscript deseq2.R                # deliverables/de/*.csv
python3 build_site_data.py      # docs/data/
```

TPM compares genes within a sample and samples within a dataset. It is not comparable
across the three datasets (tissue vs cultured cells, different library preparations).
