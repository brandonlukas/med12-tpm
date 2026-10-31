#!/usr/bin/env python3
"""Split nf-core/rnaseq merged Salmon gene TPM into one table per dataset.

Reads results_merged/salmon.merged.gene_tpm.tsv and samples.tsv, writes
deliverables/<dataset>_gene_TPM.csv and <dataset>_gene_counts.csv, plus one Excel workbook of TPM (one sheet per
dataset + a sample-info sheet) if openpyxl is installed.
"""
import json
import sys
from pathlib import Path

import pandas as pd

here = Path(__file__).resolve().parent
tpm_path = Path(sys.argv[1]) if len(sys.argv) > 1 else here / "results_merged/salmon.merged.gene_tpm.tsv"
out_dir = here / "deliverables"
out_dir.mkdir(exist_ok=True)

samples = pd.read_csv(here / "samples.tsv", sep="\t", dtype=str)
tpm = pd.read_csv(tpm_path, sep="\t")
# the counts DESeq2 was given: length-scaled, rounded to whole numbers as deseq2.R does
counts = pd.read_csv(tpm_path.with_name("salmon.merged.gene_counts_length_scaled.tsv"), sep="\t")
id_cols = [c for c in ("gene_id", "gene_name") if c in tpm.columns]

missing = set(samples["sample"]) - set(tpm.columns)
if missing:
    sys.exit(f"Samples missing from {tpm_path}: {sorted(missing)}")

sheets = {}
order = [d["id"] for d in json.load(open(here / "datasets.json"))["datasets"]]
for dataset in order:
    meta = samples[samples["dataset"] == dataset]
    # order columns by group then sample so replicates sit together
    cols = meta.sort_values(["group", "sample"])["sample"].tolist()
    df = tpm[id_cols + cols].copy()
    df[cols] = df[cols].round(3)
    df.to_csv(out_dir / f"{dataset}_gene_TPM.csv", index=False)
    c = counts[id_cols + cols].copy()
    c[cols] = c[cols].round().astype(int)
    c.to_csv(out_dir / f"{dataset}_gene_counts.csv", index=False)
    sheets[dataset] = df
    print(f"{dataset}: {len(cols)} samples, {len(df)} genes")

info = samples[["sample", "dataset", "group", "patient", "genotype", "genotype_as_listed_on_GEO_(incorrect)", "hormone_treatment", "gsm", "srr", "geo_title"]]
info.to_csv(out_dir / "sample_info.csv", index=False)

try:
    with pd.ExcelWriter(out_dir / "MED12_RNAseq_gene_TPM.xlsx") as xl:
        info.to_excel(xl, sheet_name="sample_info", index=False)
        for name, df in sheets.items():
            df.to_excel(xl, sheet_name=name, index=False)
    print(f"Wrote {out_dir / 'MED12_RNAseq_gene_TPM.xlsx'}")
except ImportError:
    print("openpyxl not installed; wrote CSVs only (pip install openpyxl for the .xlsx)")
