#!/usr/bin/env python3
"""Write the gene explorer's data into docs/data/ (served by GitHub Pages).

- meta.json: datasets, groups, samples (column order of every TPM array) and contrasts
- genes.json: [gene_id, gene_name] for every gene with nonzero TPM in some sample
- chunks/<i>.json: TPM and DESeq2 results for genes i*CHUNK .. (i+1)*CHUNK-1 of
  genes.json, so the page fetches ~100 kB per gene lookup instead of everything.

Run after the pipeline and deseq2.R:
    uv run --with pandas python build_site_data.py
"""
import json
import math
from pathlib import Path

import pandas as pd

CHUNK = 256
here = Path(__file__).resolve().parent
out = here / "docs/data"
(out / "chunks").mkdir(parents=True, exist_ok=True)

DATASETS = {
    "JCI_tissue": {
        "label": "Tissue (JCI Insight 2023)",
        "groups": ["MM", "LM_wtMED12", "LM_mutMED12"],
        "contrasts": ["LM_vs_MM", "LMmut_vs_MM", "LMwt_vs_MM", "LMmut_vs_LMwt"],
    },
    "NC_2D": {
        "label": "2D cells (Nat Commun 2023)",
        "groups": ["WT", "KO", "Mut1", "Mut2"],
        "contrasts": ["KO_vs_WT", "Mut1_vs_WT", "Mut2_vs_WT"],
    },
    "NC_3D_spheroid": {
        "label": "3D spheroids (Nat Commun 2023)",
        "groups": ["WT", "Mut1", "Mut2"],
        "contrasts": ["Mut1_vs_WT", "Mut2_vs_WT"],
    },
}


def sig(x, digits=3):
    """Round to significant digits; None for NaN so it serializes as null."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return None
    return float(f"{x:.{digits}g}")


samples = pd.read_csv(here / "samples.tsv", sep="\t", dtype=str)
samples["patient"] = samples["sample"].str.extract(r"^tissue_(\d+)_")[0]
tpm = pd.read_csv(here / "results/salmon/salmon.merged.gene_tpm.tsv", sep="\t")
tpm = tpm[tpm[samples["sample"]].sum(axis=1) > 0].reset_index(drop=True)

meta = {"chunk": CHUNK, "datasets": []}
for ds, cfg in DATASETS.items():
    s = samples[samples["dataset"] == ds].copy()
    s["group"] = pd.Categorical(s["group"], cfg["groups"])
    s = s.sort_values(["group", "sample"])
    cfg["samples"] = s["sample"].tolist()
    meta["datasets"].append({
        "id": ds, "label": cfg["label"], "groups": cfg["groups"], "contrasts": cfg["contrasts"],
        "samples": [{"sample": r.sample, "group": r.group, "patient": r.patient if isinstance(r.patient, str) else None}
                    for r in s.itertuples()],
    })
(out / "meta.json").write_text(json.dumps(meta, separators=(",", ":")))
(out / "genes.json").write_text(json.dumps(tpm[["gene_id", "gene_name"]].values.tolist(), separators=(",", ":")))

# de[(dataset, contrast)] = DataFrame indexed by gene_id
de = {}
for ds, cfg in DATASETS.items():
    for c in cfg["contrasts"]:
        d = pd.read_csv(here / f"deliverables/de/{ds}__{c}.csv").set_index("gene_id")
        de[ds, c] = d[["baseMean", "log2FoldChange", "padj"]]

for i in range(0, len(tpm), CHUNK):
    block = tpm.iloc[i:i + CHUNK]
    chunk = {}
    for row in block.to_dict("records"):
        g = row["gene_id"]
        entry = {"tpm": {}, "de": {}}
        for ds, cfg in DATASETS.items():
            entry["tpm"][ds] = [sig(row[s]) for s in cfg["samples"]]
            entry["de"][ds] = {}
            for c in cfg["contrasts"]:
                d = de[ds, c]
                if g in d.index:  # absent = removed by the expression filter
                    r = d.loc[g]
                    entry["de"][ds][c] = [sig(r.baseMean), sig(r.log2FoldChange), sig(r.padj, 2)]
        chunk[g] = entry
    (out / f"chunks/{i // CHUNK}.json").write_text(json.dumps(chunk, separators=(",", ":")))

print(f"{len(tpm)} genes in {math.ceil(len(tpm) / CHUNK)} chunks -> {out}")
