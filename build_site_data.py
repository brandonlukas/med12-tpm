#!/usr/bin/env python3
"""Write the explorer's data into docs/data/ (served by GitHub Pages).

- meta.json: datasets, groups, samples (column order of every TPM array), contrasts
- genes.json: [gene_id, gene_name, gene_type] for every gene with nonzero TPM somewhere
- chunks/<i>.json: per gene, TPM per sample, expression percentile per group and the
  DESeq2 results, for genes i*CHUNK .. (i+1)*CHUNK-1 of genes.json, so the explorer
  fetches ~100 kB per gene instead of everything
- de/<dataset>__<contrast>.json: one comparison as columns, for the DE results page
- tf/<dataset>__<contrast>.json: TF activity in one comparison (tf_activity.py), as columns
- tf/network.json: each scored TF's targets, [activated, repressed], for the target lists

Run after the pipeline and deseq2.R:
    uv run --with pandas python build_site_data.py
"""
import gzip
import json
import math
import re
from decimal import ROUND_DOWN, Decimal
from pathlib import Path

import numpy as np
import pandas as pd

CHUNK = 256
GTF = "/home/brandon/data/reference_genomes/gencode_v44/gencode.v44.primary_assembly.annotation.gtf.gz"
here = Path(__file__).resolve().parent
out = here / "docs/data"
for d in ("chunks", "de", "tf"):
    (out / d).mkdir(parents=True, exist_ok=True)

# datasets.json is the one list of datasets and comparisons; the site reads it from meta.json
CONFIG = json.load(open(here / "datasets.json"))["datasets"]
DATASETS = {d["id"]: {"groups": d["groups"], "contrasts": [c["id"] for c in d["contrasts"]]} for d in CONFIG}


def sig(x):
    """Round to 3 significant digits; None for NaN so it serializes as null."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return None
    return float(f"{x:.3g}")


def cut(x):
    """3 significant digits, rounded toward zero, for p-values and fold changes. Rounding to nearest
    could carry 0.0497 to 0.05 or 0.9996 to 1 and so change a significance call; toward zero keeps a
    value on its side of any cut-off with up to 3 significant digits, including the ones a reader types."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return None
    d = Decimal(repr(float(x)))
    return float(d.quantize(Decimal(1).scaleb(d.adjusted() - 2), rounding=ROUND_DOWN))


def dump(obj, path):
    path.write_text(json.dumps(obj, separators=(",", ":")))


gene_type = {}
with gzip.open(GTF, "rt") as f:
    for line in f:
        if "\tgene\t" in line:
            attrs = line.split("\t")[8]
            gene_type[re.search(r'gene_id "([^"]+)"', attrs)[1]] = re.search(r'gene_type "([^"]+)"', attrs)[1]

samples = pd.read_csv(here / "samples.tsv", sep="\t", dtype=str)
tpm = pd.read_csv(here / "results_merged/salmon.merged.gene_tpm.tsv", sep="\t")
tpm = tpm[tpm[samples["sample"]].sum(axis=1) > 0].reset_index(drop=True)
tpm["gene_type"] = tpm["gene_id"].map(gene_type)
index = {g: i for i, g in enumerate(tpm["gene_id"])}
coding = tpm["gene_type"] == "protein_coding"
# protein-coding genes with zero TPM in every sample are not in tpm; they count as zeros in the reference
n_coding_all = sum(t == "protein_coding" for t in gene_type.values())
undetected = n_coding_all - int(coding.sum())

# Percentile of each gene's group-mean TPM among every protein-coding gene in GENCODE v44:
# "is this gene expressed well?" judged against the genes people usually target. A fixed
# reference, so adding a dataset (which makes more genes nonzero somewhere) moves nothing.
pct = {}
meta = {"chunk": CHUNK, "n_coding": n_coding_all, "datasets": []}
for ds, cfg in DATASETS.items():
    s = samples[samples["dataset"] == ds].copy()
    s["group"] = pd.Categorical(s["group"], cfg["groups"])
    s = s.sort_values(["group", "sample"])
    cfg["samples"] = s["sample"].tolist()
    for g in cfg["groups"]:
        means = tpm[s.loc[s["group"] == g, "sample"]].mean(axis=1)
        ref = np.concatenate([np.zeros(undetected), means[coding].sort_values().to_numpy()])
        pct[ds, g] = (ref.searchsorted(means.to_numpy(), side="right") / len(ref) * 100).round().astype(int)
    d = next(d for d in CONFIG if d["id"] == ds)
    meta["datasets"].append({
        "id": ds, "name": d["name"], "note": d["note"], "paper": d["paper"], "doi": d["doi"], "groups": cfg["groups"],
        "contrasts": [{k: c.get(k) for k in ("id", "label", "short", "sides", "pair")} for c in d["contrasts"]],
        # what the hover label says about each sample; empty fields are left out
        "samples": [{k: v for k, v in {"sample": r.sample, "group": r.group, "patient": r.patient,
                     "genotype": r.genotype if r.group.startswith("LM") else None,
                     "hormone": r.hormone_treatment if d["paired"] else None,
                     "also": r.same_patient_as}.items() if isinstance(v, str) and v}
                    for r in s.itertuples()],
    })
dump(meta, out / "meta.json")
dump(tpm[["gene_id", "gene_name", "gene_type"]].values.tolist(), out / "genes.json")

# per comparison: columns for the DE page, and per gene the same rounded values for the chunks
de, scored = {}, set()
for ds, cfg in DATASETS.items():
    for c in cfg["contrasts"]:
        d = pd.read_csv(here / f"deliverables/de/{ds}__{c}.csv")
        d = d[d["gene_id"].isin(index)]
        base, lfc, padj, mle = ([f(x) for x in d[col]] for f, col in
                                ((sig, "baseMean"), (cut, "log2FC_shrunk"), (cut, "padj"), (cut, "log2FoldChange")))
        site = (np.array(padj, dtype=float) < 0.05) & (np.abs(np.array(lfc, dtype=float)) >= 1)
        assert (site == ((d["padj"] < 0.05) & (d["log2FC_shrunk"].abs() >= 1)).to_numpy()).all(), f"{ds}__{c}: rounding changed a call"
        de[ds, c] = dict(zip(d["gene_id"], zip(base, lfc, padj, mle)))
        dump({
            "g": [index[g] for g in d["gene_id"]],
            "base": base,
            "lfc": lfc,
            "mle": mle,
            "p": [cut(x) for x in d["pvalue"]],
            "padj": padj,
            "z": [sig(x) for x in d["stat"]],   # Wald statistic: what TF activity is computed from
        }, out / f"de/{ds}__{c}.json")
        t = pd.read_csv(here / f"deliverables/tf/{ds}__{c}.csv")
        dump({"tf": t["tf"].tolist(), "act": [sig(x) for x in t["activity"]], "p": [cut(x) for x in t["pvalue"]],
              "padj": [cut(x) for x in t["padj"]], "n": t["n_targets"].tolist()}, out / f"tf/{ds}__{c}.json")
        scored.update(t["tf"])

net = pd.read_csv(here / json.load(open(here / "datasets.json"))["network"])
net = net[net["source"].isin(scored)]
dump({tf: [sorted(g.loc[g["weight"] > 0, "target"]), sorted(g.loc[g["weight"] < 0, "target"])]
      for tf, g in net.groupby("source")}, out / "tf/network.json")

for i in range(0, len(tpm), CHUNK):
    chunk = {}
    for j, row in enumerate(tpm.iloc[i:i + CHUNK].to_dict("records"), start=i):
        g = row["gene_id"]
        entry = {"tpm": {}, "pct": {}, "de": {}}
        for ds, cfg in DATASETS.items():
            entry["tpm"][ds] = [sig(row[s]) for s in cfg["samples"]]
            entry["pct"][ds] = [int(pct[ds, grp][j]) for grp in cfg["groups"]]
            entry["de"][ds] = {}
            for c in cfg["contrasts"]:
                if g in de[ds, c]:  # absent = removed by the expression filter
                    entry["de"][ds][c] = list(de[ds, c][g])
        chunk[g] = entry
    dump(chunk, out / f"chunks/{i // CHUNK}.json")

print(f"{len(tpm)} genes ({int(coding.sum())} protein-coding) in {math.ceil(len(tpm) / CHUNK)} chunks -> {out}")
