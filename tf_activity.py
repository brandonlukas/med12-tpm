#!/usr/bin/env python3
"""Transcription-factor activity for every DESeq2 comparison: decoupleR's univariate linear
model (ULM) over the CollecTRI network, on DESeq2's Wald statistic, the recipe decoupleR
gives for bulk contrasts. Positive activity: the TF's targets move as if it were more active
in the comparison's first group.

The network is a frozen copy, the file datasets.json names under "network", not a fresh
download, because OmniPath updates CollecTRI and a rerun would otherwise give different numbers.
--fetch saves today's copy next to it; pointing datasets.json at it is a new data version.

Writes deliverables/tf/<dataset>__<contrast>.csv: tf, activity, pvalue, padj (BH within the
comparison), n_targets (targets measured in that comparison); and the same as one Excel workbook,
deliverables/MED12_RNAseq_TF_activity.xlsx, sheets named as in MED12_RNAseq_DESeq2.xlsx.

    uv run --with decoupler==2.1.6 --with statsmodels --with openpyxl python tf_activity.py [--fetch]
"""
import datetime
import json
import sys
from pathlib import Path

import decoupler as dc
import pandas as pd
from statsmodels.stats.multitest import multipletests

from sheet_names import check, sheet_name

here = Path(__file__).resolve().parent
cfg = json.load(open(here / "datasets.json"))
NET = here / cfg["network"]   # networks/collectri_human_<date fetched>.csv.gz
out_dir = here / "deliverables/tf"
out_dir.mkdir(parents=True, exist_ok=True)

if "--fetch" in sys.argv:
    # identical under the academic and commercial licenses (checked 2026-10-31), so free to redistribute
    net = dc.op.collectri(organism="human", license="commercial")
    new = here / f"networks/collectri_human_{datetime.date.today()}.csv.gz"
    net.to_csv(new, index=False)
    sys.exit(f"fetched CollecTRI: {len(net):,} edges, {net['source'].nunique()} TFs -> {new}; set \"network\" in datasets.json to use it")
net = pd.read_csv(NET)[["source", "target", "weight"]]

sheets = {}
for d in cfg["datasets"]:
    for c in d["contrasts"]:
        name = f'{d["id"]}__{c["id"]}'
        de = pd.read_csv(here / "deliverables/de" / f"{name}.csv").dropna(subset=["stat"])
        # a few names belong to several GENCODE genes (X/Y pseudoautosomal copies and the like):
        # keep the most expressed, since a network names genes by symbol
        de = de.sort_values("baseMean", ascending=False).drop_duplicates("gene_name")
        mat = de.set_index("gene_name")[["stat"]].T
        acts, pvals = dc.mt.ulm(data=mat, net=net, tmin=5)
        measured = net[net["target"].isin(mat.columns)].groupby("source").size()
        df = pd.DataFrame({"tf": acts.columns, "activity": acts.iloc[0].values, "pvalue": pvals.iloc[0].values})
        df = df.dropna(subset=["pvalue"])
        df["padj"] = multipletests(df["pvalue"], method="fdr_bh")[1]
        df["n_targets"] = df["tf"].map(measured).astype(int)
        df = df.sort_values("pvalue")
        df.to_csv(out_dir / f"{name}.csv", index=False, float_format="%.6g")
        sheets[sheet_name(d, c)] = df
        sig = df[df["padj"] < 0.05]
        print(f"{name}: {len(df)} TFs, {len(sig)} padj<0.05 ({(sig['activity'] > 0).sum()} up, {(sig['activity'] < 0).sum()} down)")

check(list(sheets))
readme = pd.DataFrame([
    ("Column", "Meaning"),
    ("tf", "Transcription factor (CollecTRI; AP1 and NFKB are complexes, scored from the targets of their members)"),
    ("activity", "TF activity, first group relative to the second: positive when the TF's activated targets go up and its repressed targets go down (decoupleR ULM t-value on the DESeq2 Wald statistic)"),
    ("pvalue", "p-value of the activity"),
    ("padj", "p-value adjusted for the number of TFs scored in the comparison (Benjamini-Hochberg)"),
    ("n_targets", "Targets of the TF measured in the comparison; a TF needs at least 5"),
    ("", ""),
    ("Notes", f"Network: {NET.name}, CollecTRI as downloaded on the date in its name, {len(net):,} TF-target links, frozen in the repository. "
              "In tissue, activity also reflects which cells a sample contains: myometrium has more immune and blood-vessel cells than leiomyoma. "
              "Website and methods: https://brandonlukas.github.io/med12-tpm/"),
])
with pd.ExcelWriter(here / "deliverables/MED12_RNAseq_TF_activity.xlsx", engine="openpyxl") as xl:
    readme.to_excel(xl, sheet_name="README", index=False, header=False)
    for name, df in sheets.items():
        df.to_excel(xl, sheet_name=name, index=False)
        xl.book[name].freeze_panes = "A2"
    for col, width in zip("AB", (14, 120)):
        xl.book["README"].column_dimensions[col].width = width
print(f"wrote MED12_RNAseq_TF_activity.xlsx ({len(sheets)} comparisons)")
