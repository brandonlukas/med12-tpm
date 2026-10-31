#!/usr/bin/env python3
"""Transcription-factor activity for every DESeq2 comparison: decoupleR's univariate linear
model (ULM) over the CollecTRI network, on DESeq2's Wald statistic, the recipe decoupleR
gives for bulk contrasts. Positive activity: the TF's targets move as if it were more active
in the comparison's first group.

The network is a frozen copy (networks/), not a fresh download, because OmniPath updates
CollecTRI and a rerun would otherwise give different numbers. Fetch a new copy with --fetch;
that is a new data version.

Writes deliverables/tf/<dataset>__<contrast>.csv: tf, activity, pvalue, padj (BH within the
comparison), n_targets (targets measured in that comparison).

    uv run --with decoupler==2.1.6 python tf_activity.py [--fetch]
"""
import datetime
import json
import sys
from pathlib import Path

import decoupler as dc
import pandas as pd
from statsmodels.stats.multitest import multipletests

here = Path(__file__).resolve().parent
NET = here / "networks/collectri_human.csv.gz"
out_dir = here / "deliverables/tf"
out_dir.mkdir(parents=True, exist_ok=True)

if "--fetch" in sys.argv:
    # identical under the academic and commercial licenses (checked 2026-10-31), so free to redistribute
    net = dc.op.collectri(organism="human", license="commercial")
    NET.parent.mkdir(exist_ok=True)
    net.to_csv(NET, index=False)
    print(f"fetched CollecTRI {datetime.date.today()}: {len(net):,} edges, {net['source'].nunique()} TFs -> {NET}")
net = pd.read_csv(NET)[["source", "target", "weight"]]

for d in json.load(open(here / "datasets.json"))["datasets"]:
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
        sig = df[df["padj"] < 0.05]
        print(f"{name}: {len(df)} TFs, {len(sig)} padj<0.05 ({(sig['activity'] > 0).sum()} up, {(sig['activity'] < 0).sum()} down)")
