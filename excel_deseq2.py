#!/usr/bin/env python3
"""One Excel workbook with every DESeq2 comparison, one sheet each, plus a sheet that
explains the columns and the design. Reads deliverables/de/*.csv (from deseq2.R).

    uv run --with pandas --with openpyxl python excel_deseq2.py [out.xlsx]
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd

here = Path(__file__).resolve().parent
out = Path(sys.argv[1]) if len(sys.argv) > 1 else here / "deliverables" / "MED12_RNAseq_DESeq2.xlsx"

# one sheet per comparison, in datasets.json order; Excel allows 31 characters per sheet name
CONFIG = json.load(open(here / "datasets.json"))["datasets"]
sheet_name = lambda d, c: f"{d['sheet']} " + re.sub(r", all \d+ pairs", "", c["label"]).replace(", hormone-naive", " naive")
SHEETS = [(sheet_name(d, c), f"{d['id']}__{c['id']}", f"{d['paper']}: {c['about']}") for d in CONFIG for c in d["contrasts"]]
assert all(len(name) <= 31 for name, _, _ in SHEETS) and len({n for n, _, _ in SHEETS}) == len(SHEETS)
COLUMNS = [
    ("gene_id", "Ensembl gene ID (GENCODE v44)"),
    ("gene_name", "Gene symbol"),
    ("baseMean", "Mean normalized count over all samples in the dataset"),
    ("log2FC_shrunk", "log2 fold change, first group relative to the second, with noisy estimates from low-count genes pulled toward zero (apeglm). This is the value shown on the website."),
    ("log2FoldChange", "log2 fold change, raw estimate"),
    ("lfcSE", "Standard error of the raw log2 fold change"),
    ("pvalue", "p-value (DESeq2 Wald test)"),
    ("padj", "p-value adjusted for the number of genes tested (Benjamini-Hochberg). Empty when DESeq2 set the gene aside as too lowly expressed to adjust."),
]
order = [c for c, _ in COLUMNS]

readme = pd.DataFrame(
    [("Sheet", "Comparison", "Genes tested")] + [("", "", "")]
    + [(name, design, "") for name, _, design in SHEETS]
    + [("", "", ""), ("Column", "Meaning", "")] + [(c, m, "") for c, m in COLUMNS]
    + [("", "", ""),
       ("Notes", "Fitted separately within each dataset with DESeq2 1.50.2 on Salmon gene counts (nf-core/rnaseq 3.27.0, GENCODE v44). "
                 "A gene was tested only if it had at least 10 counts in at least as many samples as the smallest group. "
                 "Rows are sorted by p-value. Website and methods: https://brandonlukas.github.io/med12-tpm/", "")],
    columns=["A", "B", "C"])

with pd.ExcelWriter(out, engine="openpyxl") as xl:
    readme.to_excel(xl, sheet_name="README", index=False, header=False)
    tested = {}
    for name, stem, _ in SHEETS:
        df = pd.read_csv(here / "deliverables" / "de" / f"{stem}.csv")[order]
        df = df.round({"baseMean": 2, "log2FC_shrunk": 4, "log2FoldChange": 4, "lfcSE": 4})
        df.to_excel(xl, sheet_name=name, index=False)
        tested[name] = int(df["pvalue"].notna().sum())
        print(f"{name:24s} {len(df):6d} rows, {tested[name]:6d} tested")
    # fill in the tested counts on the README sheet
    ws = xl.book["README"]
    for row in ws.iter_rows(min_row=3, max_row=2 + len(SHEETS)):
        row[2].value = tested.get(row[0].value, "")
    for ws in xl.book.worksheets:
        ws.freeze_panes = "A2" if ws.title != "README" else None
        for col, width in zip("ABCDEFGH", (22, 14, 12, 14, 14, 10, 12, 12)):
            ws.column_dimensions[col].width = width if ws.title != "README" else (26, 110, 14)["ABC".index(col)] if col in "ABC" else 12
print(f"wrote {out} ({out.stat().st_size // 1_000_000} MB)")
