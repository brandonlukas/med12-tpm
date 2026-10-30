#!/usr/bin/env python3
"""One Excel workbook with every DESeq2 comparison, one sheet each, plus a sheet that
explains the columns and the design. Reads deliverables/de/*.csv (from deseq2.R).

    uv run --with pandas --with openpyxl python excel_deseq2.py [out.xlsx]
"""
import sys
from pathlib import Path

import pandas as pd

here = Path(__file__).resolve().parent
out = Path(sys.argv[1]) if len(sys.argv) > 1 else here / "deliverables" / "MED12_RNAseq_DESeq2.xlsx"

# sheet name (31 characters at most), file, design in words
SHEETS = [
    ("Tissue LM vs MM",        "JCI_tissue__LM_vs_MM",        "Leiomyoma vs matched myometrium, all 10 patients (paired by patient)"),
    ("Tissue LM mut vs MM",    "JCI_tissue__LMmut_vs_MM",     "MED12-mutant leiomyoma (6) vs myometrium (paired by patient)"),
    ("Tissue LM WT vs MM",     "JCI_tissue__LMwt_vs_MM",      "MED12 wild-type leiomyoma (4) vs myometrium (paired by patient)"),
    ("Tissue LM mut vs LM WT", "JCI_tissue__LMmut_vs_LMwt",   "MED12-mutant vs wild-type leiomyoma, each relative to its own myometrium"),
    ("2D KO vs WT",            "NC_2D__KO_vs_WT",             "MED12 knock-out vs unedited cells, 2D culture (3 vs 3)"),
    ("2D Mut1 vs WT",          "NC_2D__Mut1_vs_WT",           "MED12 G44N clone 1 vs unedited cells, 2D culture (3 vs 3)"),
    ("2D Mut2 vs WT",          "NC_2D__Mut2_vs_WT",           "MED12 G44N clone 2 vs unedited cells, 2D culture (3 vs 3)"),
    ("3D Mut1 vs WT",          "NC_3D_spheroid__Mut1_vs_WT",  "MED12 G44N clone 1 vs unedited cells, spheroids (3 vs 2)"),
    ("3D Mut2 vs WT",          "NC_3D_spheroid__Mut2_vs_WT",  "MED12 G44N clone 2 vs unedited cells, spheroids (3 vs 2)"),
]
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
