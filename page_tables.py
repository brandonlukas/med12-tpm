#!/usr/bin/env python3
"""Write the parts of the pages that list datasets, from datasets.json, samples.tsv,
Salmon's per-sample output and CHANGELOG.md, so the pages cannot drift from the data:

- methods.html: the data table, the differential-expression design table, the mapping table
- downloads.html: the per-dataset CSV table, the workbook descriptions, the version list
- every page's footer: data version, GEO accessions and processing date

Each part goes between <!-- name:start --> and <!-- name:end --> markers.
Run after the pipeline and deseq2.R:  python3 page_tables.py
"""
import csv
import json
import re
from pathlib import Path

here = Path(__file__).resolve().parent
cfg = json.load(open(here / "datasets.json"))
DS, RELEASE = cfg["datasets"], cfg["release"]
samples = list(csv.DictReader(open(here / "samples.tsv"), delimiter="\t"))
URL = f"https://github.com/brandonlukas/med12-tpm/releases/download/{RELEASE}"
FORMAT = {"SR": "reverse (SR)", "ISR": "reverse (ISR)", "SF": "forward (SF)", "ISF": "forward (ISF)", "U": "unstranded (U)", "IU": "unstranded (IU)"}


def put(page, name, html):
    path = here / "docs" / page
    s = path.read_text()
    a, b = s.index(f"<!-- {name}:start -->"), s.index(f"<!-- {name}:end -->")
    path.write_text(s[:a] + f"<!-- {name}:start -->\n" + html + "\n" + s[b:])


def geo_links(d):
    return ", ".join(f'<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}">{g}</a>' for g in d["geo"])


def mapping(d):
    """Range of Salmon mapping rates and the library formats over the dataset's samples."""
    rates, formats = [], set()
    for s in (s for s in samples if s["dataset"] == d["id"]):
        q = here / cfg["batches"][s["batch"]] / "salmon" / s["sample"]
        rates.append(json.load(open(q / "aux_info/meta_info.json"))["percent_mapped"])
        formats.add(json.load(open(q / "lib_format_counts.json"))["expected_format"])
    return f"{min(rates):.1f}–{max(rates):.1f}%", ", ".join(FORMAT.get(f, f) for f in sorted(formats))


def design_rows(d):
    rows, by_design = [], {}
    for c in d["contrasts"]:
        if c["num"] == "LM":
            design = "<code>~ patient + tissue</code>"
            if c.get("subset") == "naive_mut":
                design += ", patients without hormone treatment whose LM is <i>MED12</i> mutant"
        elif d["paired"]:
            design = "<code>~ patient + group</code>, group = " + ", ".join(GROUP[g] for g in d["groups"])
        else:
            design = "<code>~ group</code>"
        by_design.setdefault(design, []).append(c["label"])
    for design, labels in by_design.items():
        rows.append(f'  <tr><td>{d["name"]}</td><td>{design}</td><td>{"; ".join(labels)}</td></tr>')
    return rows


GROUP = {"MM": "MM", "LM_wtMED12": "LM WT", "LM_mutMED12": "LM mut"}
n = len(samples)

put("methods.html", "data", '<div class="scroll"><table>\n  <tr><th>Dataset</th><th>Paper</th><th>GEO</th><th>Samples</th><th>Reads</th></tr>\n' + "\n".join(
    f'  <tr><td>{d["name"]}</td><td><a href="https://doi.org/{d["doi"]}">{d["paper"]}</a></td><td>{geo_links(d)}</td>\n'
    f'      <td>{d["methods"]}</td><td>{d["reads"]}</td></tr>' for d in DS) + "\n</table></div>")
put("methods.html", "design", '<div class="scroll"><table>\n  <tr><th>Dataset</th><th>Design</th><th>Comparisons</th></tr>\n'
    + "\n".join(r for d in DS for r in design_rows(d)) + "\n</table></div>")
put("methods.html", "mapping", '<div class="scroll"><table>\n  <tr><th>Dataset</th><th>Mapping rate</th><th>Library strandedness (Salmon)</th></tr>\n'
    + "\n".join(f'  <tr><td>{d["name"]}</td><td>{m[0]}</td><td>{m[1]}</td></tr>' for d in DS for m in [mapping(d)]) + "\n</table></div>")
put("methods.html", "count", f"{n}")

link = lambda f, text: f'<a href="{URL}/{f}">{text}</a>'
put("downloads.html", "tables", '<div class="scroll"><table class="prose-table">\n  <tr><th>Dataset</th><th>Gene TPM</th><th>Comparisons, one table each</th></tr>\n' + "\n".join(
    f'  <tr><td>{d["name"]}</td><td>{link(d["id"] + "_gene_TPM.csv.gz", "CSV")}</td>\n      <td>'
    + " ·\n          ".join(link(f'{d["id"]}__{c["id"]}.csv.gz', c["label"]) for c in d["contrasts"]) + "</td></tr>" for d in DS)
    + f'\n</table></div>\n<p>Sample table with GSM and SRR accessions, mutation and hormone treatment: {link("sample_info.csv", "sample_info.csv")}.</p>')
n_contrasts = sum(len(d["contrasts"]) for d in DS)
put("downloads.html", "workbooks", f'''<div class="dl-grid">
  <a class="dl" href="{URL}/MED12_RNAseq_gene_TPM.xlsx">
    <span class="dl-title">Gene TPM <span class="dim">Excel</span></span>
    <span>One sheet per dataset ({len(DS)}), one row per gene, one column per sample, plus a sample sheet with the
      GEO and SRA accessions, the mutations and hormone treatment.</span>
  </a>
  <a class="dl" href="{URL}/MED12_RNAseq_DESeq2.xlsx">
    <span class="dl-title">Differential expression <span class="dim">Excel</span></span>
    <span>One sheet per comparison ({n_contrasts}), every tested gene sorted by p-value, with a README sheet that explains
      the columns and each comparison's design.</span>
  </a>
</div>''')

# versions: each "## vX.Y · date" section of CHANGELOG.md, its first paragraph
log = (here / "CHANGELOG.md").read_text()
versions = re.findall(r"^## (v[\d.]+) · (\S+)\n\n(.+?)(?:\n\n|\Z)", log, re.M | re.S)
put("downloads.html", "versions", "<ul class=\"cols\">\n" + "\n".join(
    f'  <li><b>{v}</b> <span class="dim">{date}</span>: {re.sub(r"`([^`]+)`", r"<code>\1</code>", text.replace(chr(10), " "))}</li>'
    for v, date, text in versions) + "\n</ul>")

geo = ", ".join(g for d in DS for g in d["geo"])
date = versions[0][1]
footer = f"<span>Data {RELEASE} ({date}) · GEO {geo}</span>"
for page in ("index.html", "genes.html", "de.html", "methods.html", "downloads.html"):
    put(page, "footer", footer)
print(f"wrote dataset tables for {len(DS)} datasets, {n} samples, {n_contrasts} comparisons, release {RELEASE}")
