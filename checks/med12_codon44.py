#!/usr/bin/env python3
"""MED12 codon 44 read by read, for every leiomyoma (LM) sample whose FASTQs are on disk.

Codon 44 is GGT (glycine) in the reference, chrX:71,119,403-405 (GRCh38, canonical
transcript). A read counts when it spans the codon with 10 matching bases on each side,
in either orientation. G44D reads GAT, G44S AGT, G44C TGT, G44R CGT, G44V GTT. Mutations
elsewhere in exon 2 (L36R, small deletions) leave codon 44 wild-type or unread.

    python3 checks/med12_codon44.py     # writes checks/med12_codon44.tsv
"""
import collections, csv, glob, gzip, re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

here = Path(__file__).resolve().parent.parent
FW = re.compile(r"TGTAAAACAA(...)TTCAATAACC")
RV = re.compile(r"GGTTATTGAA(...)TTGTTTTACA")
RC = str.maketrans("ACGT", "TGCA")
AA = {"GGT": "G", "GAT": "D", "AGT": "S", "TGT": "C", "CGT": "R", "GTT": "V"}

def count(sample):
    c = collections.Counter()
    for fq in sorted(glob.glob(str(here / "fastq" / f"{sample}*.fastq.gz"))):
        with gzip.open(fq, "rt") as f:
            for i, line in enumerate(f):
                if i % 4 != 1: continue
                m = FW.search(line)
                if m: c[m.group(1)] += 1; continue
                m = RV.search(line)
                if m: c[m.group(1)[::-1].translate(RC)] += 1
    return sample, c

rows = [r for r in csv.DictReader(open(here / "samples.tsv"), delimiter="\t")
        if r["group"].startswith("LM") and glob.glob(str(here / "fastq" / f"{r['sample']}*.fastq.gz"))]
with ProcessPoolExecutor(12) as ex:
    found = dict(ex.map(count, [r["sample"] for r in rows]))
with open(here / "checks/med12_codon44.tsv", "w") as out:
    out.write("sample\tdataset\tgenotype_in_paper\treads_spanning\treads_by_codon\treads_say\n")
    for r in rows:
        c = found[r["sample"]]
        n = sum(c.values())
        mut = [(k, v) for k, v in c.most_common() if k != "GGT"]
        call = f"G44{AA.get(mut[0][0], '?')}" if mut and mut[0][1] >= 3 and mut[0][1] >= 0.3 * n else ("wild type at codon 44" if c["GGT"] >= 3 else "too few reads")
        by = ", ".join(f"{k} ({'G44' + AA[k] if k in AA and k != 'GGT' else 'wild type' if k == 'GGT' else 'other'}) {v}" for k, v in c.most_common())
        out.write(f"{r['sample']}\t{r['dataset']}\t{r['genotype']}\t{n}\t{by}\t{call}\n")
print(open(here / "checks/med12_codon44.tsv").read())
