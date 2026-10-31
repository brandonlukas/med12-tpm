#!/usr/bin/env python3
"""Which tissue samples come from the same person, from inherited variants in the reads.

Every myometrium (MM) sample of the tissue datasets is aligned with bwa mem (8 million reads;
bwa rather than minimap2 because two studies have 42-base reads) to the
3,000 most-expressed myometrium transcripts, one per gene; variants are called jointly with
bcftools, and each pair of samples is scored by how often their genotypes agree at sites
where both are well covered and at least one carries the variant. Two samples from one
person agree at nearly every site; two people disagree at a large share. Myometrium is
used because it carries only inherited variants, not the tumor's.

    uv run --with pandas python checks/same_patient.py   # needs the varcheck conda env (bwa, bcftools, samtools)
Writes identity/ on wdssd and checks/same_patient.tsv (every pair, sorted by agreement).
"""
import csv, glob, gzip, itertools, os, subprocess, sys
from pathlib import Path

import pandas as pd

here = Path(__file__).resolve().parent.parent
work = Path("/mnt/wdssd/pingyin/med12-tpm/identity")
env = "/home/brandon/miniconda3/envs/varcheck/bin"
N_READS = 8_000_000

samples = [r for r in csv.DictReader(open(here / "samples.tsv"), delimiter="\t")
           if r["dataset"].endswith("_tissue") and r["group"] == "MM"]

def run(cmd):
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")

# reference: per gene its most expressed transcript, for the 3,000 genes highest in MM
ref = work / "ref.fa"
if not ref.exists():
    tx = pd.concat([pd.read_csv(here / d / "salmon/salmon.merged.transcript_tpm.tsv", sep="\t", index_col=0)
                    for d in ("results", "results_v1.1")], axis=1)
    tx = tx.loc[:, ~tx.columns.duplicated()]
    gene = tx.pop("gene_id")
    mm = [s["sample"] for s in samples]
    mean = tx[mm].mean(axis=1)
    best = mean.groupby(gene).idxmax()
    keep = set(mean[best.values].sort_values(ascending=False).head(3000).index)
    with open(here / "results/genome/genome.transcripts.fa") as f, open(ref, "w") as out:
        on = False
        for line in f:
            if line.startswith(">"):
                on = line[1:].split()[0].split("|")[0] in keep
            if on: out.write(line)
    run(f"{env}/samtools faidx {ref} && {env}/bwa index {ref} 2>/dev/null")
    print(f"reference: {len(keep)} transcripts")

def reads(s):
    """The first N_READS reads (read 1 only) of a sample, as a gzipped FASTQ path or a command."""
    if s["dataset"] == "JCI_tissue":
        return f"zcat {work}/sra/{s['sample']}.fastq.gz"
    first = sorted(glob.glob(str(here / "fastq" / f"{s['sample']}*_1.fastq.gz")))
    return f"(zcat {' '.join(first)} || true) | head -n {4 * N_READS}"   # head closing early is not an error

bams = []
for s in samples:
    bam = work / f"{s['sample']}.bam"
    if not bam.exists():
        run(f"set -o pipefail; {reads(s)} | {env}/bwa mem -t 16 {ref} - 2>/dev/null"
            f" | {env}/samtools sort -@4 -o {bam}.tmp - && mv {bam}.tmp {bam} && {env}/samtools index {bam}")
        print("aligned", s["sample"], flush=True)
    bams.append(str(bam))

vcf = work / "mm.vcf.gz"
if not vcf.exists():
    run(f"{env}/bcftools mpileup --threads 16 -f {ref} -a AD,DP -d 500 -q 20 -Q 20 {' '.join(bams)} 2>/dev/null"
        f" | {env}/bcftools call -mv -Oz -o {vcf} && {env}/bcftools index {vcf}")

# genotypes at biallelic SNPs with good quality
q = subprocess.run(f"{env}/bcftools view -m2 -M2 -v snps -i 'QUAL>50' {vcf} | {env}/bcftools query -f '[%GT:%DP\\t]\\n'",
                   shell=True, capture_output=True, text=True, check=True).stdout
names = [Path(b).stem for b in bams]
G = []
for line in q.splitlines():
    row = []
    for f in line.rstrip("\t").split("\t"):
        gt, dp = f.split(":")
        row.append(None if "." in gt or dp == "." or int(dp) < 10 else gt.replace("|", "/").count("1"))
    G.append(row)
pairs = []
for i, j in itertools.combinations(range(len(names)), 2):
    sites = [(r[i], r[j]) for r in G if r[i] is not None and r[j] is not None and (r[i] or r[j])]
    agree = sum(a == b for a, b in sites)
    pairs.append((names[i], names[j], len(sites), agree / len(sites) if sites else float("nan")))
df = pd.DataFrame(pairs, columns=["sample_a", "sample_b", "sites", "agreement"]).sort_values("agreement", ascending=False)
df.to_csv(here / "checks/same_patient.tsv", sep="\t", index=False, float_format="%.3f")
print(f"{len(G)} variant sites; {len(df)} pairs")
print(df.head(10).to_string(index=False))
print("...\nmedian agreement of all pairs:", round(df["agreement"].median(), 3))
