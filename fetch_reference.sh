#!/usr/bin/env bash
# Download GENCODE v44 full (comprehensive) annotation + GRCh38 primary assembly.
# v44 is the same GENCODE release behind 10x refdata-gex-GRCh38-2024-A, so bulk
# and single-cell work share gene models; only the GTF filtering differs.
set -euo pipefail
DEST="${1:-/home/brandon/data/reference_genomes/gencode_v44}"
BASE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44
mkdir -p "$DEST"; cd "$DEST"

for f in GRCh38.primary_assembly.genome.fa.gz gencode.v44.primary_assembly.annotation.gtf.gz; do
    curl -fSL --retry 10 --retry-all-errors -C - -o "$f" "$BASE/$f"
done
curl -fsSL -o MD5SUMS "$BASE/MD5SUMS"
grep -E " (GRCh38.primary_assembly.genome.fa.gz|gencode.v44.primary_assembly.annotation.gtf.gz)$" MD5SUMS | md5sum -c -
echo "Reference ready in $DEST"
