#!/usr/bin/env bash
# Run nf-core/rnaseq (Salmon-only) on the samples in samplesheet.csv.
# Usage: [ONLY=<regex on sample>] [OUTDIR=<dir>] ./run_rnaseq.sh [extra nextflow args]
# Salmon quantifies each sample on its own, so a later batch can run into its own
# OUTDIR and be merged with merge_salmon.R; v1.1 added the Moyo and Leistico tissue
# cohorts this way:  ONLY='^(moyo|leistico)_' OUTDIR=results_v1.1 ./run_rnaseq.sh
set -euo pipefail
cd "$(dirname "$0")"

RNASEQ_VERSION="3.27.0"   # requires Nextflow >= 25.10.4
PROFILE="${PROFILE:-docker}"
OUTDIR="${OUTDIR:-results}"
ONLY="${ONLY:-.}"

# Absolute, symlink-resolved paths: fastq/, work/ and results/ may be symlinks
# (e.g. to /mnt/wdssd), and Docker bind-mounts the path it is given, so an
# unresolved symlink can dangle inside the container.
FASTQ_DIR="$(realpath fastq)"
{ head -1 samplesheet.csv; tail -n +2 samplesheet.csv | tr -d '\r' | grep -E "$ONLY"; } \
    | sed "s#,fastq/#,$FASTQ_DIR/#g" > samplesheet.abs.csv

missing=$(tail -n +2 samplesheet.abs.csv | tr ',' '\n' | grep '\.fastq\.gz$' | while read -r f; do [[ -f "$f" ]] || echo "$f"; done)
if [[ -n "$missing" ]]; then
    echo "Missing FASTQs (run ./download_fastqs.sh first):"; echo "$missing"; exit 1
fi

nextflow run nf-core/rnaseq -r "$RNASEQ_VERSION" \
    -profile "$PROFILE" \
    -params-file params.yaml \
    -c home.config \
    -w "$(realpath work)" \
    --outdir "$(realpath -m "$OUTDIR")" \
    -resume \
    "$@"

echo
echo "Done. Now run: python3 split_tpm.py"
