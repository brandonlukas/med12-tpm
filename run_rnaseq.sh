#!/usr/bin/env bash
# Run nf-core/rnaseq (Salmon-only) on all 40 samples.
# Usage: ./run_rnaseq.sh [extra nextflow args]   e.g. ./run_rnaseq.sh -profile singularity
set -euo pipefail
cd "$(dirname "$0")"

RNASEQ_VERSION="3.27.0"   # requires Nextflow >= 25.10.4
PROFILE="${PROFILE:-docker}"

# Absolute, symlink-resolved paths: fastq/, work/ and results/ may be symlinks
# (e.g. to /mnt/wdssd), and Docker bind-mounts the path it is given, so an
# unresolved symlink can dangle inside the container.
FASTQ_DIR="$(realpath fastq)"
sed "s#,fastq/#,$FASTQ_DIR/#g" samplesheet.csv > samplesheet.abs.csv

missing=$(tail -n +2 samplesheet.abs.csv | tr ',' '\n' | grep '\.fastq\.gz$' | while read -r f; do [[ -f "$f" ]] || echo "$f"; done)
if [[ -n "$missing" ]]; then
    echo "Missing FASTQs (run ./download_fastqs.sh first):"; echo "$missing"; exit 1
fi

nextflow run nf-core/rnaseq -r "$RNASEQ_VERSION" \
    -profile "$PROFILE" \
    -params-file params.yaml \
    -c home.config \
    -w "$(realpath work)" \
    --outdir "$(realpath results)" \
    -resume \
    "$@"

echo
echo "Done. Now run: python3 split_tpm.py"
