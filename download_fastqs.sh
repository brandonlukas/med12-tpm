#!/usr/bin/env bash
# Download runs from NCBI's AWS mirror (sra-pub-run-odp) and convert to FASTQ.
# ENA served ~0.1 MB/s per file here, the AWS mirror ~30 MB/s. Each run is
# verified by vdb-validate and by its read count against samples.tsv.
# Writes fastq/<name>.fastq.gz (single-end) or fastq/<name>_{1,2}.fastq.gz
# (paired-end), the paths in samplesheet.csv; <name> is the sample, or
# <sample>_<run> for a sample sequenced in several runs. Finished runs are skipped.
# Usage: ./download_fastqs.sh [parallel_jobs] [dataset ...]   (default 4 jobs, all datasets)
set -euo pipefail
cd "$(dirname "$0")"
JOBS="${1:-4}"
DATASETS="${*:2}"
FASTQ="$(realpath fastq)"
TMP="$FASTQ/.sra_tmp"
mkdir -p "$TMP"

fetch() {
    sample="$1"; srr="$2"; layout="$3"; reads="$4"
    if [[ "$layout" == PAIRED ]]; then outs=("${sample}_1" "${sample}_2"); else outs=("$sample"); fi
    if [[ -f "$FASTQ/${outs[-1]}.fastq.gz" ]]; then echo "ok (cached)  $sample"; return 0; fi
    d="$TMP/$srr"; rm -rf "$d"; mkdir -p "$d"
    curl -fsSL --retry 10 --retry-all-errors --retry-delay 15 \
        -o "$d/$srr.sra" "https://sra-pub-run-odp.s3.amazonaws.com/sra/$srr/$srr"
    vdb-validate "$d/$srr.sra" > "$d/validate.log" 2>&1 || { echo "VALIDATE FAILED $sample" >&2; return 1; }
    fasterq-dump --split-3 --threads 6 --temp "$d" --outdir "$d" "$d/$srr.sra" > /dev/null
    if [[ "$layout" == PAIRED ]]; then srcs=("${srr}_1" "${srr}_2"); else srcs=("$srr"); fi
    for i in "${!srcs[@]}"; do
        n=$(( $(wc -l < "$d/${srcs[$i]}.fastq") / 4 ))
        if [[ "$n" != "$reads" ]]; then
            echo "READ COUNT MISMATCH $sample ${srcs[$i]}: $n, samples.tsv says $reads" >&2; return 1
        fi
    done
    # write .tmp then rename: the last output's existence marks the run as done
    for i in "${!srcs[@]}"; do
        pigz -p 6 -c "$d/${srcs[$i]}.fastq" > "$FASTQ/${outs[$i]}.fastq.gz.tmp"
        mv "$FASTQ/${outs[$i]}.fastq.gz.tmp" "$FASTQ/${outs[$i]}.fastq.gz"
    done
    rm -rf "$d"
    echo "ok           $sample ($n reads)"
}
export -f fetch
export FASTQ TMP

# samples.tsv: one row per sample; srr and read_count list a sample's runs, comma-separated
runs() {
    python3 - "$DATASETS" <<'PY'
import csv, sys
want = sys.argv[1].split()
for r in csv.DictReader(open("samples.tsv"), delimiter="\t"):
    if want and r["dataset"] not in want:
        continue
    srrs, counts = r["srr"].split(","), r["read_count"].split(",")
    for srr, n in zip(srrs, counts):
        print(r["sample"] if len(srrs) == 1 else f"{r['sample']}_{srr}", srr, r["layout"], n)
PY
}
runs | xargs -P "$JOBS" -L 1 bash -c 'fetch "$0" "$1" "$2" "$3"'

rmdir "$TMP"
echo "All $(runs | wc -l) runs downloaded and verified."
