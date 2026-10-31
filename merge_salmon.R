#!/usr/bin/env Rscript
# Put the gene tables of the nf-core/rnaseq batches side by side: the batches are
# samples.tsv's batch column, their folders in datasets.json. Each batch keeps its own tables: TPM is
# per sample, but length-scaled counts depend on the samples tximport saw, so
# re-running tximport over everything would shift v1.0's released DE results.
# Every dataset lies within one run, and DESeq2 never crosses datasets.
#
# Usage: Rscript merge_salmon.R   (writes results_merged/)
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
setwd(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))

out <- "results_merged"
samples <- fread("samples.tsv")
dirs <- unname(unlist(fromJSON("datasets.json")$batches[unique(samples$batch)]))   # unnamed, or cbind prefixes the column names
dir.create(out, showWarnings = FALSE)

for (f in c("gene_tpm", "gene_counts", "gene_counts_length_scaled")) {
    tabs <- lapply(dirs, function(d) fread(file.path(d, "salmon", paste0("salmon.merged.", f, ".tsv"))))
    # one index and one GTF for every run, so the same genes in the same order
    for (t in tabs[-1]) stopifnot(identical(t$gene_id, tabs[[1]]$gene_id), identical(t$gene_name, tabs[[1]]$gene_name))
    merged <- do.call(cbind, c(tabs[1], lapply(tabs[-1], function(t) t[, -(1:2)])))
    stopifnot(!anyDuplicated(names(merged)), setequal(setdiff(names(merged), c("gene_id", "gene_name")), samples$sample))
    fwrite(merged, file.path(out, paste0("salmon.merged.", f, ".tsv")), sep = "\t")
    message(f, ": ", ncol(merged) - 2, " samples")
}
