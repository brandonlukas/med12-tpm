#!/usr/bin/env Rscript
# DESeq2 within each dataset; no comparison crosses datasets (tissue vs cells,
# different library preps and read lengths). Datasets and comparisons: datasets.json.
#   Tissue (paired): ~ patient + group, so each LM is compared with its own patient's
#     MM. LM_mut vs LM_wt is estimable because every patient has an MM baseline: it
#     asks whether the LM-vs-MM shift differs by MED12 genotype. LM vs MM pooled over
#     genotype uses ~ patient + tissue, on all pairs or on a subset ("naive_mut": the
#     patients without hormone treatment whose LM is MED12 mutant).
#   Cells: ~ group, each clone vs WT.
# Input: length-scaled counts from nf-core (tximport countsFromAbundance =
# "lengthScaledTPM"), which tximport recommends for count-matrix DESeq2.
# Output: deliverables/de/<dataset>__<contrast>.csv (log2FoldChange unshrunken,
# log2FC_shrunk from apeglm)
# Usage: Rscript deseq2.R [counts.tsv] [dataset ...]
suppressPackageStartupMessages({
    library(DESeq2)
    library(data.table)
    library(jsonlite)
})

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
setwd(here)
args <- commandArgs(TRUE)
counts_path <- if (length(args)) args[1] else "results_merged/salmon.merged.gene_counts_length_scaled.tsv"
only <- args[-1]
out_dir <- "deliverables/de"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- fread("samples.tsv")
counts <- fread(counts_path)
genes <- counts[, .(gene_id, gene_name)]
have <- intersect(samples$sample, names(counts))
mat <- round(as.matrix(counts[, have, with = FALSE]))
rownames(mat) <- genes$gene_id

samples[, tissue := ifelse(group == "MM", "MM", "LM")]
plan <- fromJSON("datasets.json", simplifyDataFrame = FALSE)$datasets

# log2FC_shrunk: apeglm, which needs the contrast as a model coefficient, so the
# factor is releveled to make the denominator the reference. p-values come from
# the unshrunken Wald test either way.
write_de <- function(dds, var, num, den, dataset, name) {
    if (levels(dds[[var]])[1] != den) {
        dds[[var]] <- relevel(dds[[var]], den)
        dds <- nbinomWaldTest(dds, quiet = TRUE)
    }
    coef <- sprintf("%s_%s_vs_%s", var, num, den)
    res <- results(dds, name = coef)
    shr <- lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE)
    dt <- data.table(gene_id = rownames(res), as.data.frame(res))
    dt[, log2FC_shrunk := shr$log2FoldChange]
    dt <- merge(genes, dt, by = "gene_id", sort = FALSE)
    dt[, `:=`(dataset = dataset, contrast = name)]
    setorder(dt, pvalue, na.last = TRUE)
    fwrite(dt, file.path(out_dir, sprintf("%s__%s.csv", dataset, name)))
    cat(sprintf("%-15s %-14s padj<0.05: %5d\n", dataset, name, sum(dt$padj < 0.05, na.rm = TRUE)))
}

fit <- function(meta, design) {
    dds <- DESeqDataSetFromMatrix(mat[, meta$sample], meta, design)
    # a gene needs >= 10 counts in at least as many samples as the smallest group
    keep <- rowSums(counts(dds) >= 10) >= min(table(meta$group))
    DESeq(dds[keep, ], quiet = TRUE)
}

subsets <- list(
    naive_mut = function(meta) {
        keep <- meta$patient[meta$group == "LM_mutMED12" & meta$hormone_treatment == "none"]
        meta[meta$patient %in% keep, ]
    }
)

for (p in plan) {
    ds <- p$id
    if (length(only) && !ds %in% only) next
    meta <- as.data.frame(samples[dataset == ds])
    rownames(meta) <- meta$sample
    meta$group <- factor(meta$group, levels = unlist(p$groups))
    if (p$paired) {
        meta$patient <- factor(meta$patient)
        meta$tissue <- factor(meta$tissue, levels = c("MM", "LM"))
    }
    dds <- NULL   # the ~ group fit, shared by every group-vs-group comparison
    for (k in p$contrasts) {
        m <- meta
        if (!is.null(k$subset)) m <- droplevels(subsets[[k$subset]](meta))
        if (k$num == "LM") {
            write_de(fit(m, ~ patient + tissue), "tissue", "LM", "MM", ds, k$id)
        } else {
            if (is.null(dds)) dds <- fit(meta, if (p$paired) ~ patient + group else ~ group)
            write_de(dds, "group", k$num, k$den, ds, k$id)
        }
    }
}
