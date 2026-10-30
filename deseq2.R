#!/usr/bin/env Rscript
# DESeq2 within each dataset; no comparison crosses datasets (tissue vs cells,
# different library preps and read lengths).
#   JCI_tissue: ~ patient + group, so each LM is compared with its own patient's
#     MM. LM_mut vs LM_wt is estimable here because every patient has an MM
#     baseline: it asks whether the LM-vs-MM shift differs by MED12 genotype.
#     LM vs MM pooled over genotype uses ~ patient + tissue.
#   NC_2D, NC_3D_spheroid: ~ group, each clone vs WT.
# Input: length-scaled counts from nf-core (tximport countsFromAbundance =
# "lengthScaledTPM"), which tximport recommends for count-matrix DESeq2.
# Output: deliverables/de/<dataset>__<contrast>.csv
# Usage: Rscript deseq2.R [counts.tsv]
suppressPackageStartupMessages({
    library(DESeq2)
    library(data.table)
})

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
setwd(here)
args <- commandArgs(TRUE)
counts_path <- if (length(args)) args[1] else "results/salmon/salmon.merged.gene_counts_length_scaled.tsv"
out_dir <- "deliverables/de"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- fread("samples.tsv")
counts <- fread(counts_path)
genes <- counts[, .(gene_id, gene_name)]
mat <- round(as.matrix(counts[, samples$sample, with = FALSE]))
rownames(mat) <- genes$gene_id

samples[, patient := ifelse(dataset == "JCI_tissue", sub("^tissue_(\\d+)_.*", "\\1", sample), NA)]
samples[, tissue := ifelse(group == "MM", "MM", "LM")]

# dataset -> list(levels of group, list(contrast name = c(numerator, denominator)))
plan <- list(
    JCI_tissue = list(
        levels = c("MM", "LM_wtMED12", "LM_mutMED12"),
        contrasts = list(
            LMmut_vs_MM   = c("LM_mutMED12", "MM"),
            LMwt_vs_MM    = c("LM_wtMED12", "MM"),
            LMmut_vs_LMwt = c("LM_mutMED12", "LM_wtMED12")
        )
    ),
    NC_2D = list(
        levels = c("WT", "KO", "Mut1", "Mut2"),
        contrasts = list(KO_vs_WT = c("KO", "WT"), Mut1_vs_WT = c("Mut1", "WT"), Mut2_vs_WT = c("Mut2", "WT"))
    ),
    NC_3D_spheroid = list(
        levels = c("WT", "Mut1", "Mut2"),
        contrasts = list(Mut1_vs_WT = c("Mut1", "WT"), Mut2_vs_WT = c("Mut2", "WT"))
    )
)

write_de <- function(dds, contrast, dataset, name) {
    res <- results(dds, contrast = contrast)
    dt <- data.table(gene_id = rownames(res), as.data.frame(res))
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

for (ds in names(plan)) {
    p <- plan[[ds]]
    meta <- as.data.frame(samples[dataset == ds])
    rownames(meta) <- meta$sample
    meta$group <- factor(meta$group, levels = p$levels)
    if (ds == "JCI_tissue") {
        meta$patient <- factor(meta$patient)
        meta$tissue <- factor(meta$tissue, levels = c("MM", "LM"))
        write_de(fit(meta, ~ patient + tissue), c("tissue", "LM", "MM"), ds, "LM_vs_MM")
        dds <- fit(meta, ~ patient + group)
    } else {
        dds <- fit(meta, ~ group)
    }
    for (nm in names(p$contrasts)) write_de(dds, c("group", p$contrasts[[nm]]), ds, nm)
}
