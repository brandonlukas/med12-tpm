#!/usr/bin/env Rscript
# PCA per dataset on variance-stabilized counts (DESeq2 vst, top 500 most variable genes),
# the standard sample-level check. Writes docs/data/pca.json: per dataset, the samples' PC1
# and PC2 scores and the variance each PC explains.
suppressPackageStartupMessages({ library(DESeq2); library(data.table); library(jsonlite) })
setwd(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
samples <- fread("samples.tsv"); counts <- fread("results/salmon/salmon.merged.gene_counts_length_scaled.tsv")
out <- list()
for (ds in c("JCI_tissue", "NC_2D", "NC_3D_spheroid")) {
  meta <- as.data.frame(samples[dataset == ds]); rownames(meta) <- meta$sample
  mat <- round(as.matrix(counts[, meta$sample, with = FALSE])); rownames(mat) <- counts$gene_id
  dds <- DESeqDataSetFromMatrix(mat[rowSums(mat) >= 10, ], meta, ~ 1)
  v <- assay(vst(dds, blind = TRUE))
  top <- head(order(rowVars(v), decreasing = TRUE), 500)
  p <- prcomp(t(v[top, ]), scale. = FALSE)
  ve <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
  out[[ds]] <- list(explained = ve[1:4], samples = lapply(seq_len(nrow(meta)), function(i) list(
    sample = meta$sample[i], group = meta$group[i], pc1 = round(p$x[i, 1], 2), pc2 = round(p$x[i, 2], 2))))
  cat(sprintf("%-15s PC1 %4.1f%%  PC2 %4.1f%%  PC3 %4.1f%%\n", ds, ve[1], ve[2], ve[3]))
  print(data.frame(group = meta$group, PC1 = round(p$x[, 1], 1), PC2 = round(p$x[, 2], 1)), row.names = TRUE)
}
write_json(out, "docs/data/pca.json", auto_unbox = TRUE, pretty = FALSE)
