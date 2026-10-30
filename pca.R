#!/usr/bin/env Rscript
# PCA per dataset exactly as DESeq2's plotPCA does it: vst-transformed counts, the 500 most
# variable genes, prcomp without scaling, variance explained from sdev^2. Genes with fewer than
# 10 counts in total are dropped first, as the DESeq2 vignette recommends. Writes
# docs/data/pca.json: per dataset, the samples' PC1 and PC2 scores and the variance explained.
suppressPackageStartupMessages({ library(DESeq2); library(data.table); library(jsonlite) })
setwd(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
samples <- fread("samples.tsv"); counts <- fread("results/salmon/salmon.merged.gene_counts_length_scaled.tsv")
out <- list()
for (ds in c("JCI_tissue", "NC_2D", "NC_3D_spheroid")) {
  meta <- as.data.frame(samples[dataset == ds]); rownames(meta) <- meta$sample
  mat <- round(as.matrix(counts[, meta$sample, with = FALSE])); rownames(mat) <- counts$gene_id
  dds <- DESeqDataSetFromMatrix(mat[rowSums(mat) >= 10, ], meta, ~ 1)
  pd <- plotPCA(vst(dds, blind = TRUE), intgroup = "group", ntop = 500, returnData = TRUE)
  ve <- round(100 * attr(pd, "percentVar"), 1)
  pd <- pd[meta$sample, ]
  out[[ds]] <- list(explained = ve, samples = lapply(seq_len(nrow(meta)), function(i) list(
    sample = meta$sample[i], group = meta$group[i], pc1 = round(pd$PC1[i], 2), pc2 = round(pd$PC2[i], 2))))
  cat(sprintf("%-15s PC1 %4.1f%%  PC2 %4.1f%%\n", ds, ve[1], ve[2]))
}
write_json(out, "docs/data/pca.json", auto_unbox = TRUE, pretty = FALSE)
