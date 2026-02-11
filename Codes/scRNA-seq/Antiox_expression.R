
title: "Auria_scRNA_DEG"
author: "Wenqing Chen"
date: "2025-04-10"

# Load packages

library(ggplot2)
library(dplyr)
library(tibble)
library(readr)
library(Seurat)
library(MoMAColors)
library(ggpubr)
library(patchwork)

# ---- SFigure 2c ----
rds_data <- readRDS("/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/TNBC180225.rds")

rds_data$oxstress_group <- ifelse(rds_data@meta.data$oxstress > 0.15, "High", "Low")

cancer_rds_data <- subset(rds_data, celltype_major == "Cancer Epithelial")


genes <- c("GCLC", "NQO1", "TXNRD1")

plots <- lapply(genes, function(g){
  VlnPlot(
    cancer_rds_data,
    features = g,
    group.by = "oxstress_group",
    pt.size = 0,
    cols = c("steelblue", "firebrick")
  ) + 
    stat_compare_means(method = "wilcox.test", label = "p.signif") + # p.format
    theme(legend.position = "none") +
    ggtitle(g)
})

p_all <- wrap_plots(plots, ncol = 3)


ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/ox_markers_violin_tumor_sig.pdf", plot = p_all, width = 6, height = 5, dpi = 300)


