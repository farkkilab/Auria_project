
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

# ---- SFigure 1D ----
rds_data <- readRDS("/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/TNBC180225.rds")

rds_data$oxstress_group <- ifelse(rds_data@meta.data$oxstress > 0.15, "High", "Low")

cancer_rds_data <- subset(rds_data, celltype_major == "Cancer Epithelial")

p1 <- VlnPlot(cancer_rds_data, features = c("GCLC", "NQO1", "TXNRD1"),
        group.by = "oxstress_group",
        pt.size = 0,  
        cols = c("steelblue", "firebrick"))
ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/ox_markers_violin_tumor.pdf", plot = p1, width = 6, height = 5, dpi = 300)


