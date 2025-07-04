
# Created: 30/06/25
# Last modification: 30/06/25
# Author(s): Wenqing Chen

# ---- Load libraries ----

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(rstatix)


scSeurat <- readRDS("/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/scSeurat.rds")
rds_data <- readRDS("/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/TNBC180225.rds")

scSeurat_TNBC <- subset(scSeurat, subset = subtype == "TNBC" )

table(rds_data@meta.data$Patient)

table(scSeurat_TNBC@meta.data$Patient)
