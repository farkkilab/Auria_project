
# Created: 03/07/25
# Last modification: 03/07/25
# Author(s): Wenqing Chen

library(Hmisc)
library(Rmisc)

library(corrplot)
library(ggcorrplot)
library(RColorBrewer)
library(grDevices)
library(dplyr)


# ---- Load data ----
data = read.csv('/Users/wenqchen/Desktop/Projects/Auria/Data/t-CycIF/Data_AllMarkers_filtered_RCNs.csv')
data_tumor = data %>% filter(celltype == "Tumor")


target_markers <- c("Ki67_2", "GCLC", "NQO1", "TXNRD1", "GLUT1", "H2ax", "pRB", "pSTAT3", "Vimentin", "pS6_235", "H3K27me3", "E.Cadherin")
filtered_data <- data_tumor[, target_markers, drop = FALSE]
colnames(filtered_data)[colnames(filtered_data) == "Ki67_2"] <- "Ki67"


rcorr_result <- rcorr(as.matrix(filtered_data), type = "spearman")

cor_matrix <- rcorr_result$r
p_matrix <- rcorr_result$P 

pdf("/Users/wenqchen/Desktop/Projects/Auria/Plots/Correlation/More_markers/tumor_correlation_plot.pdf", width = 10, height = 8) 
corrplot(cor_matrix, 
         method = "circle",
         type = "upper",
         col = c("#313772", "#2c4ca0", "#326db6", "#478ecc", "#75b5dc",  "#d2e2ef", "#fee3ce", "#eabaa1", "#dc917b", "#d16d5b", "#c44438", "#b7282e"),
         tl.col = "black", 
         tl.cex = 1
)
dev.off()

