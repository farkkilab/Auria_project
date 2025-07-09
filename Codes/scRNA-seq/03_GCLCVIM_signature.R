### GCLCVIM signature
# Author(s): Sara Palomino, Wenqing Chen

# ---- Load libraries ----

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(rstatix)
library(UCell)
library(BuenColors)

### ---- TNBC ----
TNBC <- saveRDS(scSeurat_TNBC, "scSeurat_TNBC_def.rds")

### ----  GCLCVIM signature ---- 
genes <- c("GCLC", "VIM")

# Function to compute the median of non-zero values
compute_median <- function(marker_exp) {
  median_value <- median(marker_exp)
  return(median_value)
}

# Initialize lists to store gene expression and median thresholds
gene_exp_median <- list()
thresholds_median <- list()

# Define genes (Ensure you specify exactly two genes in 'genes')
for (i in 1:length(genes)){  # Only loop for two genes
  exp <- GetAssayData(TNBC, assay = "RNA", slot = "data")[genes[i], ]
  gene_exp_median[[i]] <- exp
  thresholds_median[[i]] <- compute_median(exp)
}

# Plot histograms for the two genes
hist(gene_exp_median[[1]], main = paste("Histogram of", genes[1]), xlab = genes[1])
hist(gene_exp_median[[2]], main = paste("Histogram of", genes[2]), xlab = genes[2])

# Print computed median thresholds
print(thresholds_median)

# Classify cells based on median thresholds for two genes
categories_median <- ifelse(gene_exp_median[[1]] > thresholds_median[[1]] & gene_exp_median[[2]] > thresholds_median[[2]], "++", 
                            ifelse(gene_exp_median[[1]] > thresholds_median[[1]] & gene_exp_median[[2]] <= thresholds_median[[2]], "+-",  
                                   ifelse(gene_exp_median[[1]] <= thresholds_median[[1]] & gene_exp_median[[2]] > thresholds_median[[2]], "-+",  
                                          "--")))


# Add categories as metadata to the Seurat object
TNBC$GCLCVIM <- as.factor(categories_median)

table(TNBC$GCLCVIM, TNBC$Patient)

color_palette <-  c("#1f77b4", "#bcbd22","#d62728", "#ff7f0e","#2ca02c",  "#636363","#","#17becf",  "#8c564b","#ad494a","#9467bd","#e377c2")

DimPlot(TNBC, group.by = "GCLCVIM", order = T, cols = color_palette, pt.size = 0.3, raster = F)


# Differential expression analysis to get the signature
TNBC_tumor <- subset(TNBC, subset = normal_cell_call == 'cancer')

Idents(TNBC_tumor) <- "GCLCVIM"

# GCLC+VIM+ vs all in tumor
contrast.GCLCVIM <- FindMarkers(TNBC_tumor, ident.1 = "++",  logfc.threshold=log2(1.25), min.pct = 0.75)
contrast.GCLCVIM$genes <- rownames(contrast.GCLCVIM)
contrast.GCLCVIM$delta <- contrast.GCLCVIM$pct.1-contrast.GCLCVIM$pct.2
contrast.GCLCVIM_list <- rownames(contrast.GCLCVIM)[(abs(contrast.GCLCVIM$delta)) >= 0.49 & contrast.GCLCVIM$avg_log2FC>0 ]
signature2 <- contrast.4_list

Idents(TNBC) <- "normal_cell_call"
contrast.tumor <- FindMarkers(TNBC, ident.1 = "cancer", logfc.threshold=log2(1.5), min.pct = 0.75)
contrast.tumor$genes <- rownames(contrast.tumor)
contrast.tumor$delta <- contrast.tumor$pct.1-contrast.tumor$pct.2

contrast.tumor_list <- contrast.tumor$genes[(abs(contrast.tumor$delta)) >= 0.49 & contrast.tumor$avg_log2FC>0]


signature_GCLCVIMTUMOR <- c(signature2, contrast.tumor_list)


TNBC <- AddModuleScore_UCell(TNBC, features = list(signature_GCLCVIMTUMOR))
TNBC@meta.data$signature_GCLCVIMTUMOR <- TNBC@meta.data$signature_1_UCell
TNBC@meta.data$signature_1_UCell <- NULL


TNBC_tumor <- subset(TNBC, subset = normal_cell_call == 'cancer')
Idents(TNBC_tumor) <- "GCLCVIM"


VlnPlot(TNBC, features = 'signature_GCLCVIMTUMOR', group.by = "Patient", alpha=0.3)
VlnPlot(TNBC_tumor, features = 'signature_GCLCVIMTUMOR', group.by = "Patient", alpha=0.3)

# UMAP of GCLC_VIM
get_umap_signature <- function(data, color_by_pos, title){
  data <-  data[order(data[[color_by_pos]]),]
  p <- ggplot(data, aes(x=UMAP_1, y= UMAP_2, color= .data[[color_by_pos]])) + 
    geom_point(size=0.6, alpha=0.6) +
    scale_color_gradientn(colors = jdb_palette("solar_extra")) + ggtitle(title) + 
    theme_classic() + theme(legend.position = "right", legend.key.size = unit(0.5, "cm"),panel.background = element_blank(),
                            axis.line = element_blank(),       # Remove axis lines
                            axis.ticks = element_blank(),      # Remove axis ticks
                            axis.text = element_blank(),       # Remove axis text
                            axis.title = element_blank()       # Remove axis titles # Smaller axis labels
                            
    )
  return(p)
} 

tmp <-  FetchData(AddModuleScore_UCell(TNBC, features = list( signature_GCLCVIMTUMOR)), vars=c('UMAP_1', 'UMAP_2', 'Patient','celltype_major',"signature_GCLCVIMTUMOR"))
p <- get_umap_signature(tmp2, "signature_GCLCVIMTUMOR", " ")


saveRDS(scSeurat_TNBC, "scSeurat_TNBC_def.rds")
saveRDS(signature_GCLCVIMTUMOR, "GCLCVIM_TUMOR.rds")