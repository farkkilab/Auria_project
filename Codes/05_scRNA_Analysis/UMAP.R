# Created: 02/07/25
# Last modification: 02/07/25
# Author(s): Wenqing Chen

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(rstatix)
library(UCell)
library(BuenColors)


# ---- Figure 3G ----
# UMAP of GCLC_VIM

gclcvim_tumor <- readRDS("/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/GCLCVIM_TUMOR.rds")
TNBC <- readRDS("/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/TNBC180225.rds")

TNBC <- subset(x=TNBC, subset = Patient != "CID4465")
table(TNBC@meta.data$Patient)


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

tmp2 <-  FetchData(AddModuleScore_UCell(TNBC, features = list( gclcvim_tumor)), vars=c('UMAP_1', 'UMAP_2', 'Patient','celltype_major',"signature_GCLCVIMTUMOR"))

p <- get_umap_signature(tmp2, "signature_GCLCVIMTUMOR", " ")
ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/GCLC_VIM_TNBC_umap.pdf",
       width =5.5, height = 3.5, dpi = 300)



# ---- Figure 1D ----

TNBC <- readRDS("/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/TNBC180225.rds")

TNBC <- subset(x=TNBC, subset = Patient != "CID4465")
table(TNBC@meta.data$Patient)

TNBC.tumor <- subset(TNBC, subset = celltype_major == "Cancer Epithelial")

patient.list <- SplitObject(TNBC.tumor, split.by = "Patient")  # Replace "Patient" with the metadata column name

patient.list <- lapply(patient.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
  return(x)
})

anchors <- FindIntegrationAnchors(object.list = patient.list, dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)


DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:10)
integrated <- FindNeighbors(integrated, dims = 1:10)
integrated <- FindClusters(integrated, resolution = 0.2)

DimPlot(integrated, reduction = "umap", group.by = "Patient")

DimPlot(integrated, reduction = "umap", group.by = "GCLCVIM",  order = T, pt.size = 0.3)
DimPlot(integrated, reduction = "umap", group.by = "cancer_type", order = T, pt.size = 0.3)
DimPlot(integrated, reduction = "umap", group.by = "Patient", order = T, pt.size = 0.3)

DimPlot(integrated, reduction = "umap", group.by = "treatment", order = T, pt.size = 0.3)

DimPlot(integrated, reduction = "umap", split.by  = "Patient", order = T, pt.size = 0.3)


DimPlot(integrated, reduction = "umap")
ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/TNBC_umap.pdf",
       width =5, height = 4, dpi = 300)

# Extract UMAP coordinates and metadata
umap_data <- as.data.frame(integrated@reductions$umap@cell.embeddings)
umap_data$Patient <- integrated@meta.data$Patient
umap_data$treatment <- integrated@meta.data$treatment

umap_data$GCLCVIM_Signature <- integrated@meta.data$signature_GCLCVIMTUMOR

umap_data$Oxstress_Signature <- integrated@meta.data$oxstress


# UMAP Plot, faceted by Patient
# GCLCVIM_Signature
ggplot(umap_data, aes(x = umap_1, y = umap_2, color = GCLCVIM_Signature)) +
  geom_point(size = 0.5, alpha = 0.7) +  # Adjust point size and transparency
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +  # Color scale
  theme_minimal() +
  ggtitle("GLCLC Signature Expression by Patient") +
  facet_wrap(~Patient, scales = "free")  # Facet per patient

ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/GCLC_VIM_TNBC_umap_patients.pdf",
       width =10, height = 7, dpi = 300)


# Oxstress_Signature
ggplot(umap_data, aes(x = umap_1, y = umap_2, color = Oxstress_Signature)) +
  geom_point(size = 0.5, alpha = 0.7) +  # Adjust point size and transparency
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +  # Color scale
  theme_minimal() +
  ggtitle("oxstress Signature Expression by Patient") +
  facet_wrap(~Patient, scales = "free")  # Facet per patient

ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/oxstress_TNBC_umap_patients.pdf",
       width =10, height = 7, dpi = 300)

