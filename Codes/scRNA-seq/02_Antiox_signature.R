### Antiox signature
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

# ---- Color list ----
patient_colors <- c(
  # ER+
  "CID4461" = "#f7dfb8", "CID4463" = "#ffda5f", "CID4471" = "#ffba51",
  "CID4530N" = "#fd8302", "CID4535" = "#a26605", "CID3941" = "#f2a4a4",
  "CID3948" = "#ff8989", "CID4067" = "#EF3B2C", "CID4290A" = "#A50F15",
  
  # TNBC
  "CID4495" = "#d4ebd1", "CID44971" = "#A1D99B","CID44991" = "#41AB5D",
  "CID4513" = "#9ECAE1", "CID4515" = "#6BAED6", "CID4523" = "#3182BD",
  "CID3963" = "#08519C"
)

subtype_color <- c("ER+" = "#BFE4EE", "HER2+" = "#FBD178", "TNBC" = "#F9C3BF")

celltype_major_color <- c(
  "Endothelial" = "#abdce7", "CAFs" = "#b6c3f1", "PVL" = "#a3d7bf",
  "B-cells" = "#ddd69d", "T-cells" = "#bcd59b", "Myeloid" = "#dbbce2",
  "Plasmablasts" = "#99b6d7", "Cancer Epithelial" = "#f1cdce" )


# ---- ER+ & TNBC file ----
## ---- Load data ----
### ---- ER+ & TNBC ----
scSeurat <- readRDS("/Volumes/Lola/ONEDRIVE/Projects/AURIA/00.Data/gex/scSeurat_def.rds")

### ----  Antiox signature ---- 
oxstress <- read.csv("/Volumes/Lola/ONEDRIVE/Projects/AURIA/00.Data/gex/OxStressSig.csv")

# Compute the signature
scSeurat <- AddModuleScore_UCell(scSeurat, features = list(oxstress))
scSeurat@meta.data$oxstress <- scSeurat@meta.data$signature_1_UCell 
scSeurat@meta.data$signature_1_UCell <- NULL


table(scSeurat@meta.data$Patient, scSeurat@meta.data$subtype)
table(scSeurat@meta.data$Patient, scSeurat@meta.data$celltype_major)

# Calculate the mean oxstress
df_oxstress <- scSeurat@meta.data %>% select(Patient, celltype_major, subtype, oxstress)

df_oxstress <- df_oxstress %>%
  group_by(Patient) %>%
  mutate(mean_oxstress_patient = mean(oxstress, na.rm = TRUE)) %>%
  ungroup()

df_oxstress <- df_oxstress %>%
  group_by(Patient, celltype_major) %>%
  mutate(mean_oxstress_patient_celltype = mean(oxstress, na.rm = TRUE)) %>%
  ungroup()

summary_df <- df_oxstress %>%
  distinct(Patient, celltype_major, subtype, mean_oxstress_patient, mean_oxstress_patient_celltype)

# because here we only need tumor cells, so we should use mean_oxstress_patient_celltype, not mean_oxstress_patient
summary_df_tumor <- subset(summary_df, (celltype_major == "Cancer Epithelial"))
summary_df_tumor <- summary_df_tumor %>% distinct(Patient, subtype, mean_oxstress_patient_celltype)

# Only TNBC patients [per celltype per patient plot: only TNBC]
summary_df_TNBC <- subset(summary_df, subset = subtype == "TNBC")


### ---- TNBC data ----
# recompute the signature when only TNBC
# Subset only TNBC patients
Idents(sc_seurat) <- "subtype"
scSeurat_TNBC  <- subset(x=sc_seurat, idents ="TNBC")
table(scSeurat_TNBC@meta.data$Patient)

scSeurat_TNBC$oxstress <- NULL

scSeurat_TNBC <- AddModuleScore_UCell(scSeurat_TNBC, features = list(oxstress))
scSeurat_TNBC@meta.data$oxstress <- scSeurat_TNBC@meta.data$oxstress_tumor 
scSeurat_TNBC@meta.data$signature_1_UCell <- NULL

scSeurat_TNBC@meta.data %>%
  group_by(Patient) %>%
  summarise(mean_oxstress = mean(oxstress, na.rm = TRUE))


scSeurat_TNBC <- FetchData(AddModuleScore_UCell(scSeurat_TNBC, features = list( oxstress)), vars=c('UMAP_1', 'UMAP_2', 'Patient','celltype_major',"oxstress"))
tmp <-  FetchData(AddModuleScore_UCell(scSeurat_TNBC, features = list( oxstress)), vars=c('UMAP_1', 'UMAP_2', 'Patient','celltype_major',"oxstress"))

### Plot it in a UMAP
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
get_umap_signature(tmp, "oxstress", " ")


saveRDS(scSeurat_TNBC, "scSeurat_TNBC_def.rds")


