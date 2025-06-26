# NOMEDA plots exploration (oxstress signature). 
# Created: 20/01/25
# Last modification: 27/01/25
# Author(s): Sara Palomino


library(Seurat)
library(dplyr)
library(scCustomize)
library(ggplot2)
library(UCell)
library(BuenColors)
library(genekitr)
library(clusterProfiler)
library(fgsea)
library(ComplexHeatmap)
library(tidyr)
library(org.Hs.eg.db)
library(scales)
library(ComplexHeatmap)

## Load seurat object
scSeurat <- readRDS("/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/scSeurat.rds")

table(scSeurat@meta.data$Patient)
table(scSeurat@meta.data$Patient, scSeurat@meta.data$subtype)

# Subset only TNBC patients
Idents(scSeurat) <- "subtype"
TNBC  <- subset(x=scSeurat, idents ="TNBC")
table(TNBC@meta.data$Patient)

# There are two patients with no cancer cells, I would remove them
# >> CID3946  & CID44041  
Idents(TNBC) <- "Patient"
TNBC <- subset(x=TNBC, idents ="CID3946", invert=TRUE )
TNBC <- subset(x=TNBC, idents ="CID44041", invert=TRUE )

# Also remove Normal epithelial cells
Idents(TNBC) <- "celltype_major"
TNBC <- subset(x=TNBC, idents ="Normal Epithelial", invert=TRUE )

DimPlot(TNBC, group.by = "celltype_major")


## Add metadata for further exploration
clinical_metadata <- data.frame(
  patient = c("CID3963", "CID4465", "CID4495", "CID44971", "CID44991", "CID4513", "CID4515", "CID4523"),
  age = c(61, 54,	63,	49,	47,	73,	67,	52),
  cancer_type  = c("IDC",	"IDC",	"IDC",	"IDC",	"IDC",	"MBC",	"IDC",	"MBC"),
  ki67   =  c(0.43,	0.7,	0.8,	0.4,	0.65,	0.75,	0.6,	0.9),
  treatment   =  c("Treated",	"Naïve",	"Naïve",	"Naïve",	"Naïve",	"Treated",	"Naïve",	"Treated"), 
  tumour_cells = c(222 , 124 , 1184, 894  ,4018  ,   1058  ,   2169  ,   1167 ))

TNBC@meta.data$treatment <- as.character(TNBC@meta.data$Patient)
TNBC@meta.data$age <- as.character(TNBC@meta.data$Patient)
TNBC@meta.data$cancer_type <- as.character(TNBC@meta.data$Patient)

for(i in 1:length(clinical_metadata$patient)){
  patient <- clinical_metadata$patient
  treatment <- clinical_metadata$treatment
  age <- clinical_metadata$age
  cancer_type <- clinical_metadata$cancer_type
  
  TNBC@meta.data$treatment[TNBC@meta.data$Patient==patient[i]] <- treatment[i]
  TNBC@meta.data$age[TNBC@meta.data$Patient==patient[i]] <- age[i]
  TNBC@meta.data$cancer_type[TNBC@meta.data$Patient==patient[i]] <- cancer_type[i]
  
}


DimPlot(TNBC, group.by = "treatment")


# Save into a vector all genes expressed in the TNBC object
TNBC_allgenes <- rownames(TNBC)

# Load the oxidative stress signature provided by Tuulia and Nomeda
oxstress <- read.csv("/Volumes/Lola/ONEDRIVE/Projects/AURIA/00.Data/gex/OxStressSig.csv")

# Check if all the genes in the oxstress are expressed in the TNBC ("GPX6"  "GSTT1" "GSTT2" are not present!)
oxstress$Oxidative.Stress[!oxstress$Oxidative.Stress %in% TNBC_allgenes]

# Compute the signature
TNBC <- AddModuleScore_UCell(TNBC, features = list(oxstress))
TNBC@meta.data$oxstress <- TNBC@meta.data$signature_1_UCell 
TNBC@meta.data$signature_1_UCell <- NULL


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


tmp_oxstress <-  FetchData(AddModuleScore_UCell(TNBC, features = list( oxstress)), vars=c('UMAP_1', 'UMAP_2', 'Patient','celltype_major',"oxstress"))
UMAP_OXSTRESS <- get_umap_signature(tmp_oxstress, "oxstress", " ")

UMAP_celltypes <- DimPlot(TNBC, group.by = "celltype_major") +  ggtitle(NULL) + NoAxes() + 
  theme(legend.text =   element_text(size=10), legend.key.size = unit(0.5, "cm"))


##############
# Exploring expression in tumor cells
TNBC.tumor <- subset(TNBC, subset = normal_cell_call == 'cancer')
table(TNBC.tumor@meta.data$Patient)

VlnPlot(TNBC.tumor, features = "oxstress", raster=F, group.by = "Patient", alpha = 0.4) +
  xlab("Patient") + ylab("Expression") + ggtitle("Oxidative Stress (Tumour)")  +
  theme(legend.position = "none", axis.text = element_text(size = 10) , axis.line = element_line(size = 0.3),  # Thinner axis lines
        axis.title = element_text(size = 12))

TNBC.tumor <- ScaleData(TNBC.tumor, features =  rownames(TNBC.tumor))

# Heatmap of the markers
TNBC.tumor_scaled <- ScaleData(object = TNBC.tumor, features = oxstress$Oxidative.Stress)

Idents(TNBC.tumor_scaled) <- "normal_cell_call"
# Compute average expression per patient
avg_exp_tumor <- AverageExpression(TNBC.tumor_scaled, features = oxstress$Oxidative.Stress, group.by = "Patient", return.seurat = TRUE)

DoHeatmap(avg_exp_tumor, features = oxstress$Oxidative.Stress, label = FALSE) + 
  scale_fill_gradient2(low = c("#3361A5", "#248AF3", "#14B3FF", "#88CEEF"), 
                       mid = "#C1D5DC", 
                       high = c("#EAD397", "#FDB31A" ,"#E42A2A" ,"#A31D1D"), 
                       midpoint = 0, 
                       guide = "colourbar", 
                       aesthetics = "fill")

library(ComplexHeatmap)
Idents(TNBC.tumor_scaled) <- "Patient"
avg_exp_tumor <- as.data.frame(AverageExpression(TNBC.tumor_scaled))

avg_pat_oxstress <- avg_exp_tumor[rownames(avg_exp_tumor) %in% oxstress$Oxidative.Stress,]
Heatmap(avg_pat_oxstress)

Z_score <- as.data.frame( t(scale(t(avg_pat_oxstress), center = TRUE, scale = TRUE) ))
names(Z_score) <- unique(TNBC.tumor_scaled$Patient)

Heatmap(Z_score, row_names_gp = gpar(fontsize = 8), cluster_columns = F , cluster_rows= F )


###################################################
# Define the antioxidant genes
antioxidant_genes <- oxstress$Oxidative.Stress[oxstress$Oxidative.Stress %in% TNBC_allgenes]

# Get metadata to filter by patient
metadata <- TNBC.tumor@meta.data %>%
  as.data.frame()%>%
  mutate(CellID = rownames(.))  

# Sample between 100 and 1000 cells per patient
metadata <- metadata %>%
  group_by(Patient) %>%
  mutate(sample_size = min(150, max(100, n()))) %>%  # Precompute sample sizes
  ungroup()

# Initialize an empty vector to store sampled cell names
sampled_cells <- c()
# Loop over each patient
for (patient_id in unique(metadata$Patient)) {
  # Get all cell IDs for this patient
  patient_cells <- metadata$CellID[metadata$Patient == patient_id]
  # Determine how many cells to sample
  num_cells <- length(patient_cells)
  sample_size <-  unique( metadata$sample_size[metadata$Patient == patient_id])# Ensure between 100 and 1000
  # Sample the correct number of cells
  sampled_cells <- c(sampled_cells, sample(patient_cells, sample_size, replace = FALSE))
}

# Verify number of sampled cells per patient
table(metadata[metadata$CellID %in% sampled_cells, "Patient"])

# Subset the Seurat object to retain only sampled cells
TNBC.tumor_subset <- subset(TNBC.tumor, cells = sampled_cells)


Idents(TNBC.tumor_subset) <- "Patient"
avg_exp_tumor <- as.data.frame(AverageExpression(TNBC.tumor_subset))

avg_pat_oxstress <- avg_exp_tumor[rownames(avg_exp_tumor) %in% oxstress$Oxidative.Stress,]

Z_score <- as.data.frame( t(scale(t(avg_pat_oxstress), center = TRUE, scale = TRUE) ))
names(Z_score) <- unique(TNBC.tumor_scaled$Patient)

Z_score[is.nan(as.matrix(Z_score))] <- 0
ht <- Heatmap(Z_score, row_names_gp = gpar(fontsize = 8), cluster_columns = T ,cluster_rows= T,  na_col = "grey" )

# Draw heatmap and retrieve column order
ht_drawn <- draw(ht)  # Draw the heatmap first
column_ordered <- column_order(ht_drawn)[[1]]  # Extract ordered column indices

# Get ordered column names that will be used afterwards
ordered_colnames <- unique(colnames(mat))[column_ordered]


Idents(TNBC.tumor_subset) <- "celltype_major"
DotPlot(TNBC.tumor_subset, features = antioxidant_genes)


# heatmap showing subsampled single cells from each patient and annotating the patients

TNBC.tumor_subset <- ScaleData(object = TNBC.tumor_subset, features = antioxidant_genes)
Idents(TNBC.tumor_subset) <- "Patient"

DoHeatmap(TNBC.tumor_subset, features = antioxidant_genes, label = FALSE) + 
  scale_fill_gradient2(low = c("#3361A5", "#248AF3", "#14B3FF", "#88CEEF"), 
                       mid = "#C1D5DC", 
                       high = c("#EAD397", "#FDB31A" ,"#E42A2A" ,"#A31D1D"), 
                       midpoint = 0, 
                       guide = "colourbar", 
                       aesthetics = "fill")

library(ComplexHeatmap)


mat <- as.data.frame(TNBC.tumor_subset@assays$RNA@layers$scale.data)
rownames(mat) <- antioxidant_genes

patients <- data.frame(Patient =TNBC.tumor_subset$Patient)
names(mat) <- patients$Patient


patient_colors <- c("#1f77b4", "#d62728","#2ca02c",  "#ff7f0e", "#9467bd", "#636363","#bcbd22","#17becf")
names(patient_colors) <- (unique(patient_ids_factor))

column_ha = HeatmapAnnotation(Patient = (names(mat) ),
                              col = list(Patient = patient_colors)  )

Heatmap(as.matrix(mat), top_annotation = column_ha, 
        cluster_rows = TRUE, 
        cluster_columns = FALSE, 
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        column_title = "Oxstress genes in tumour cells (downsampled Nmax=150)",
        row_names_gp = gpar(fontsize = 6),
        name="Expression", show_column_names = FALSE  )



# "#1f77b4" "#d62728" "#2ca02c" "#ff7f0e" "#9467bd" "#636363" "#bcbd22" "#17becf" 



# Step 1: Extract patient ID from column names
patient_ids <- names(mat) 
# Step 2: Convert to factor using the provided order
patient_ids_factor <- factor(patient_ids, levels = ordered_patients)

# Step 3: Reorder columns so that cells from the same patient are grouped together in the ordered sequence
mat_ordered <- mat[, order(patient_ids_factor)]
cols <- sub("\\..*", "", colnames(mat_ordered))

names(mat_ordered) <- cols

patient_colors <- c("#1f77b4", "#d62728","#2ca02c",  "#ff7f0e", "#9467bd", "#636363","#bcbd22","#17becf")
names(patient_colors) <- (unique(patient_ids_factor))
  
column_ha = HeatmapAnnotation(Patient = (cols),
                              col = list(Patient = patient_colors)  )

Heatmap(as.matrix(mat_ordered), top_annotation = column_ha, 
        cluster_rows = TRUE, 
        cluster_columns = FALSE, 
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        column_title = "Oxstress genes in tumour cells (downsampled Nmax=150)",
        row_names_gp = gpar(fontsize = 6),
        name="Expression", show_column_names = FALSE  )


library(pheatmap)

TNBC.tumor_scaled <- ScaleData(object = TNBC.tumor, features = oxstress$Oxidative.Stress)
# Extract expression data for antioxidant genes
expr_matrix <- GetAssayData(TNBC.tumor_scaled, assay = "RNA", layer = "data")[antioxidant_genes, ]

# Remove genes that are not expressed (all zeros)
expr_matrix <- expr_matrix[rowSums(expr_matrix) > 0, ]

# Compute correlation matrix
cor_matrix <- cor(t(as.matrix(expr_matrix)), method = "pearson")  # Can also use "spearman"


prueba <- (as.data.frame(expr_matrix))
# Plot correlation matrix using pheatmap
pheatmap(cor_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         breaks = seq(-1, 1, length.out = 51),  # Properly maps the color scale to correlations
         display_numbers = F, 
         main = "Antioxidant Gene Correlation in Tumor Cells")

pheatmap(cor_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Smooth transition
         breaks = seq(-1, 1, length.out = 101),  # Ensure color mapping is uniform
         display_numbers = FALSE,  # Remove numbers to keep it clean
         border_color = NA,  # Removes the grid lines
         main = "Antioxidant Gene Correlation in Tumor Cells",
         fontsize = 10,  # Reduce font size for better readability
         fontsize_row = 8,  # Adjust row label font size
         fontsize_col = 8,  # Adjust column label font size
         legend_breaks = c(-1, -0.5, 0, 0.5, 1),  # Improve legend readability
         legend_labels = c("-1 (Neg)", "-0.5", "0 (None)", "0.5", "1 (Pos)")
)


# Identify Highly Expressed Genes in a Signature
# Get the average expression for each gene
avg_exp <- rowMeans(GetAssayData(TNBC, assay = "RNA", slot = "data")[antioxidant_genes, ])

# Rank genes by expression
top_genes <- sort(avg_exp, decreasing = TRUE)

# Print top highly expressed genes
head(top_genes, 10)  # Show top 10

setwd("/Volumes/Lola/ONEDRIVE/Projects/AURIA/02.Results/")
pdf("oxstress_Dotplot.pdf", height = 14, width = 5)
oxstress_Dotplot <- DotPlot(TNBC.tumor, antioxidant_genes) +
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +  RotatedAxis()+ xlab('') +  ylab('') + 
  scale_size(breaks = c(0,25,50, 75, 100), limits=c(0,100), name = "Percentage") + labs(title=" ") + 
  theme(axis.text = element_text(size = 8) ,legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_text(size = 8)) +
  coord_flip()
oxstress_Dotplot
dev.off()                                     

                                                            

