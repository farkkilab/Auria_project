
### scRNASeq preparation AURIA project
# Author(s): Sara Palomino, Wenqing Chen


## https://singlecell.broadinstitute.org/single_cell/study/SCP1039/a-single-cell-and-spatially-resolved-atlas-of-human-breast-cancers#study-summary

# ---- Load libraries ----
library(Seurat)


#### Prepare and explore data  #### 

expression_matrix <- ReadMtx(
  mtx = "gene_sorted-matrix.mtx", features = "features.tsv",
  cells = "barcodes.tsv.gz"
)

seurat_object <- CreateSeuratObject(counts = expression_matrix)

seurat_object = UpdateSeuratObject(object = seurat_object)

# see the sample names
table(Idents(seurat_object))

# read in meta (contains cell types)
metadata <- read.csv("Whole_miniatlas_meta.csv", header = T)
metadata <- metadata[-1,]
rownames(metadata) <- metadata$NAME

# add meta to the seurat object
seurat_object <- AddMetaData(seurat_object, meta, col.name = NULL)


# add umap coordinates
# convert to matrix
UMAP_coordinates <- read.csv("Whole_miniatlas_umap.coords.tsv", sep="\t")
UMAP_coordinates <- UMAP_coordinates[-1,]
rownames(UMAP_coordinates) <- UMAP_coordinates$NAME
colnames(UMAP_coordinates) <- c("NAME", "UMAP_1", "UMAP_2")
UMAP_coordinates$UMAP_1 <- as.numeric(UMAP_coordinates$UMAP_1)
UMAP_coordinates$UMAP_2 <- as.numeric(UMAP_coordinates$UMAP_2)
UMAP_coordinates <- UMAP_coordinates[,-1]
UMAP_coordinates_mat <- as(UMAP_coordinates, "matrix")

# Create DimReducObject and add to object
seurat_object[['UMAP']] <- CreateDimReducObject(embeddings = UMAP_coordinates_mat, key = "UMAP_", global = T, assay = "RNA")


# set the  layer 
count. <- GetAssay(object = seurat_object[["RNA"]], layer = "counts")
seurat_object <- SetAssay(
  object = seurat_object,
  layer = "",
  new. = count.,
  assay = "RNA"
)


# Subset patient type
Idents(seurat_object) <- "subtype"
seurat_object  <- subset(x=seurat_object, idents ="HER2+", invert=TRUE) # we keep TNBC and ER+
table(seurat_object$subtype)

table(seurat_object@meta.data$Patient, seurat_object@meta.data$celltype_major)

# There are patients with few cancer cells, remove them (<150)
# >> CID3946, CID4040, CID4398, & CID44041  

Idents(seurat_object) <- "Patient"
seurat_object <- subset(x=seurat_object, idents ="CID3946", invert=TRUE )
seurat_object <- subset(x=seurat_object, idents ="CID44041", invert=TRUE )
seurat_object <- subset(x=seurat_object, idents ="CID4040", invert=TRUE )
seurat_object <- subset(x=seurat_object, idents ="CID4398", invert=TRUE )
seurat_object <- subset(x=seurat_object, idents ="CID4465", invert=TRUE )

table(seurat_object$Patient)
length(table(seurat_object$Patient))

# Also remove Normal epithelial cells
Idents(seurat_object) <- "celltype_major"
seurat_object <- subset(x=seurat_object, idents ="Normal Epithelial", invert=TRUE )

# Rescale and Renormalize
DefaultAssay(seurat_object) <- "RNA"
seurat_object <- NormalizeData(seurat_object)
seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))

seurat_object <- FindVariableFeatures(object = seurat_object)
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
DimPlot(seurat_object, reduction = "pca")
ElbowPlot(seurat_object)

seurat_object <- FindNeighbors(seurat_object, dims = 1:15)
seurat_object <- RunUMAP(seurat_object, dims = 1:15)


saveRDS("scSeurat_def.rds", seurat_object)
