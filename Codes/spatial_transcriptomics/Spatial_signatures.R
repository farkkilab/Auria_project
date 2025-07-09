
### spatial transcriptomics preparation AURIA project
# Author(s): Sara Palomino, Wenqing Chen


## https://singlecell.broadinstitute.org/single_cell/study/SCP1039/a-single-cell-and-spatially-resolved-atlas-of-human-breast-cancers#study-summary
# 01: SETUP---------------------------------------
# ---- Load libraries ----
library(ggplot2)
library(patchwork)
library(dplyr)
library(Seurat)
library(Matrix)
library(hdf5r)
library(arrow)
library(zeallot)
library(STutility)
library(magrittr)
library(magick)
library(RColorBrewer)
library(ggpubr)
library(BuenColors)




# 02: PREPARE INPUT DATA --------------------------------------------------
setwd("/Volumes/Lola/ONEDRIVE/Projects/AURIA/00.Data/spatial")
infiles_h5 <- list.files("/Volumes/Lola/ONEDRIVE/Projects/AURIA/00.Data/spatial/filtered_count_matrices/", full.names = T, recursive = T)
infiles <- list.files("/Volumes/Lola/ONEDRIVE/Projects/AURIA/00.Data/spatial/spatial/", full.names = T, recursive = T)
metadata_sample <- list.files("/Volumes/Lola/ONEDRIVE/Projects/AURIA/00.Data/spatial/metadata/", full.names = T, recursive = T)
temp_csv <- read.csv("Meta_Data.csv")

infiles_h5 <- list.files("D:/ONEDRIVE/Projects/AURIA/00.Data/spatial/filtered_count_matrices/", full.names = T, recursive = T)
infiles <- list.files("D:/ONEDRIVE/Projects/AURIA/00.Data/spatial/spatial/", full.names = T, recursive = T)
metadata_sample <- list.files("D:/ONEDRIVE/Projects/AURIA/00.Data/spatial/metadata/", full.names = T, recursive = T)
temp_csv <- read.csv("Meta_Data.csv")



# We are going to focus only in the TNBC samples.
TNBCsamples <- temp_csv$Clinical_Case[temp_csv$subtype=="TNBC"]

# We need as input .h5 files (they should be in the data folder, but they are not, so let's create them, also seurat objects.)

for(sample in TNBCsamples){
  print(sample)
  #dir <- paste0("/Volumes/Lola/ONEDRIVE/Projects/AURIA/00.Data/spatial/filtered_count_matrices/", sample, sep="", "_filtered_count_matrix")
  dir <- paste0("D:/ONEDRIVE/Projects/AURIA/00.Data/spatial/filtered_count_matrices/", sample, sep="", "_filtered_count_matrix")
  
  counts <- Read10X(dir, gene.column = 1,
                    cell.column = 1,
                    unique.features = TRUE,
                    strip.suffix = FALSE)
  #setwd(paste0("/Volumes/Lola/ONEDRIVE/Projects/AURIA/00.Data/spatial/filtered_count_matrices/", sample, sep="", "_filtered_count_matrix"))
  #write10xCounts(paste0(sample, sep="_","matrix.h5"), counts, type = "HDF5")
  # ===
  data = Seurat::CreateSeuratObject(
    counts = counts , 
    project = 'test', 
    assay = 'Spatial')
  
  data$slice = 1 
  data$region = 'test' 
  
  #imgpath = paste0("/Volumes/Lola/ONEDRIVE/Projects/AURIA/00.Data/spatial/spatial/", sample,sep="_", "spatial")
  imgpath = paste0("D:/ONEDRIVE/Projects/AURIA/00.Data/spatial/spatial/", sample,sep="_", "spatial")
  
  img = Seurat::Read10X_Image(image.dir = imgpath)  
  
  Seurat::DefaultAssay(object = img) <- 'Spatial'  
  
  img = img[colnames(x = data)]  
  data[['image']] = img  
  
  #setwd(paste0("/Volumes/Lola/ONEDRIVE/Projects/AURIA/00.Data/spatial/filtered_count_matrices/", sample, sep="", "_filtered_count_matrix"))
  #saveRDS(data, paste0(sample,"_seurat.rds", sep=""))
}

# We need a reference table with all the information and directories from each sample:
infoTable <- NULL
for(sample in TNBCsamples){
  print(sample)
  infiles_sample <- grep(sample,infiles,value = T)
  infiles_sample_h5 <- grep(sample,infiles_h5,value = T)
  
  temp_df <- data.frame(samples=as.character(grep("h5",infiles_sample_h5,value = T)),
                        imgs=as.character(grep("tissue_hires_image",infiles_sample,value = T)),
                        spotfiles=as.character(grep("tissue_positions",infiles_sample,value = T)),
                        json=as.character(grep("json",infiles_sample,value = T)),
                        patientid=sample,
                        subtype=temp_csv[temp_csv$Clinical_Case==sample,"subtype"])
  infoTable <- rbind(infoTable, temp_df)
}

for(col in colnames(infoTable)){
  infoTable[,col] <- as.character(infoTable[,col])
}

infoTable <- infoTable %>% distinct()

# Also the metadata, in this case, we have the annotatios already, great!
infoMetadata <- NULL
for(sample in TNBCsamples){
  print(sample)
  
  metadata = data.frame(metadata=as.character(grep(sample,metadata_sample,value = T)))
  
  infoMetadata <- rbind(infoMetadata, metadata)
}

infoTable <- cbind(infoTable,infoMetadata)

write.csv(infoTable, "infoTable.csv")

# Metadata
# List all files in the directory
metadata_sample <- list.files("/Volumes/Lola/ONEDRIVE/Projects/AURIA/00.Data/spatial/metadata/", full.names = T, recursive = T)

# Initialize an empty list
metadata.list <- list()
# Loop through the files and read/append
for (file in metadata_sample) {
  print(paste("Reading file:", file))
  
  # Read the file (use read.csv or read.delim based on your file type)
  metadata <- read.csv(file, stringsAsFactors = FALSE)
  
  # Combine into the master data frame
  metadata.list[[file]] <- metadata
}

names(metadata.list) <- temp_csv$Clinical_Case

# Keep only TNBC patients
metadata.list <- metadata.list[intersect(names(metadata.list), TNBCsamples)]


# 03: LOAD AND PROCESS DATA ---------------------------------------------------------------
# Create a list with all the patients
se.list <- lapply(unique(infoTable$patientid), function(s) {
  print(s)
  se <- InputFromTable(infoTable[infoTable$patientid==s,,drop=F], platform = "Visium")
  se[["Spatial"]] <- se[["RNA"]]
  # Switch the default assay to Spatial
  DefaultAssay(se) <- "Spatial"
  se <- LoadImages(se,time.resolve = FALSE) 
  
})  
names(se.list) <- TNBCsamples

# Add the image information to each patient
for(sample in TNBCsamples){
  print(sample)
  imgpath = paste0("/Volumes/Lola/ONEDRIVE/Projects/AURIA/00.Data/spatial/spatial/", sample,sep="_", "spatial")
  img = Seurat::Read10X_Image(image.dir = imgpath)  
  
  Seurat::DefaultAssay(object = img) <- 'Spatial'  
  
  img = img[colnames(x = se.list[[sample]])]  
  se.list[[sample]][['image']] = img  
} 


# Now we have to add the metadata to the Seurat object for each patient and in the correct order!
for(sample in TNBCsamples){
  if (names(se.list[sample]) == names(metadata.list[sample])) {
    print(paste0(sample,": First check: sample matches metadata"))
  }
  
  meta.se <- data.frame(X=sub("_.*$", "", (colnames(se.list[[sample]]))), patientid=sample)
  metadata <- merge(meta.se, metadata.list[[sample]], by=c("X", "patientid"))
  
  metadata <-  metadata[order(match(paste(metadata$X, metadata$patientid), paste(meta.se$X, meta.se$patientid))), ]
  if( identical(metadata$X, meta.se$X)) {
    print(paste0(sample,": Second check >> ready to add the metadata!"))
  }
  
  se.list[[sample]] <- AddMetaData(se.list[[sample]], metadata = metadata, col.name = "Classification")
}



#### >>"CID4465"  ####
sample_CID4465 <- "CID4465"

se_CID4465 <- se.list[[sample_CID4465]]

se_CID4465 <- SCTransform(se_CID4465) %>%
  LoadImages(time.resolve = FALSE) 

se_CID4465 <- RunPCA(se_CID4465, assay = "SCT", verbose = FALSE)
se_CID4465 <- FindNeighbors(se_CID4465, reduction = "pca", dims = 1:30)
se_CID4465 <- FindClusters(se_CID4465, verbose = FALSE)
se_CID4465 <- RunUMAP(se_CID4465, reduction = "pca", dims = 1:30)

# Plot annotations
FeatureOverlay(se_CID4465,
               features = "Classification", 
               # pt.size = 0.75,
               # pt.alpha = 0.75,
               sampleids = 1,
               type = "raw")


DimPlot(se_CID4465, group.by = "Classification")

SpatialDimPlot(se_CID4465, group.by = "Classification")

## Add coordinates
head(se_CID4465@images$image)

# The positions were not properly added, so I'm doing it again:
spotfile <- "/Volumes/Lola/ONEDRIVE/Projects/AURIA/00.Data/spatial/spatial/CID4465_spatial/tissue_positions_list.csv"
spots <- read.csv(spotfile, header = FALSE)
# Rename columns for clarity (adjust based on your file format)
colnames(spots) <- c("barcode", "tissue", "row", "col", "imagerow", "imagecol")
# Inspect the first few rows
head(spots)

spots <- spots[spots$barcode %in% sub("_.*$", "", (colnames(se_CID4465))), ]
identical(spots$barcode,sub("_.*$", "", (colnames(se_CID4465))))

# Ensure the barcodes align
rownames(spots) <- colnames(se_CID4465)

# Add the spatial coordinates to the Seurat object
se_CID4465@images$image@coordinates <- spots

# Confirm that coordinates have been added
head(se_CID4465@images$image@coordinates)

## ADD signatures
### ----  Antiox signature ---- 
oxstress <- read.csv("/Volumes/Lola/ONEDRIVE/Projects/AURIA/00.Data/gex/OxStressSig.csv")

# > OXSTRESS
sign_name <- "oxstress"
se_CID4465 <- AddModuleScore_UCell(se_CID4465, features = sign_list, name = '_ucell')

FeaturePlot(se_CID4465, features =paste0(sign_name, '_ucell'), label = TRUE, repel = TRUE)

FeatureOverlay(se_CID4465, features =paste0(sign_name, '_ucell'), cols = c("navyblue", "cyan", "yellow", "red", "darkred"))

a <- SpatialFeaturePlot(se_CID4465, features = "NQO1", slot = "counts")
b <- SpatialFeaturePlot(se_CID4465, features = "GCLC", slot = "counts")
c <- SpatialFeaturePlot(se_CID4465, features = "VIM", slot = "counts")

wrap_plots(a, b,c)

Idents(se_CID4465) <- "Classification"
VlnPlot(se_CID4465, features = paste0(sign_name, '_ucell'))

ggplot(as.data.frame(se_CD4465@meta.data), aes(x = get("Classification"), y =get(paste0(sign_name, '_ucell')))) +
  geom_violin() +
  xlab("Classification") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, colour = "black") +
  ylab(sign_name)

p1 <- SpatialFeaturePlot(se_CID4465, features =paste0(sign_name, '_ucell'), ncol=ceiling(length(sign_list)/3)) + theme(legend.position = "right")
p2 <- SpatialDimPlot(se_CID4465, group.by = "Classification")

# > GCLCVIM
GCLCVIM <- readRDS("/Volumes/Lola/ONEDRIVE/Projects/AURIA/02.Results/GCLCVIM_TUMOR.rds")
GCLCVIM_list <- list(GCLCVIM)  
names(GCLCVIM_list) <-"GCLCVIM"
sign_name <- "GCLCVIM"

se_CID4465 <- AddModuleScore_UCell(se_CID4465, features = GCLCVIM_list, name = '_ucell')
FeaturePlot(se_CID4465, features =paste0(sign_name, '_ucell'), label = TRUE, repel = TRUE)

FeatureOverlay(se_CID4465, features =paste0(sign_name, '_ucell'), cols = c("navyblue", "cyan", "yellow", "red", "darkred"))

Idents(se_CID4465) <- "Classification"
VlnPlot(se_CID4465, features = paste0(sign_name, '_ucell'))

ggplot(as.data.frame(se_CID4465@meta.data), aes(x = get("Classification"), y =get(paste0(sign_name, '_ucell')))) +
  geom_violin() +
  xlab("Classification") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, colour = "black") +
  ylab(sign_name)

p1 <- SpatialFeaturePlot(se_CID4465, features =paste0(sign_name, '_ucell'), ncol=ceiling(length(sign_list)/3)) + theme(legend.position = "right")
p2 <- SpatialDimPlot(se_CID4465, group.by = "Classification")

wrap_plots(p1, p2)


#### >> "CID44971F" ####

sample_CID44971 <- "CID44971"

se_CID44971 <- se.list[[sample_CID44971]]

se_CID44971 <- SCTransform(se_CID44971) %>%
  LoadImages(time.resolve = FALSE) 

se_CID44971 <- RunPCA(se_CID44971, assay = "SCT", verbose = FALSE)
se_CID44971 <- FindNeighbors(se_CID44971, reduction = "pca", dims = 1:30)
se_CID44971 <- FindClusters(se_CID44971, verbose = FALSE)
se_CID44971 <- RunUMAP(se_CID44971, reduction = "pca", dims = 1:30)

# Plot annotations
FeatureOverlay(se_CID44971,
               features = "Classification", 
               # pt.size = 0.75,
               # pt.alpha = 0.75,
               sampleids = 1,
               type = "raw")


DimPlot(se_CID44971, group.by = "Classification")


## Add coordinates
head(se_CID44971@images$image)

# The positions were not properly added, so I'm doing it again:
spotfile <- "/Volumes/Lola/ONEDRIVE/Projects/AURIA/00.Data/spatial/spatial/CID44971_spatial/tissue_positions_list.csv"
spots <- read.csv(spotfile, header = FALSE)
# Rename columns for clarity (adjust based on your file format)
colnames(spots) <- c("barcode", "tissue", "row", "col", "imagerow", "imagecol")
# Inspect the first few rows
head(spots)

spots <- spots[spots$barcode %in% sub("_.*$", "", (colnames(se_CID44971))), ]
identical(spots$barcode,sub("_.*$", "", (colnames(se_CID44971))))

# Ensure the barcodes align
rownames(spots) <- colnames(se_CID44971)

# Add the spatial coordinates to the Seurat object
se_CID44971@images$image@coordinates <- spots

# Confirm that coordinates have been added
head(se_CID44971@images$image@coordinates)

## ADD signatures

# > OXSTRESS
sign_name <- "oxstress"
se_CID44971 <- AddModuleScore_UCell(se_CID44971, features = sign_list, name = '_ucell')

FeaturePlot(se_CID44971, features =paste0(sign_name, '_ucell'), label = TRUE, repel = TRUE)

FeatureOverlay(se_CID44971, features =paste0(sign_name, '_ucell'), cols = c("navyblue", "cyan", "yellow", "red", "darkred"))

a<- SpatialFeaturePlot(se_CID44971, features = "NQO1", slot = "counts",cols = c("navyblue", "cyan", "yellow", "red", "darkred"))
b<- SpatialFeaturePlot(se_CID44971, features = "GCLC", slot = "counts")
c <- SpatialFeaturePlot(se_CID44971, features = "GCLC", slot = "counts")

wrap_plots(a, b,c)

Idents(se_CID44971) <- "Classification"
VlnPlot(se_CID44971, features = paste0(sign_name, '_ucell'))

ggplot(as.data.frame(se_CID44971@meta.data), aes(x = get("Classification"), y =get(paste0(sign_name, '_ucell')))) +
  geom_violin() +
  xlab("Classification") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, colour = "black") +
  ylab(sign_name)

p1 <- SpatialFeaturePlot(se_CID44971, features =paste0(sign_name, '_ucell'), ncol=ceiling(length(sign_list)/3)) + theme(legend.position = "right")
p2 <- SpatialDimPlot(se_CID44971, group.by = "Classification")

# > GCLCVIM
sign_name <- "GCLCVIM"
se_CID44971 <- AddModuleScore_UCell(se_CID44971, features = GCLCVIM_list, name = '_ucell')

FeatureOverlay(se_CID44971, features =paste0(sign_name, '_ucell'), cols = c("navyblue", "cyan", "yellow", "red", "darkred"))

Idents(se_CID44971) <- "Classification"
VlnPlot(se_CID44971, features = paste0(sign_name, '_ucell'))

ggplot(as.data.frame(se_CID44971@meta.data), aes(x = get("Classification"), y =get(paste0(sign_name, '_ucell')))) +
  geom_violin() +
  xlab("Classification") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, colour = "black") +
  ylab(sign_name)

p1 <- SpatialFeaturePlot(se_CID44971, features =paste0(sign_name, '_ucell'), ncol=ceiling(length(sign_list)/3)) + theme(legend.position = "right")
p2 <- SpatialDimPlot(se_CID44971, group.by = "Classification")

wrap_plots(p1, p2)


#####
# Add the signature to the deconvolution results of each sample 

deconvolution_results <- read.csv("/Users/spalominoe/Python_jobs/deconvolution_results.csv")

table(gsub(".*?-1_", "-1_", deconvolution_results$X))

#### se_CID44971 #### 
head(names(se_CID44971@active.ident))

head(deconvolution_results)
rownames(deconvolution_results) <- deconvolution_results$X

coord <- (names(se_CID44971@active.ident))

coord <- sub("_1$", "_CID44971-9-9", coord)

deconv_CID44971 <- deconvolution_results[ rownames(deconvolution_results) %in% coord,]

signatures <- se_CID44971$metadata


signatures <- data.frame(coord=  (names(se_CID44971@active.ident)) , oxstress=se_CID44971@meta.data$oxstress_ucell, GCLC_VIM=se_CID44971@meta.data$GCLCVIM_ucell)
signatures$coord <- sub("_1$", "_CID44971-9-9", signatures$coord)

signatures <- signatures[signatures$coord %in% rownames(deconv_CID44971),]

write.csv(signatures, "signatures_CID44971.csv")


#### se_CID4465 #### 
head(names(se_CID4465@active.ident))

head(deconvolution_results)
rownames(deconvolution_results) <- deconvolution_results$X

coord <- (names(se_CID4465@active.ident))

coord <- sub("_1$", "_CID4465-6-6", coord)

deconv_CID44971 <- deconvolution_results[ rownames(deconvolution_results) %in% coord,]

signatures <- se_CID4465$metadata


signatures <- data.frame(coord=  (names(se_CID4465@active.ident)) , oxstress=se_CID4465@meta.data$oxstress_ucell, GCLC_VIM=se_CID4465@meta.data$GCLCVIM_ucell, GCLCVIM_tumor=se_CID4465@meta.data$GCLCVIM_tumor_ucell)
signatures$coord <- sub("_1$", "_CID4465-6-6", signatures$coord)

signatures <- signatures[signatures$coord %in% rownames(deconv_CID44971),]

write.csv(signatures, "signatures_CID4465.csv")


signatures_TNBC<- rbind(signatures_CID4465,signatures_CID44971)

write.csv(signatures_TNBC "signatures_TNBC_spatial.csv")



