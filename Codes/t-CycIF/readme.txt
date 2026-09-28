clinical_data_per_core.csv
Column names in this file are:

patient_id-AB19_1654: patient id that consists of TMA block name and patient id, for example AB19-1654_P00001.  If there is value NA, the 
data of the patient is not anymore available due to the 
biobank refusal.  
imageid: TMA slide ID
Annotation: Histological region sampled, e.g., Tumor center, lymph node metastasis, etc.
Stage: disease stage
core_imageid: core number together with image in which core can be found; e.g., "core1_19_02_01A" indicates the first core in image 19_02_01A
status: status at end of follow up; 0=alive , 1=deceased
time: time from sample to the end of follow-up, in years
rounded_age: patient's age at time of sampling, rounded in years

Data_AllMarkers_filtered_RCNs.csv
The columns "imageid", "core_imageid", "Annotation", "Stage", "patient_id_AB19_1654" correspond to the columns with the same names in the clinical_data_per_core.csv and can be used to map the clinical data to the single cell data.

The following columns provide the marker intensities for all segmented cells:
"S9.6", "HER3", "panCK", "SMA", "CD45", "Ki67_2", "CD4", "FOXP3", "CD8a", "CD11b", "CD68", "CD206", "H2ax", "Vimentin", "H3K27me3", "E.Cadherin", "GCLC", "TXNRD1", "NQO1", "pERK", "TXN", "pSTAT3", "pS6_235", "pRB", "GLUT1"

The location of each segmented cell is given by "X_centroid", "Y_centroid".

Additional column names in this file are:
CellID: resulting segmented cell identity
cores.x: identifier for the TMA cores (e.g., core77, core 42)
celltype: the designated cell type for each segmented cells (tumor, CD8+ T cells, CD4+ T cells, regulatory T cells, CD68+ macrophages, CD206+ macrophages, and other CD45+ immune cells, SMA+ stromal cells, VIM+ stromal cells, and other panCK-CD45- cells) 
Area: the area of the segmented cell in pixels
Eccentricity: the segmented cell's eccentricity
antioxidant_NQO1: the expression classification status of the NQO1 antioxidant marker for the cell
antioxidant_GCLC: the expression classification status of the GCLC antioxidant marker for the cell
antioxidant_TXNRD1: the expression classification status of the TXNRD1 antioxidant marker for the cell
antioxidant_GLUT1: the expression classification status of the GLUT1 marker for the cell
NQO1_GCLC_VIM: a combined phenotype string indicating the expression status (positive + or negative -) of NQO1, GCLC, and Vimentin (VIM)
celltype_2: a sub-clustered cell type annotation, often combining the base cell type with specific functional marker phenotypes
celltype_1: grouped cell type classification
Tumor_state: copy column of celltype_1
GCLC_VIM: a combined annotation merging the cell type with its corresponding GCLC and Vimentin expression profiles
CellCategory: major cell types including tumor, storm and immune
SC45px_neigh_kmeans_10: the spatial cluster ID assigned to the cell using K-means clustering (where k=10) based on the cellular composition within its 45-pixel radius spatial neighborhood
RCNs: Recurrent Cellular Neighborhood assignment of cell
