
# Created: 14/05/25
# Last modification: 06/06/25
# Author(s): Wenqing Chen

library(Seurat)
library(dplyr)
library(ggplot2)

custom_colors <- c("++" = '#F0776D', '-+' = '#a9d9bb')
#'+-' = '#FBD2CB', '--' = '#337bac'

rds_data <- readRDS("/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/TNBC180225.rds")

anastasis_signature <- list(c('EGR1', 'FOSB', 'SNAI1', 'NR4A2', 'ATF3', 'DUSP1', 'JUN', 'EDN1', 'PPP1R15A', 'ADM',
                              'ZFP36', 'NR4A1', 'GADD45B', 'SOCS3', 'OTUD1', 'DUSP10', 'NR1D1', 'FOS', 'IER2', 'LOC284454',
                              'IER5', 'LOC644656', 'ID2', 'RSRP1', 'KLF11', 'NRARP', 'PIM1', 'CITED2', 'KLF10', 'KMT2E-AS1',
                              'CSRNP1', 'TRIM52', 'ING1', 'RHOB', 'CREBRF', 'ADRB2', 'JUNB', 'LINC00673', 'PHF13', 'ZBTB10', 'CEBPB'))


source("~/Documents/GitHub/Auria/Codes/05_scRNA_Analysis/signature_comparison.R")
plot_signature_violin(
  seurat_obj = rds_data,
  signature_genes = anastasis_signature,
  signature_name = "anastasis",
  score_colname = "signature_anastasis",
  group_var = "GCLCVIM",
  group_colors = custom_colors,
  subset_celltype = "Cancer Epithelial",
  group_subset = c("++", "-+"),
  p.adjust.method = "BH",
  save_path = "/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/TNBC_signatures.rds",
  violin_save_path = "/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Signature_score/Anastasis_score.pdf",
  width = 4, height = 6
)

PCNA_signature <- list(c('LMNB1','NFE2','KIF20A','GTSE1','FECH','RRM2','MCM4','SMC4',
                              'ALAS2','CCNA2','CCNB1','GYPB','LSM6','RPIA','CHAF1A','POLE2',
                              'FBXO7','RPP30','NCAPG2','KEL','MCM2','GYPA','ASF1B','TOP2A',
                              'TRMT5','NUP37','NCAPD2','TACC3','RAD51AP1','FEN1','LYL1',
                              'KLF1','VRK1','RACGAP1','DNAJC9','CDCA3','PGD','TROAP','TRIM58',
                              'ADAMTS13','KLF15','CKS2','RFC4','TFDP1','NUP210','GATA1','RFWD3',
                              'HMGB2','ESPL1','CENPA','PRC1','SNRPD1','RFC3','RHCE','ZWINT',
                              'LBR','CDKN3','CDCA4','FOXM1','SHCBP1','MICB','GINS2','FBXO5',
                              'CKLF','NCAPD3','BPGM','RHD','TPX2','CKS1B','BIRC5','MCM5','MAD2L1',
                              'TRIM10','MCM6','PSMD9','GTPBP2','CDC20','KIF4A','CCNB2','RPA3',
                              'MCM3','TAL1','PPIH','KIF2C','PTTG1','KIF22','EPB42','CDCA8','HMBS',
                              'TIMELESS','AURKB','DTL','RHAG','PPBP','TYMS','CDT1','SNRPB','GINS1',
                              'UBE2C','SPTA1','AURKA','PLEK','LIG1','SNF8','ARID3A','PCNA','MKI67',
                              'MELK','NUSAP1','MCM7','NUDT1','TCF3','OIP5','BUB1B','H3F3A','PF4',
                              'WHSC1','APOBEC3B','KIAA0101','CDC2','HMGN2','SFRS2','CDC45L','DDX39',
                              'ORC6L','ERAF','C21orf45','MLF1IP','BZRPL1'))

rds_data_anastasis <- readRDS("/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/TNBC_signatures.rds")

plot_signature_violin(
  seurat_obj = rds_data_anastasis,
  signature_genes = PCNA_signature,
  signature_name = "proliferation",
  score_colname = "signature_proliferation",
  group_var = "GCLCVIM",
  group_colors = custom_colors,
  subset_celltype = "Cancer Epithelial",
  group_subset = c("++", "-+"),
  p.adjust.method = "BH",
  save_path = "/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/TNBC_signatures.rds",
  violin_save_path = "/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Signature_score/proliferation_score.pdf",
  width = 4, height = 6
)
