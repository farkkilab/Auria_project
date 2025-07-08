
# Created: 02/07/25
# Last modification: 02/07/25
# Author(s): Wenqing Chen

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
  "CID3948" = "#ff8989", "CID4067" = "#EF3B2C", "CID4290A" = "#A50F15"
  
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
# Load all scRNA data from the paper
### ---- ER+ & TNBC ----
scSeurat <- readRDS("/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/scSeurat_def.rds")

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
scSeurat_TNBC <- subset(scSeurat, subset = subtype == "TNBC" )
table(scSeurat_TNBC@meta.data$Patient, scSeurat_TNBC@meta.data$subtype)

## ---- Plots ----
### ---- Figure 1B ----
# Box plot per patient
# ER+ and TNBC
subtype_ordered_patients <- summary_df_tumor |>
  dplyr::distinct(Patient, subtype) |>
  dplyr::mutate(subtype = factor(subtype, levels = c("ER+", "TNBC"))) |>
  dplyr::arrange(subtype) |>
  dplyr::pull(Patient)
summary_df_tumor$Patient <- factor(summary_df_tumor$Patient, levels = subtype_ordered_patients)

ggplot(summary_df_tumor, aes(x = subtype, y = mean_oxstress_patient_celltype)) +
  geom_boxplot(aes(fill = subtype), outlier.shape = NA, color = "black") +
  geom_jitter(aes(color = Patient), width = 0.2, size = 2, alpha = 0.8) +
  scale_color_manual(values = patient_colors) +
  scale_fill_manual(values = subtype_color) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, # wilcoxon sum rank test (no paired samples)
                     comparisons = list(c("ER+", "TNBC")),
                     label = "p.signif") +
  labs(
    title = "Oxidative Stress per Patient Subtype",
    x = "Patient Subtype",
    y = "Mean Oxidative Stress Score"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )
ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/oxstress_per_patient_tumor.pdf", width = 4, height = 6, dpi = 300)

### ---- SFigure 1B ----
# per celltype in TNBC
subtype_ordered_patients <- summary_df_TNBC |>
  dplyr::distinct(Patient, subtype) |>
  dplyr::mutate(subtype = factor(subtype, levels = c("TNBC"))) |>
  dplyr::arrange(subtype) |>
  dplyr::pull(Patient)
summary_df_TNBC$Patient <- factor(summary_df_TNBC$Patient, levels = subtype_ordered_patients)

my_comparisons <- list(c("B-cells", "CAFs"), c("B-cells", "Cancer Epithelial"), c("B-cells", "Endothelial"),
                       c("B-cells", "Myeloid"), c("B-cells", "Plasmablasts"), c("B-cells", "PVL"), c("B-cells", "T-cells"),
                       c("CAFs", "Cancer Epithelial"), c("CAFs", "Endothelial"), c("CAFs", "Myeloid"),
                       c("CAFs", "Plasmablasts"), c("CAFs", "PVL"), c("CAFs", "T-cells"),
                       c("Cancer Epithelial", "Endothelial"), c("Cancer Epithelial", "Myeloid"),
                       c("Cancer Epithelial", "Plasmablasts"), c("Cancer Epithelial", "PVL"), c("Cancer Epithelial", "T-cells"),
                       c("Endothelial", "Myeloid"), c("Endothelial", "Plasmablasts"), c("Endothelial", "PVL"), c("Endothelial", "T-cells"),
                       c("Myeloid", "Plasmablasts"), c("Myeloid", "PVL"), c("Myeloid", "T-cells"),
                       c("Plasmablasts", "PVL"), c("Plasmablasts", "T-cells"),
                       c("PVL", "T-cells")
)

comparisons <- combn(unique(summary_df_TNBC$celltype_major), 2, simplify = FALSE)
kruskal_res <- summary_df_TNBC %>%
  kruskal_test(mean_oxstress_patient_celltype ~ celltype_major)

dunn_res <- summary_df_TNBC %>%
  dunn_test(mean_oxstress_patient_celltype ~ celltype_major, p.adjust.method = "BH")
dunn_sig <- dunn_res %>% 
  filter(p.adj <= 0.05) %>%
  mutate(y.position = seq(0.22, 0.27, length.out = n()))

ggplot(summary_df_TNBC, aes(x = celltype_major, y = mean_oxstress_patient_celltype)) +
  geom_boxplot(aes(fill = celltype_major), outlier.shape = NA, color = "black") +
  geom_jitter(aes(color = Patient), width = 0.2, size = 2, alpha = 0.8) +
  scale_color_manual(values = patient_colors) +
  scale_fill_manual(values = celltype_major_color) +
  geom_signif(
    data = dunn_sig,
    aes(xmin = group1, xmax = group2, annotations = p.adj.signif, y_position = y.position),
    manual = TRUE
  ) +
  labs(
    title = "Oxidative Stress per Cell Type (Per Patient)",
    x = "Cell Type",
    y = "Mean Oxidative Stress Score"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )
ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/oxstress_per_patient_celltype_TNBC_sig.pdf", width = 7, height = 7, dpi = 300)

### ---- Figure 1C ----
# violin plot of oxstress
# per cell in TNBC
scSeurat_TNBC_tumor <- subset(scSeurat_TNBC, (celltype_major == "Cancer Epithelial"))
oxstress_df <- data.frame(
  Patient = scSeurat_TNBC_tumor@meta.data$Patient,
  OxStress = scSeurat_TNBC_tumor@meta.data$oxstress
)

ggplot(oxstress_df, aes(x = Patient, y = OxStress, fill = Patient)) +
  geom_violin(trim = FALSE, alpha = 0.5, width = 0.9) +
  #geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.3, color = "black") +
  geom_jitter(aes(color = Patient),
              width = 0.35, size = 0.05, alpha = 0.3, show.legend = FALSE) + 
  theme_classic() +
  labs(title = "Oxidative Stress (Tumor)", x = "Patients", y = "Expression") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/Oxstress_expression_Tumor_TNBC.pdf",
       width = 6, height = 4, dpi = 300)


# ---- TNBC file ----
## ---- Load data ----
TNBC <- readRDS("/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/TNBC180225.rds")
TNBC <- subset(TNBC, !(Patient %in% c("CID4465")))
table(TNBC@meta.data$Patient, TNBC@meta.data$celltype_major)

## ---- Plots ----
### ---- Figure 3H ----
# violin plot of GCLC_VIM 
# per cell in TNBC
TNBC_tumor <- subset(TNBC, (celltype_major == "Cancer Epithelial"))
GCLCVIM_df <- data.frame(
  Patient = TNBC_tumor@meta.data$Patient,
  GCLCVIM_tumor = TNBC_tumor@meta.data$signature_GCLCVIMTUMOR
)

ggplot(GCLCVIM_df, aes(x = Patient, y = GCLCVIM_tumor, fill = Patient)) +
  geom_violin(trim = FALSE, alpha = 0.5, width = 1.0) +
  #geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.3, color = "black") +
  geom_jitter(aes(color = Patient),
              width = 0.25, size = 0.05, alpha = 0.3, show.legend = FALSE) + 
  theme_classic() +
  labs(title = "GCLC_VIM (Tumor)", x = "Patients", y = "Expression") + 
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/GCLC_VIM_Tumor_TNBC.pdf",
       width = 4, height = 4, dpi = 300)


## ---- correlation ----

### ---- SFigure 3E ----
# only tumor cells

oxstress_GV_tumor_df <- subset(oxstress_GV_df, subset = celltype_major == "Cancer Epithelial")

# hexagonal binning plot
ggplot(oxstress_GV_tumor_df, aes(x = oxstress, y = signature_GCLCVIMTUMOR)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c() +
  labs(
    title = "All tumor cells of TNBC: OxStress vs GCLCVIM",
    x = "OxStress score", y = "GCLCVIM score"
  ) +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = "p-value < 2.2e-16", hjust = 1.1, vjust = 1.2, size = 5)
ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/GCLCVIM_Oxstress_TNBC_tumor.pdf",
       width =10, height = 7, dpi = 300)

