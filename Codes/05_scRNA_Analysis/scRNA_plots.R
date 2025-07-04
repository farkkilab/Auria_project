
# Created: 14/05/25
# Last modification: 06/06/25
# Author(s): Wenqing Chen

# ---- Load libraries ----

#remotes::install_github("REditorSupport/languageserver")
#remotes::install_github("nx10/httpgd")

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)

#install.packages("rstatix")
library(rstatix)


# ---- Figure 1 -----------------------------------------------------------
# ---- Load data ----
## ---- all data ----
# Load all scRNA data from the paper
scSeurat <- readRDS("/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/scSeurat.rds")

Idents(scSeurat) <- "celltype_major"
scSeurat <- subset(x=scSeurat, idents="Normal Epithelial", invert=TRUE )
# Remove HER2+ patients, only keep ER+ and TNBC
scSeurat_other <- subset(scSeurat, subset = subtype != "HER2+")

table(scSeurat@meta.data$celltype_major)
table(scSeurat_other@meta.data$Patient)
table(scSeurat@meta.data$Patient, scSeurat@meta.data$subtype)
table(scSeurat@meta.data$Patient, scSeurat@meta.data$celltype_major)
table(scSeurat@meta.data$Patient, scSeurat@meta.data$celltype_minor)
table(scSeurat@meta.data$celltype_major, scSeurat@meta.data$normal_cell_call)

## ---- TNBC data ----
scSeurat_TNBC <- subset(scSeurat_other, subset = subtype == "TNBC" )
table(scSeurat_TNBC@meta.data$Patient, scSeurat_TNBC@meta.data$celltype_major)
table(scSeurat_TNBC@meta.data$celltype_major, scSeurat_TNBC@meta.data$celltype_minor)


# ---- oxstress signatures scores ------------------------------------------
## ---- Data ----
# Select relevant metadata columns
#df_oxstress <- scSeurat@meta.data %>% select(Patient, celltype_major, subtype, oxstress)

df_oxstress <- scSeurat_other@meta.data %>% select(Patient, celltype_major, subtype, oxstress)

# Calculate the mean oxstress signatures scores per patient
df_oxstress <- df_oxstress %>%
  group_by(Patient) %>%
  mutate(mean_oxstress_patient = mean(oxstress, na.rm = TRUE)) %>%
  ungroup()

# Calculate the mean oxstress signatures scores per celltype per patient
df_oxstress <- df_oxstress %>%
  group_by(Patient, celltype_major) %>%
  mutate(mean_oxstress_patient_celltype = mean(oxstress, na.rm = TRUE)) %>%
  ungroup()

# All patients with ER+ and TNBC
summary_df <- df_oxstress %>%
  distinct(Patient, celltype_major, subtype, mean_oxstress_patient, mean_oxstress_patient_celltype)


## all patients (10 patients) [per patient plot]
summary_df_perpatient <- df_oxstress %>% distinct(Patient, subtype, mean_oxstress_patient)
## remove some patients (8 patients) [per patient plot]
summary_df_perpatient <- subset(summary_df_perpatient, !(Patient %in% c("CID4465", "CID3946", 'CID44041')))

## select cancer cells (10 patients) [per patient plot: only tumor]
summary_df_tumor <- subset(summary_df, (celltype_major == "Cancer Epithelial"))
summary_df_tumor <- summary_df_tumor %>% distinct(Patient, subtype, mean_oxstress_patient)
summary_df_tumor <- subset(summary_df_tumor, !(Patient %in% c("CID4465")))

# Only TNBC patients [per celltype per patient plot: only TNBC]
summary_df_TNBC <- subset(summary_df, subset = subtype == "TNBC")
summary_df_TNBC <- subset(summary_df_TNBC, !(Patient %in% c("CID4465", "CID3946", 'CID44041')))


## ---- Box Plots ----
# Color list
patient_colors <- c(
  # ER+
  "CID4461" = "#f7dfb8",
  "CID4463" = "#ffda5f",
  "CID4471" = "#ffba51",
  "CID4530N" = "#fd8302",
  "CID4535" = "#a26605",
  #"CID4040" = "#CD5118",
  "CID3941" = "#f2a4a4",
  "CID3948" = "#ff8989",
  "CID4067" = "#EF3B2C",
  "CID4290A" = "#A50F15",
  #"CID4398" = "#67000D",
  
  # HER2+
  #"CID3586" = "#DEEBF7",
  #"CID3921" = "#9ECAE1",
  #"CID45171" = "#6BAED6",
  #"CID3838" = "#3182BD",
  #"CID4066" = "#08519C",
  
  # TNBC
  #"CID44041" = "#E5F5E0",
  "CID4465" = "#C7E9C0",
  "CID4495" = "#d4ebd1",
  "CID44971" = "#A1D99B",
  "CID44991" = "#41AB5D",
  "CID4513" = "#9ECAE1",
  "CID4515" = "#6BAED6",
  "CID4523" = "#3182BD",
  #"CID3946" = "#238B45",
  "CID3963" = "#08519C"
)

subtype_color <- c("ER+" = "#BFE4EE", "HER2+" = "#FBD178", "TNBC" = "#F9C3BF")


celltype_major_color <- c(
  "Endothelial" = "#abdce7", "CAFs" = "#b6c3f1", "PVL" = "#a3d7bf",
  "B-cells" = "#ddd69d", "T-cells" = "#bcd59b", "Myeloid" = "#dbbce2",
  "Plasmablasts" = "#99b6d7", "Cancer Epithelial" = "#f1cdce" )

celltype_major_color_2 <- c(
  "Endothelial" = "#abdce7", "CAFs" = "#b6c3f1", "PVL" = "#a3d7bf",
  "B-cells" = "#ddd69d", "T-cells" = "#bcd59b", "Myeloid" = "#dbbce2",
  "Plasmablasts" = "#99b6d7", "Cancer Epithelial" = "#f1cdce" )


# Box plots

### ---- Per patient ----
# Order the patients label
subtype_ordered_patients <- summary_df_perpatient |>
  dplyr::distinct(Patient, subtype) |>
  dplyr::mutate(subtype = factor(subtype, levels = c("ER+", "HER2+", "TNBC"))) |>
  dplyr::arrange(subtype) |>
  dplyr::pull(Patient)
summary_df_perpatient$Patient <- factor(summary_df_perpatient$Patient, levels = subtype_ordered_patients)

ggplot(summary_df_perpatient, aes(x = subtype, y = mean_oxstress_patient)) +
  geom_boxplot(aes(fill = subtype), outlier.shape = NA, color = "black") +
  geom_jitter(aes(color = Patient), width = 0.2, size = 2, alpha = 0.8) +
  scale_color_manual(values = patient_colors) +
  scale_fill_manual(values = subtype_color) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("ER+", "HER2+"), c("ER+", "TNBC"), c("HER2+", "TNBC")),
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
ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/oxstress_per_patient.pdf", width = 6, height = 5, dpi = 300)
#ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/oxstress_per_patient_filter.pdf", width = 6, height = 5, dpi = 300)

# stat try
stat.test <- summary_df_perpatient %>%
  rstatix::kruskal_test(mean_oxstress_patient ~ subtype)

summary_df_perpatient %>%
  dunn_test(mean_oxstress_patient ~ subtype, p.adjust.method = "BH")

#### ---- only tumor in ER+ and TNBC ----
# Order the patients label
# !!!!! Use !!!!!
subtype_ordered_patients <- summary_df_tumor |>
  dplyr::distinct(Patient, subtype) |>
  dplyr::mutate(subtype = factor(subtype, levels = c("ER+", "TNBC"))) |>
  dplyr::arrange(subtype) |>
  dplyr::pull(Patient)
summary_df_tumor$Patient <- factor(summary_df_tumor$Patient, levels = subtype_ordered_patients)

ggplot(summary_df_tumor, aes(x = subtype, y = mean_oxstress_patient)) +
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
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )
ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/oxstress_per_patient_tumor.pdf", width = 4, height = 6, dpi = 300)


### ---- Per celltype per patient ----
#### ---- All patients ----
# Order the patients label
subtype_ordered_patients <- summary_df |>
  dplyr::distinct(Patient, subtype) |>
  dplyr::mutate(subtype = factor(subtype, levels = c("ER+", "HER2+", "TNBC"))) |>
  dplyr::arrange(subtype) |>
  dplyr::pull(Patient)
summary_df$Patient <- factor(summary_df$Patient, levels = subtype_ordered_patients)


ggplot(summary_df, aes(x = celltype_major, y = mean_oxstress_patient_celltype)) +
  geom_boxplot(aes(fill = celltype_major), outlier.shape = NA, color = "black") +
  geom_jitter(aes(color = Patient), width = 0.2, size = 2, alpha = 0.8) +
  scale_color_manual(values = patient_colors) +
  scale_fill_manual(values = celltype_major_color) +
  labs(
    title = "Oxidative Stress per Cell Type (Per Patient)",
    x = "Cell Type",
    y = "Mean Oxidative Stress Score"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )
ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/oxstress_per_patient_celltype.pdf", width = 7, height = 7, dpi = 300)

comparisons <- combn(unique(summary_df$celltype_major), 2, simplify = FALSE)
kruskal_res <- summary_df %>%
  kruskal_test(mean_oxstress_patient_celltype ~ celltype_major)
dunn_res <- summary_df %>%
  dunn_test(mean_oxstress_patient_celltype ~ celltype_major, p.adjust.method = "BH")


#### ---- TNBC patients ----
# Per celltype per patient
# Order the patients label
subtype_ordered_patients <- summary_df_TNBC |>
  dplyr::distinct(Patient, subtype) |>
  dplyr::mutate(subtype = factor(subtype, levels = c("TNBC"))) |>
  dplyr::arrange(subtype) |>
  dplyr::pull(Patient)
summary_df_TNBC$Patient <- factor(summary_df_TNBC$Patient, levels = subtype_ordered_patients)

#########################################
unique(summary_df_TNBC$celltype_major)

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

ggplot(summary_df_TNBC, aes(x = celltype_major, y = mean_oxstress_patient_celltype)) +
  geom_boxplot(aes(fill = celltype_major), outlier.shape = NA, color = "black") +
  geom_jitter(aes(color = Patient), width = 0.2, size = 2, alpha = 0.8) +
  scale_color_manual(values = patient_colors) +
  scale_fill_manual(values = celltype_major_color) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, # wilcoxon sum rank test (no paired samples)
                     comparisons = my_comparisons, p.adjust.method = "holm",
                     label = "p.signif", tip.length = 0.02, hide.ns = TRUE,
                     step.increase = 0.1) +
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

######################################
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

ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/oxstress_per_patient_celltype_TNBC.pdf", width = 7, height = 7, dpi = 300)
ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/oxstress_per_patient_celltype_TNBC_sig.pdf", width = 7, height = 7, dpi = 300)



##### ----Cell counts ----

#calculate how many cells in each celltype per patient
table(scSeurat_TNBC@meta.data$Patient, scSeurat_TNBC@meta.data$celltype_major)
table(scSeurat_TNBC@meta.data$Patient, scSeurat_TNBC@meta.data$celltype_minor)
table(scSeurat_TNBC@meta.data$celltype_major, scSeurat_TNBC@meta.data$celltype_minor)
unique(scSeurat_TNBC@meta.data$celltype_minor)

scSeurat_TNBC_filter <- subset(scSeurat_TNBC, !(Patient %in% c("CID4465", "CID3946", 'CID44041')))


df <- as.data.frame(table(scSeurat_TNBC_filter@meta.data$Patient, scSeurat_TNBC_filter@meta.data$celltype_major))
colnames(df) <- c("Patient", "CellType", "Count")
# Plot
ggplot(df, aes(x = Patient, y = Count, fill = CellType)) +
  geom_bar(stat = "identity", outlier.shape = NA, color = "#354E6B") +
  #geom_col(position = "dodge") +
  scale_fill_manual(values = celltype_major_color) +
  theme_classic() +
  labs(title = "Cell Counts per Cell Type per Patient",
       x = "Patient",
       y = "Cell Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/cellCounts_per_patient_celltype_TNBC.pdf", width = 5, height = 4, dpi = 300)

##### ----violin plot ----
scSeurat_TNBC_filter <- subset(scSeurat_TNBC, !(Patient %in% c("CID4465", "CID3946", 'CID44041')))
scSeurat_TNBC_filter_tumor <- subset(scSeurat_TNBC_filter, (celltype_major == "Cancer Epithelial"))
oxstress_df <- data.frame(
  Patient = scSeurat_TNBC_filter_tumor@meta.data$Patient,
  OxStress = scSeurat_TNBC_filter_tumor@meta.data$oxstress
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







# ---- Figure 2 -----------------------------------------------------------
# ---- load data ----
rds_data <- readRDS("/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/TNBC180225.rds")

# ---- GCLC VIM expression ----


# Parameters
markers <- c("GCLC", "VIM")
group_col <- "GCLCVIM"          

# Get normalized data
expr_data <- GetAssayData(rds_data, assay = "RNA", layer = "data")

# Subset expression for the markers
expr_subset <- expr_data[markers, ]
# Transpose and convert to data frame (cells as rows)
expr_df <- as.data.frame(t(as.matrix(expr_subset)))
# Add group info
expr_df[[group_col]] <- rds_data@meta.data[[group_col]]
expr_df <- expr_df %>%
  mutate(
    GCLC_Group = ifelse(str_sub(GCLCVIM, 1, 1) == "-", "GCLC_Neg", "GCLC_Pos"),
    VIM_Group = ifelse(str_sub(GCLCVIM, 2, 2) == "-", "VIM_Neg", "VIM_Pos")
  )



# GCLC expression in VIM group
GCLC_expr_df <- expr_df[c("GCLC", "VIM_Group")]

GCLC_expr_long <- GCLC_expr_df %>%
  pivot_longer(cols = "GCLC", names_to = "gene", values_to = "expression")

GCLC_summary  <- GCLC_expr_long %>%
  group_by(VIM_Group, gene) %>%
  summarise(
    avg_expression = mean(expression),
    pct_expressed = mean(expression > 0)  * 100
  ) %>%
  ungroup()

# stats test
#install.packages("coin")
# wilcox rank-sum test for expression
wilcox_test <- wilcox.test(GCLC ~ VIM_Group, data = GCLC_expr_df)
wilcox_p <- signif(wilcox_test$p.value, 3)
#wilcox_eff  <- wilcox_effsize(GCLC ~ VIM_Group, data = GCLC_expr_df)
# fisher test for percentage
GCLC_per_df <- GCLC_expr_df
GCLC_per_df <- GCLC_per_df %>%
  mutate(expressed = ifelse(GCLC > 0, "Yes", "No"))
table_expr <- table(GCLC_per_df$expressed, GCLC_per_df$VIM_Group)
fisher_test <- fisher.test(table_expr)
fisher_p <- signif(fisher_test$p.value, 3)

# Bubble plot
df_GCLC_summary <- data.frame(GCLC_summary)

p <- ggplot(df_GCLC_summary, aes(x = gene, y = VIM_Group)) +
  geom_point(aes(size = pct_expressed, fill = avg_expression),
             shape = 21, color = "black", stroke = 0.8) +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient(low = "white", high = "red") + 
  theme_minimal() +
  labs(
    title = "",
    x = "",
    y = "",
    size = "% Expressed",
    fill = "Avg Expression"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),                           
    axis.line = element_line(color = "black", size = 0.8)       
  )
combine_p <- p + annotate("text", 
             x = 0.5, y = 2.5, 
             label = paste0("Expression p = ", wilcox_p, 
                            "\nPercentage p = ", fisher_p),
             hjust = 0, size = 3)
ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/GCLC_in_VIMGroup.pdf",
       plot = combine_p,
       width = 4, height = 5, dpi = 300)


# VIM expression in GCLC group
VIM_expr_df <- expr_df[c("VIM", "GCLC_Group")]
VIM_expr_long <- VIM_expr_df %>%
  pivot_longer(cols = "VIM", names_to = "gene", values_to = "expression")
VIM_summary <- VIM_expr_long %>%
  group_by(GCLC_Group, gene) %>%
  summarise(
    avg_expression = mean(expression),
    pct_expressed = mean(expression > 0) * 100
  ) %>%
  ungroup()

# stats test
# wilcox rank-sum test for expression
wilcox_test <- wilcox.test(VIM ~ GCLC_Group, data = VIM_expr_df)
wilcox_p <- signif(wilcox_test$p.value, 3)
# fisher test for percentage
VIM_per_df <- VIM_expr_df
VIM_per_df <- VIM_per_df %>%
  mutate(expressed = ifelse(VIM > 0, "Yes", "No"))
table_expr <- table(VIM_per_df$expressed, VIM_per_df$GCLC_Group)
fisher_test <- fisher.test(table_expr)
fisher_p <- signif(fisher_test$p.value, 3)

# Bubble plot
df_VIM_summary <- data.frame(VIM_summary)

p <- ggplot(df_VIM_summary, aes(x = gene, y = GCLC_Group)) +
  geom_point(aes(size = pct_expressed, fill = avg_expression),
             shape = 21, color = "black", stroke = 0.8) +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    x = "",
    y = "",
    size = "% Expressed",
    fill = "Avg Expression"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),                           
    axis.line = element_line(color = "black", size = 0.8)       
  )
combine_p <- p + annotate("text", 
                          x = 0.5, y = 2.5, 
                          label = paste0("Expression p = ", wilcox_p, 
                                         "\nPercentage p = ", fisher_p),
                          hjust = 0, size = 3)
ggsave("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/Oxstress/VIM_in_GCLCGroup.pdf",
       plot = combine_p,
       width = 5, height = 5, dpi = 300)



# ---- Volcano plot (GSEA) ----
library(ggplot2)
library(ggrepel)
library(dplyr)

## ---- GCLC VIM ----
plot_df <- read.csv('/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/GSEA/Oxstress_group_GSEA.csv')

plot_df <- plot_df %>%
  mutate(`-log10(p.adj)` = -log10(ifelse(FDR.q.val == 0, 1e-4, FDR.q.val)))

plot_df <- subset(plot_df, subset = Group != "Oxstress_High_Oxstress_Low" )
unique_groups <- unique(plot_df$Group)

for (grp in unique_groups) {
  group_data <- plot_df %>% filter(Group == grp)
  
  top_pathways <- group_data %>%
    filter(significance != "not_sig") %>%
    arrange(desc(`-log10(p.adj)`)) %>%
    slice_head(n = 16)
  
  p <- ggplot(group_data, aes(x = NES, y = `-log10(p.adj)`)) +
    geom_point(aes(color = significance), size = 3, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    geom_label_repel(
      data = top_pathways,
      aes(label = gsub("HALLMARK_", "", Term), fill = significance),
      size = 2.5,
      box.padding = 0.3,
      max.overlaps = 100,
      segment.size = 0.2,
      segment.color = "grey50",
      fontface = "bold"
    ) +
    scale_fill_manual(values = c("up" = "#FBEAEA", "down" = "#EAF3FB")) +
    scale_color_manual(values = c("up" = "#D62728", "down" = "#1F77B4", "not_sig" = "grey")) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste("Group:", grp),
      x = "Normalized Enrichment Score (NES)",
      y = "-log10(FDR q-value)",
      color = "Significance"
    )
  
  print(p)
  ggsave(
    filename = paste0("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/GSEA_volcano_plot_",grp,".pdf"),
    plot = p,
    width = 5.5, height = 8.5, dpi = 300
  )
}

## ---- Oxstress ----
plot_df <- read.csv('/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/GSEA/Oxstress_group_GSEA.csv')
plot_df <- plot_df %>%
  mutate(`-log10(p.adj)` = -log10(ifelse(FDR.q.val == 0, 1e-4, FDR.q.val)))

plot_df <- subset(plot_df, subset = Group == "Oxstress_High_Oxstress_Low" )
unique_groups <- unique(plot_df$Group)

for (grp in unique_groups) {
  group_data <- plot_df %>% filter(Group == grp)
  
  top_pathways <- group_data %>%
    filter(significance != "not_sig") %>%
    arrange(desc(`-log10(p.adj)`)) %>%
    slice_head(n = 25)
  
  p <- ggplot(group_data, aes(x = NES, y = `-log10(p.adj)`)) +
    geom_point(aes(color = significance), size = 3, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    geom_label_repel(
      data = top_pathways,
      aes(label = gsub("HALLMARK_", "", Term), fill = significance),
      size = 2.5,
      box.padding = 0.3,
      max.overlaps = 100,
      segment.size = 0.2,
      segment.color = "grey50",
      fontface = "bold"
    ) +
    scale_fill_manual(values = c("up" = "#FBEAEA", "down" = "#EAF3FB")) +
    scale_color_manual(values = c("up" = "#D62728", "down" = "#1F77B4", "not_sig" = "grey")) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste("Group:", grp),
      x = "Normalized Enrichment Score (NES)",
      y = "-log10(FDR q-value)",
      color = "Significance"
    )
  
  print(p)
  ggsave(
    filename = paste0("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/GSEA_volcano_plot_",grp,".pdf"),
    plot = p,
    width = 7.5, height = 7.5, dpi = 300
  )
}












