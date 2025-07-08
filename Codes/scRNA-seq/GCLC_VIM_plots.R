
# Created: 14/05/25
# Last modification: 06/06/25
# Author(s): Wenqing Chen

# ---- Load libraries ----

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(rstatix)


# ---- load data ----
rds_data <- readRDS("/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/TNBC180225.rds")

# GCLC VIM expression
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


# ---- Figure 2I ----
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

# ---- Figure 2J ----
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

