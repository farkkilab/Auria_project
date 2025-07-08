
# Created: 14/05/25
# Last modification: 06/06/25
# Author(s): Wenqing Chen

# ---- Load libraries ----
library(ggplot2)
library(ggrepel)
library(dplyr)

# ---- Volcano plot (GSEA) ----
## ---- Figure 3A & 3B ----
# GCLC VIM
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

## ---- Figure 1E ----
# Oxstress groups
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

