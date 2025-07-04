
# Created: 03/07/25
# Last modification: 03/07/25
# Author(s): Wenqing Chen

# ---- Volcano plot ----
library(ggplot2)
library(ggrepel)
library(dplyr)


## ---- Oxstress ----

# Define the pathways to label
label_pathways <- c(
  'INTERFERON_GAMMA_RESPONSE', 'INTERFERON_ALPHA_RESPONSE', 'IL2_STAT5_SIGNALING',
  'MYC_TARGETS_V1', 'GLYCOLYSIS', 'OXIDATIVE_PHOSPHORYLATION', 'FATTY_ACID_METABOLISM'
)

plot_df <- read.csv('/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/GSEA/Oxstress_group_GSEA.csv')
plot_df <- plot_df %>%
  mutate(`-log10(p.adj)` = -log10(ifelse(FDR.q.val == 0, 1e-4, FDR.q.val)))

plot_df <- subset(plot_df, subset = Group == "Oxstress_High_Oxstress_Low" )
unique_groups <- unique(plot_df$Group)

for (grp in unique_groups) {
  group_data <- plot_df %>% filter(Group == grp)
  
  top_pathways <- group_data %>%
    filter(gsub("HALLMARK_", "", Term) %in% label_pathways)
  
  p <- ggplot(group_data, aes(x = NES, y = `-log10(p.adj)`)) +
    geom_point(aes(color = significance), size = 3, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    geom_label_repel(
      data = top_pathways,
      aes(label = gsub("HALLMARK_", "", Term), fill = significance),
      size = 1.6,
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
      #title = paste("Group:", grp),
      x = "Normalized Enrichment Score (NES)",
      y = "-log10(FDR q-value)",
      color = "Significance"
    ) +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    )
  
  print(p)
  ggsave(
    filename = paste0("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/volcano/new_2_GSEA_volcano_plot_",grp,".pdf"),
    plot = p,
    width = 5, height = 4, dpi = 300
  )
}



## ---- GCLC VIM ----

### ---- GCLC+VIM+ VS GCLC+VIM- ----
label_pathways <- c(
  'EPITHELIAL_MESENCHYMAL_TRANSITION', 'HEDGEHOG_SIGNALING', 'INFLAMMATORY_RESPONSE',
  'TGF_BETA_SIGNALING', 'TNFA_SIGNALING_VIA_NFKB')

plot_df <- read.csv('/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/GSEA/Oxstress_group_GSEA.csv')

plot_df <- plot_df %>%
  mutate(`-log10(p.adj)` = -log10(ifelse(FDR.q.val == 0, 1e-4, FDR.q.val)))

plot_df <- subset(plot_df, subset = Group == "GCLC+VIM+_GCLC+VIM-" )
unique_groups <- unique(plot_df$Group)

for (grp in unique_groups) {
  group_data <- plot_df %>% filter(Group == grp)
  
  top_pathways <- group_data %>%
    filter(gsub("HALLMARK_", "", Term) %in% label_pathways)
  
  p <- ggplot(group_data, aes(x = NES, y = `-log10(p.adj)`)) +
    geom_point(aes(color = significance), size = 3, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    geom_label_repel(
      data = top_pathways,
      aes(label = gsub("HALLMARK_", "", Term), fill = significance),
      size = 1.6,
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
      #title = paste("Group:", grp),
      x = "Normalized Enrichment Score (NES)",
      y = "-log10(FDR q-value)",
      color = "Significance"
    ) +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    )
  
  print(p)
  ggsave(
    filename = paste0("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/volcano/new_GSEA_volcano_plot_",grp,".pdf"),
    plot = p,
    width = 3.5, height = 3.5, dpi = 300
  )
}

### ---- GCLC+VIM+ VS GCLC-VIM+ ----
label_pathways <- c(
  'ALLOGRAFT_REJECTION', 'INTERFERON_GAMMA_RESPONSE', 'IL6_JAK_STAT3_SIGNALING')

plot_df <- read.csv('/Users/wenqchen/Desktop/Projects/Auria/Data/scRNA/GSEA/Oxstress_group_GSEA.csv')

plot_df <- plot_df %>%
  mutate(`-log10(p.adj)` = -log10(ifelse(FDR.q.val == 0, 1e-4, FDR.q.val)))

plot_df <- subset(plot_df, subset = Group == "GCLC+VIM+_GCLC-VIM+" )
unique_groups <- unique(plot_df$Group)

for (grp in unique_groups) {
  group_data <- plot_df %>% filter(Group == grp)
  
  top_pathways <- group_data %>%
    filter(gsub("HALLMARK_", "", Term) %in% label_pathways)
  
  p <- ggplot(group_data, aes(x = NES, y = `-log10(p.adj)`)) +
    geom_point(aes(color = significance), size = 3, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    geom_label_repel(
      data = top_pathways,
      aes(label = gsub("HALLMARK_", "", Term), fill = significance),
      size = 1.6,
      box.padding = 0.3,
      max.overlaps = 100,
      segment.size = 0.2,
      segment.color = "grey50",
      fontface = "bold"
    ) +
    scale_fill_manual(values = c("up" = "#FBEAEA", "down" = "#EAF3FB")) +
    scale_color_manual(values = c("up" = "#D62728", "down" = "#1F77B4", "not_sig" = "grey")) +
    theme_minimal(base_size = 12) +
      labs(
        #title = paste("Group:", grp),
        x = "Normalized Enrichment Score (NES)",
        y = "-log10(FDR q-value)",
        color = "Significance"
      ) +
      theme(
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8)
      )
  
  print(p)
  ggsave(
    filename = paste0("/Users/wenqchen/Desktop/Projects/Auria/Plots/scRNA/volcano/new_GSEA_volcano_plot_",grp,".pdf"),
    plot = p,
    width = 3.5, height = 3.5, dpi = 300
  )
}


