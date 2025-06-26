
# install.packages("ggpubr")

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(readr)


plot_signature_violin <- function(
    seurat_obj,
    signature_genes,
    signature_name = "signature",
    score_colname = "signature_score",
    group_var,
    subset_celltype = NULL,
    group_subset = NULL,
    group_colors = NULL,
    p.adjust.method = "BH",
    save_path = NULL,
    violin_save_path = NULL,
    width, height
) {
  # Step 1: Add module score
  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = signature_genes,
    name = signature_name
  )
  
  # Automatically get the score column (e.g., "signature1")
  score_col_raw <- paste0(signature_name, "1")
  
  # Rename to user-defined score column name
  seurat_obj@meta.data[[score_colname]] <- seurat_obj@meta.data[[score_col_raw]]
  seurat_obj@meta.data[[score_col_raw]] <- NULL
  
  # Save the result if a path is provided
  if (!is.null(save_path)) {
    write_rds(seurat_obj, file = save_path)
  }
  
  # Step 2: Optionally subset celltype (e.g., "Cancer Epithelial")
  if (!is.null(subset_celltype)) {
    seurat_obj <- subset(seurat_obj, subset = celltype_major == subset_celltype)
  }
  
  # Step 3: Fetch data
  df <- FetchData(seurat_obj, vars = c(score_colname, group_var))
  df[[group_var]] <- as.character(df[[group_var]]) 
  
  # Step 4: Optionally subset groups
  if (!is.null(group_subset)) {
    df <- df[df[[group_var]] %in% group_subset, ]
    df[[group_var]] <- droplevels(as.factor(df[[group_var]]))
  }
  
  # Step 5: Pairwise Wilcoxon test
  pairwise_res <- pairwise.wilcox.test(
    x = df[[score_colname]],
    g = df[[group_var]],
    p.adjust.method = p.adjust.method
  )
  # Print and save the test result
  print(pairwise_res)

  
  # Step 6: Convert p-values to annotation format
  # Convert p-values to annotation format
  pval_df <- as.data.frame(as.table(pairwise_res$p.value)) %>%
    na.omit() %>%
    mutate(
      group1 = as.character(Var1),
      group2 = as.character(Var2),
      label = cut(Freq,
                       breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("***", "**", "*", "ns")),
      y.position = seq(
        from = max(df[[score_colname]], na.rm = TRUE) + 0.2,
        by = 0.1,
        length.out = dplyr::n()
           ))
  
  # Generate color palette based on number of groups
  group_levels <- unique(df[[group_var]])
  if (!is.null(group_colors)) {
    # Use user defined colors
    fill_colors <- group_colors[group_levels]
  } else {
    # By default, the hue palette is used. 
    fill_colors <- scales::hue_pal()(length(group_levels))
    names(fill_colors) <- group_levels
  }

  # Prepare comparisons list
  comparisons_list <- split(pval_df[, c("group1", "group2")], seq(nrow(pval_df)))
  comparisons_list <- lapply(comparisons_list, unlist) 
  
  # Step 7: Violin plot with stars
  plot <- ggplot(df, aes(x = .data[[group_var]], y = .data[[score_colname]], fill = .data[[group_var]])) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    scale_fill_manual(values = fill_colors) +
    ggsignif::geom_signif(
      comparisons = comparisons_list,
      annotations = pval_df$label,
      y_position = pval_df$y.position,
      tip_length = 0.01,
      textsize = 5
    ) +
    labs(x = group_var, y = score_colname) +
    theme_minimal()
  if (!is.null(violin_save_path)) {
    ggsave(violin_save_path, plot = plot, width = width, height = height)
  }
  
  return(plot)
}
