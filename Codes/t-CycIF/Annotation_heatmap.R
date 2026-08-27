

# Created: 19/08/25
# Last modification: 23/06/26
# Author(s): Wenqing Chen

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

clinical_data <- read.csv('/Users/wenqchen/Desktop/Projects/Auria/Data/t-CycIF/clinical_data_per_patient.csv')

clinical_data$Stage <- factor(
  clinical_data$Stage, 
  levels = c("1", "2", "3", "4", "unknown")
)

regions <- c("Tumor_center", "Invasive_border", "Inflammation", "Lymph_node_met")

clinical_data <- clinical_data %>%
  mutate(across(all_of(regions), ~ factor(.x, levels = c("Yes", "No")))) %>%
  arrange(Tumor_center, Invasive_border, Inflammation, Lymph_node_met)


col_list = list(
  Stage = c("1"="#58B668","2"="#046B38","3"="red","4"="darkred","unknown"="#808080"),
  Rounded_age = colorRamp2(c(20,50,100), c("#073162","#89BED9","#F5F8FA")),
  Survival_time = colorRamp2(c(0,10,20), c("#C41B7E","white","#FCE1EF")),
  Survival_status = c("0"="#cec3e7","1"="#9473cc"),
  Tumor_center = c("Yes"="red","No"="pink"),
  Invasive_border = c("Yes"="red","No"="pink"),
  Inflammation = c("Yes"="red","No"="pink"),
  Lymph_node_met = c("Yes"="red","No"="pink")
)

ha_sorted <- HeatmapAnnotation(
  df = clinical_data %>% select(
    Stage, Rounded_age, Survival_time, Survival_status,
    Tumor_center, Invasive_border, Inflammation, Lymph_node_met
  ),
  col = col_list,
  annotation_name_side = "left",
  gap = unit(1, "mm"),
  gp = gpar(col = "black")
)

n <- nrow(clinical_data)

pdf("/Users/wenqchen/Desktop/Projects/Auria/Clinical/final/annotation_heatmap_patient_sorted_by_regions.pdf", width = 30, height = 3)

Heatmap(
  matrix(0, nrow = 1, ncol = n), 
  top_annotation = ha_sorted,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  name = "dummy",
  height = unit(0, "mm"),
  col = c("0" = "white"),
  rect_gp = gpar(col = "white", lwd = 0)
)

dev.off()



