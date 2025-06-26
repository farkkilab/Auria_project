

# Created: 09/06/25
# Last modification: 09/06/25
# Author(s): Wenqing Chen



library(dplyr)
library(tidyr)
#install.packages(c("survival", "survminer"))
library(survival)
library(survminer)
#install.packages("forestplot")
library(forestplot)
#install.packages("meta")
library(meta)

data <- read.csv('/Users/wenqchen/Desktop/Projects/Auria/Data/t-CycIF/Data_AllMarkers_filtered_RCNs.csv')

df_clinical_GV <- read.csv('/Users/wenqchen/Desktop/Projects/Auria/Data/t-CycIF/clinical_data_GCLC_VIM_quarter.csv')

clinical_data <- read.csv('/Users/wenqchen/Desktop/Projects/Auria/Data/t-CycIF/clinical_data_per_core.csv')


# ---- Cox Harzard Model -----------------------------------------------------------
## ---- Multiple cox ----
### ---- Functional marker ----

data$celltype_EMT <- data$celltype
data[which(data$celltype == "Tumor" & data$Vimentin >= 0.5), "celltype_EMT"] <- "EMT"
data$celltype_merged <- data$celltype
data[which(data$celltype == "Immune" | data$celltype == "CD4.T.cells" | data$celltype == "CD8.T.cells" | data$celltype == "FOXP3.CD4.Tregs" | data$celltype == "CD68.Macrophages" |  data$celltype == "CD206.Macrophages" ), "celltype_EMT"] <- "Immune"
#data_stages <- data[-which(data$imageid == "T3" | data$imageid == "M1" | data$imageid == "INF3"),]
data_stages <- data
unique(data_stages$celltype_merged)
for (i in c("Tumor center", "Invasive border", "Inflammation")){
  for (n in c("Tumor", "EMT")){
    freq_spat_expr_tc <- data_stages[which(data_stages$Annotation == i & data_stages$celltype_EMT == n),] %>%
      group_by(patient_id_AB19_1654, celltype_merged) %>% summarize(GCLC = median(GCLC), NQO1 = median(NQO1), 
                                                                    TXNRD1 = median(TXNRD1), Ki67 = median(Ki67_2))
    
    freq_spat_expr_tc <- as.data.frame(freq_spat_expr_tc)
    freq_spat_expr_long <- pivot_longer(freq_spat_expr_tc, cols=GCLC:Ki67, values_to="expr")
    
    clinical <- clinical_data
    #colnames(clinical)[1] <- "patient_id_AB19_1654"
    #colnames(clinical)[12] <- "Stage"
    
    list_cox <- list()
    for (j in unique(freq_spat_expr_long$name)){
      
      colnames(freq_spat_expr_long)[1] <- "patient_id_AB19_1654"
      
      proportions <- freq_spat_expr_long[which(freq_spat_expr_long$name == j),]
      
      
      proportions <- merge(proportions, clinical[, c("patient_id_AB19_1654", "time", "status", "Stage", "rounded_age")], by="patient_id_AB19_1654")
      proportions <- unique(proportions)
      proportions$status <- as.numeric(proportions$status)
      
      tile_number = 5
      proportions <- proportions %>%
        mutate(expr = ntile(expr, tile_number))
      
      #proportions$expr <- as.factor(proportions$expr)
      proportions <- proportions[-which(proportions$Stage == "unknown"),]
      print(paste0(j))
      res.cox <- coxph(Surv(time, status) ~ expr + Stage + tt(rounded_age), data =  proportions, tt=function(x, t, ...) {
        age <- x + t
        cbind(cage=age, cage2= (age-50)^2, cage3= (age-50)^3)
      })
      table1 <- res.cox %>% 
        broom::tidy(exp = TRUE)
      table1$var <- j
      ci <- exp(confint(res.cox))
      ci <- as.data.frame(ci)
      table1 <- cbind(table1, ci)
      list_cox[[j]] <- table1
      print(summary(res.cox))
      #cox.zph(res.cox)
    }
    
    inf <- bind_rows(list_cox)
    
    inf <- inf[which(inf$term == "expr"),]
    write.csv(inf, paste0("/Users/wenqchen/Desktop/Projects/Auria/Clinical/",i,"_func_markers_in_",n,"_cox_regression_table_ntile_",tile_number,".csv"))
    
  }
}


# Forest plot
for (m in c("Tumor center", "Invasive border", "Inflammation")){
  for(n in unique(data_stages$celltype_EMT)[c(2)]){
    forest_data <- read.csv(paste0("/Users/wenqchen/Desktop/Projects/Auria/Clinical/final/",m,"_func_markers_in_",n,"_cox_regression_table_ntile_",tile_number,".csv"))
    colnames(forest_data)[c(8:9)] <- c("low", "high")
    colnames(forest_data)[7] <- "Subgroup"
    forest_data$fdr <- p.adjust(forest_data$p.value)
    forest_data$` ` <- paste(rep(" ", 20), collapse = " ")
    
    tabletext <- rbind(
      c("Marker", "HR", "95% CI", "p-value"),  # header
      cbind(
        forest_data$Subgroup,
        sprintf("%.2f", forest_data$estimate),
        sprintf("[%.2f, %.2f]", forest_data$low, forest_data$high),
        sprintf("%.3f", forest_data$p.value)
      )
    )
    
    colors <- ifelse(forest_data$p.value < 0.05, "royalblue", "royalblue")
    
    p <- forestplot(
      labeltext = tabletext,
      mean  = c(NA, forest_data$estimate),
      lower = c(NA, forest_data$low),
      upper = c(NA, forest_data$high),
      zero = 1,
      xlog = TRUE,
      boxsize = 0.2,
      line.margin = 0.1,
      graph.pos = 2, 
      col = fpColors(box = colors, line = "darkgray", zero = "black"),
      txt_gp = fpTxtGp(
        label = gpar(fontsize = 12),
        ticks = gpar(fontsize = 10),
        xlab  = gpar(fontface = "bold")
      ),
      xlab = "Hazard Ratio (log scale)",
      title = "Forest Plot of Functional Markers"
    )
    
    # Print plot
    pdf(paste0("/Users/wenqchen/Desktop/Projects/Auria/Clinical/final/functional_markers_",m,"_",n,"_new_celltype_calls_ntile_",tile_number,".pdf"), width=6, height=4)
    plot(p)
    dev.off()
  }}

# ---- Single cox ----

df_clinical_GV <- read.csv('/Users/wenqchen/Desktop/Projects/Auria/Data/t-CycIF/clinical_data_GCLC_VIM_quarter.csv')
df_clinical_GV$G_V_Pos <- df_clinical_GV$Tumor.GCLC.VIM..3
ntile = 5
df_clinical_GV <- df_clinical_GV %>% mutate(G_V_Pos = ntile(G_V_Pos, ntile))

df_sub <- subset(df_clinical_GV, subset = G_V_Pos %in% c('1', ntile))

fit <- survfit(Surv(time, status) ~ G_V_Pos, data = df_sub)

surv_plot <- ggsurvplot(
  fit, data = df_sub, 
  size = 1,                 # change line size
  palette = c("#E7B800", "#2E9FDF"),# custom color palettes
  FontSize = 8,
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  # Add risk table
  risk.table = TRUE,
  risk.table.col = 'strata',      # Risk table color by groups
  legend.title = "Tumor_GCLC+VIM+",
  legend.labs = c("High", "Low"),
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_light()      # Change ggplot2 theme
)

# Remove legend from risk table
surv_plot$table <- surv_plot$table + theme(legend.position = "none") + labs(y = NULL)
combined_plot <- surv_plot$plot / surv_plot$table + plot_layout(heights = c(3.3, 0.7))

# 
ggsave(
  filename = paste0("/Users/wenqchen/Desktop/Projects/Auria/Clinical/G_V_Pos_ntile_",ntile,"_survival_curve.pdf"),
  plot = combined_plot,
  width = 6, height = 5, dpi = 300
)

## ---- Multivariate Cox ----
df_clinical_GV <- df_clinical_GV[-which(df_clinical_GV$Stage == "unknown"),]
res.cox <- coxph(Surv(time, status) ~ G_V_Pos + Stage + rounded_age, data = df_clinical_GV)
summary(res.cox)

# Create the new data  
G_V_Pos_df <- with(df_clinical_GV,
               data.frame(G_V_Pos = c('1', ntile),
                          rounded_age = rep(mean(rounded_age, na.rm = TRUE), 2),
                          Stage = c('1', '1')
                          ))

# Survival curves
fit <- survfit(res.cox, newdata = G_V_Pos_df)
ggsurvplot(fit, data = G_V_Pos_df, conf.int = TRUE, pval = TRUE,
           legend.labs=c("High", "Low"),
           ggtheme = theme_bw())

ggsurvplot(
  fit,
  size = 1,                 # change line size
  palette = c("#E7B800", "#2E9FDF"),# custom color palettes
  FontSize = 8,
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  # Add risk table
  risk.table = TRUE,
  risk.table.col = 'strata',      # Risk table color by groups
  legend.title = "Tumor_GCLC+VIM+",
  legend.labs = c("High", "Low"),
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)


res.cox <- coxph(Surv(time, status) ~ G_V_Pos + Stage + tt(rounded_age), data = df_clinical_GV, tt=function(x, t, ...) {
  age <- x + t
  cbind(cage=age, cage2= (age-50)^2, cage3= (age-50)^3)
})
table1 <- res.cox %>% 
  broom::tidy(exp = TRUE)
#table1$var <- j
ci <- exp(confint(res.cox))
ci <- as.data.frame(ci)
table1 <- cbind(table1, ci)
#list_cox[[j]] <- table1
print(summary(res.cox))

forest_data = table1
forest_data <- subset(forest_data, subset = term != 'Stage3')
forest_data <- subset(forest_data, subset = term != 'tt(rounded_age)cage')
colnames(forest_data)[c(6:7)] <- c("low", "high")
colnames(forest_data)[1] <- "Subgroup"
forest_data$fdr <- p.adjust(forest_data$p.value)
forest_data$` ` <- paste(rep(" ", 20), collapse = " ")

tabletext <- rbind(
  c("Subgroup", "HR", "95% CI", "p-value"),  # header
  cbind(
    forest_data$Subgroup,
    sprintf("%.2f", forest_data$estimate),
    sprintf("[%.2f, %.2f]", forest_data$low, forest_data$high),
    sprintf("%.3f", forest_data$p.value)
  )
)

colors <- ifelse(forest_data$p.value < 0.05, "firebrick", "royalblue")


forestplot(
  labeltext = tabletext,
  mean  = c(NA, forest_data$estimate),
  lower = c(NA, forest_data$low),
  upper = c(NA, forest_data$high),
  zero = 1,
  xlog = TRUE,
  boxsize = 0.2,
  line.margin = 0.1,
  graph.pos = 2, 
  col = fpColors(box = colors, line = "darkgray", zero = "black"),
  txt_gp = fpTxtGp(
    label = gpar(fontsize = 12),
    ticks = gpar(fontsize = 10),
    xlab  = gpar(fontface = "bold")
  ),
  xlab = "Hazard Ratio (log scale)",
  title = "Forest Plot of Functional Markers"
)




















