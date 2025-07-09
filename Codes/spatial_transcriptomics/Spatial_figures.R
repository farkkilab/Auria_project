
### spatial plots AURIA project
# Author(s): Sara Palomino, Wenqing Chen

# ---- Load libraries ----
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("/Users/spalominoe/Python_jobs/01_scripts/figures/")

# ---- Color list ----
color_vector <- c("#0072B2", "#B3CDE3","#FBC15E","#8DD3C7",  "#FB8072", "#CAB2D6", "#A6611A",  "#F781BF","#999999",       
                  "#BDBF00","#00CED1", "#A6CEE3"  )

######### sample_CID44971 ######### 
CID44971 <- read.csv("CID44971_spatial.csv")


CID44971_aggregated_data <- as.data.frame(t(CID44971 %>%
                                              group_by(oxstress_yn) %>%
                                              summarise(across(ends_with('total'), mean))))

cell_types <- gsub("_total$", "", rownames(CID44971_aggregated_data))
rownames(CID44971_aggregated_data) <- cell_types


CID44971_df <- data.frame(percent=c(CID44971_aggregated_data[-1,1],CID44971_aggregated_data[-1,2]))
CID44971_df$Sample <- c(rep("no oxstress",12), rep("oxtress", 12))
CID44971_df$cell_type <- c(rownames(CID44971_aggregated_data)[-1],rownames(CID44971_aggregated_data) [-1])


CID44971 <- as.data.frame(CID44971[,-c(13:16)])

# Reshape the data to long format
CID44971_data_long <- CID44971 %>%
  pivot_longer(cols = -c(oxstress_yn,GCLCVIM_yn, oxstress, GCLCVIM) , names_to = "cell_type", values_to = "percent")


#### OXSTRESS
# Compute the log fold change (logFC) for each cell type
oxstress_logFC_results <- CID44971_data_long %>%
  group_by(cell_type, oxstress_yn) %>%
  summarise(mean_percent = mean(percent)) %>%
  spread(key = oxstress_yn, value = mean_percent) %>%
  mutate(logFC = log2(`1` / `0`)) %>%
  select(cell_type, logFC)

# Perform t-test (or Wilcoxon test if data is not normally distributed)
stat_test_results_oxstsress <- CID44971_data_long %>%
  group_by(cell_type) %>%
  summarise(p_value = wilcox.test(percent ~ oxstress_yn)$p.value)

# Now add the adjusted p-value to the results
stat_test_results_oxstsress$adjusted_p_value <- p.adjust(stat_test_results_oxstsress$p_value, method = "BH")

# Combine logFC and adjusted p-values into a final results dataframe
final_results_oxstsress <- left_join(oxstress_logFC_results, stat_test_results_oxstsress, by = "cell_type")

CID44971_final_results_oxstsress <- final_results_oxstsress

CID44971_final_results_oxstsress$slide <- "Slide 1"


######### sample_CID4465 ######### 
# Reshape the data to long format
CID4465_data_long <- CID4465 %>%
  pivot_longer(cols = -c(oxstress_yn,GCLCVIM_yn, oxstress, GCLCVIM) , names_to = "cell_type", values_to = "percent")


#### OXSTRESS
# Compute the log fold change (logFC) for each cell type
oxstress_logFC_results <- CID4465_data_long %>%
  group_by(cell_type, oxstress_yn) %>%
  summarise(mean_percent = mean(percent)) %>%
  spread(key = oxstress_yn, value = mean_percent) %>%
  mutate(logFC = log2(`1` / `0`)) %>%
  select(cell_type, logFC)

# Perform t-test (or Wilcoxon test if data is not normally distributed)
stat_test_results_oxstsress <- CID4465_data_long %>%
  group_by(cell_type) %>%
  summarise(p_value = wilcox.test(percent ~ oxstress_yn)$p.value)


# Now add the adjusted p-value to the results
stat_test_results_oxstsress$adjusted_p_value <- p.adjust(stat_test_results_oxstsress$p_value, method = "BH")

# Combine logFC and adjusted p-values into a final results dataframe
final_results_oxstsress <- left_join(oxstress_logFC_results, stat_test_results_oxstsress, by = "cell_type")

CID4465_final_results_oxstsress <- final_results_oxstsress
CID4465_final_results_oxstsress$slide <- "Slide 2"



### FIGURES ####

# > Barplot
CID44971_df$slide <- "slide 1"
CID44971_df$status <- "slide 1 oxstress"
CID44971_df$status[CID44971_df$Sample=="no oxstress"] <- "slide 1 no oxstress"

CID4465_df$slide <- "slide 2"
CID4465_df$status <- "slide 2 oxstress"
CID4465_df$status[CID4465_df$Sample=="no oxstress"] <- "slide 2 no oxstress"
CID4465_df$percent <- as.numeric(CID4465_df$percent)

slides_1_2 <- rbind(CID44971_df,CID4465_df)

slides_1_2$status <- factor(slides_1_2$status, levels = c("slide 2 no oxstress",    "slide 2 oxstress",
                                                          "slide 1 no oxstress",    "slide 1 oxstress" ))


barplot_horizontal <- ggplot(slides_1_2, aes(x = percent, y = status, fill = factor(cell_type))) + 
  geom_bar(stat = "identity") + 
  geom_col(position = position_stack(reverse = TRUE)) + 
  scale_fill_manual(values = color_vector) +  # Corrected here to pass color_palette to 'values'
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(colour = "black", size = 10),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.position = "bottom",
    panel.background = element_rect(fill = 'transparent'),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  labs(x = "Cell counts (%)", y = "", fill = "Cell Type")

pdf("barplot_horizontal.pdf", height = 3)
barplot_horizontal
dev.off()

# > Cleveland
CID4465_final_results_oxstsress
CID44971_final_results_oxstsress

stats <- rbind(CID44971_final_results_oxstsress, CID4465_final_results_oxstsress)
stats$cell_type <- gsub("_total$", "", stats$cell_type)
stats$cell_type[stats$cell_type=="Myoepithelia"] <- "Myoepithelial"
stats$cell_type[stats$cell_type=="B_cells"] <- "B cells"

stats$cell_type <- factor(stats$cell_type, levels = rev(unique(stats$cell_type)))


stats$colors <- "#f183a1" # slide 1
stats$colors[stats$slide=="Slide 2"] <- "#6a9cd4"


cleveland_oxstress <- ggplot(stats, aes(logFC, (cell_type), color=slide)) +
  ggtitle("Cell type %: Oxstress hotspots vs non oxstress") +
  geom_point(aes(size=adjusted_p_value),alpha =0.7) +
  guides(size = guide_legend(reverse=F))+
  scale_size("P-value (FDR)",range=c(4,1), breaks=c(0.01,0.05,0.1, 0.5)) +
  scale_colour_manual(values=setNames(stats$colors,stats$slide))+
  xlim(-4, 4)+ geom_vline(xintercept = c(-0.58,0.58), linetype="dotted", color = "black", size=0.5)+
  labs(x = "Log2FC", y = " ", fill = "") + 
  theme(
    axis.title =element_text(size=11),
    axis.text.x=element_text(angle=45, hjust=1,size=10),
    axis.text.y=element_text(size=11),
    plot.title = element_text(face = "bold", size=12),
    legend.text =   element_text(size=10), legend.key.size = unit(0.4, "cm"), legend.title = element_text(size=11))

cleveland_oxstress

setwd("/Users/spalominoe/Desktop/AURIA FIGURES/")
pdf("cleveland_oxstress.pdf", height = 5, width = 6)
cleveland_oxstress
dev.off()


# > corplot
names(CID44971)
CID44971_tocor <- CID44971[,-c(13,14,16)]
names(CID44971_tocor)    
str(CID44971_tocor)
names(CID44971_tocor) <- gsub("_total$", "", names(CID44971_tocor))

CID44971_tocor[is.na(CID44971_tocor)] <- "0"
CID44971_tocor$oxstress <- as.numeric(CID44971_tocor$oxstress)
str(CID44971_tocor)
names(CID44971_tocor) <- gsub("_total$", "", names(CID44971_tocor))

M = cor(CID44971_tocor)
stats_cor <- cor.mtest(M)
p_values <- stats_cor$p

corrplot(as.matrix(M[,13])) # colorful number
corrplot(as.matrix(M[,13]),  col=rev(COL2('RdBu', 200)),  diag = FALSE)


# Define the significance level
significance_level <- 0.05

# Add asterisks to the correlations that are statistically significant (p-value < 0.05)
significance_matrix <- ifelse(p_values < significance_level, "*", "")

prueba <- as.matrix(M[,13])

# Plot the correlation matrix with asterisks
CID44971_tocor <- corrplot(prueba, col=rev(COL2('RdBu', 200)),
                           method = "circle",
                           sig.level = significance_level,  # Define the significance level
                           insig = "label_sig", type = 'lower', diag = FALSE, title="Slide 1")  # Show asterisk for significant correlations)  # Add coefficients to the plot)  # Adjust number size for correlation coefficients





names(CID4465)
CID4465_tocor <- CID4465[,-c(13:18,20)]
names(CID4465_tocor)    
str(CID4465_tocor)
names(CID4465_tocor) <- gsub("_total$", "", names(CID4465_tocor))

CID4465_tocor[is.na(CID4465_tocor)] <- "0"
CID4465_tocor$oxstress <- as.numeric(CID4465_tocor$oxstress)
str(CID4465_tocor)
names(CID4465_tocor) <- gsub("_total$", "", names(CID4465_tocor))


M = cor(CID4465_tocor)
stats_cor <- cor.mtest(M)
p_values <- stats_cor$p

corrplot(as.matrix(M[,13])) # colorful number
corrplot(as.matrix(M[,13]),  col=rev(COL2('RdBu', 200)),  diag = FALSE)



agg_tocor <- rbind(CID44971_tocor,CID4465_tocor)


M = cor(agg_tocor)
stats_cor <- cor.mtest(M)
p_values <- stats_cor$p

corrplot(M) # colorful number
corrplot(M[], order = 'AOE', col=rev(COL2('RdBu', 200)))


# Define the significance level
significance_level <- 0.05

# Add asterisks to the correlations that are statistically significant (p-value < 0.05)
significance_matrix <- ifelse(p_values < significance_level, "*", "")


