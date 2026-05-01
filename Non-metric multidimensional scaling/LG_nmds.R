# Load packages
library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)

# 1. Data preparation
ko_data <- read.csv("LG_nmds.csv", row.names = 1, check.names = FALSE)

# Ensure that the Sample column is read as character when loading group data
group_data <- read.csv("LG_group.csv", header = TRUE) %>%
  mutate(Sample = as.character(Sample))  # Key modification

# 2. Calculate relative abundance
ko_rel <- t(t(ko_data) / colSums(ko_data)) * 100

# 3. Filter low-abundance KOs
ko_filtered <- ko_rel[rowSums(ko_rel > 0.1) >= 5, ]

# 4. NMDS analysis
set.seed(123)
nmds_result <- metaMDS(t(ko_filtered), distance = "bray", k = 2)

# 5. Extract coordinates and merge with group data; data types now match
nmds_points <- scores(nmds_result, "sites") %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  left_join(group_data, by = "Sample")

# 6. Visualization
ggplot(nmds_points, aes(x = NMDS1, y = NMDS2, color = Group)) +
  geom_point(size = 4) +
  geom_text_repel(
    aes(label = Sample), 
    box.padding = 0.5,
    show.legend = FALSE,
    size = 5  # Increase sample label font size
  ) +
  stat_ellipse(
    aes(group = Group),
    level = 0.8,
    linetype = 2,
    linewidth = 0.8  # Adjust ellipse line width
  ) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "NMDS of KO Relative Abundance by Group",
    subtitle = paste("Stress =", round(nmds_result$stress, 3)),
    color = "Experimental Group"
  ) +
  theme_minimal() +
  theme(
    # Adjust axis and plot text sizes
    axis.title = element_text(size = 16, face = "bold"),  # Axis titles
    axis.text = element_text(size = 16),  # Axis tick labels
    plot.title = element_text(size = 16, face = "bold"),  # Main title
    plot.subtitle = element_text(size = 16),  # Subtitle
    legend.title = element_text(size = 16),  # Legend title
    legend.text = element_text(size = 16)  # Legend text
  )

# PERMANOVA test for between-group differences
adonis2(t(ko_filtered) ~ Group, data = group_data, permutations = 999)

# Similarity analysis using ANOSIM
anosim(t(ko_filtered), group_data$Group)

ggsave("LG_NMDS_with_Groups.pdf", width = 8, height = 6)
write.csv(nmds_points, "LG_NMDS_coordinates_with_groups.csv")

# Plot a boxplot of between-group distances
library(ggpubr)
dist_matrix <- vegdist(t(ko_filtered))
group_df <- data.frame(
  Distance = as.vector(dist_matrix),
  Comparison = outer(group_data$Group, group_data$Group, 
                     function(x,y) paste(x,"vs",y))[lower.tri(dist_matrix)]
)
ggboxplot(group_df, x = "Comparison", y = "Distance") +
  stat_compare_means()


--------------------------------------------------------------------
  # Load packages
  library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(reshape2)  # Used for data reshaping

# 1. Data preparation
ko_data <- read.csv("LG_nmds.csv", row.names = 1, check.names = FALSE)

# Ensure that the Sample column is read as character when loading group data
group_data <- read.csv("LG_group.csv", header = TRUE) %>%
  mutate(Sample = as.character(Sample))

# 2. Calculate relative abundance
ko_rel <- t(t(ko_data) / colSums(ko_data)) * 100

# 3. Filter low-abundance KOs
ko_filtered <- ko_rel[rowSums(ko_rel > 0.1) >= 5, ]

# 4. Calculate the Bray-Curtis distance matrix, the core metric of beta diversity
set.seed(123)
bray_dist <- vegdist(t(ko_filtered), method = "bray", diag = TRUE, upper = TRUE)
bray_matrix <- as.matrix(bray_dist)

# 5. Save the Bray-Curtis distance matrix
write.csv(bray_matrix, "LG_BrayCurtis_DistanceMatrix.csv")

# 6. Output a statistical summary of the distance matrix
cat("\n===== Statistical summary of the Bray-Curtis distance matrix =====\n")
cat("Number of samples:", ncol(ko_filtered), "\n")
cat("Distance matrix dimensions:", dim(bray_matrix), "\n")
cat("Minimum distance:", min(bray_matrix[upper.tri(bray_matrix)]), "\n")
cat("Maximum distance:", max(bray_matrix[upper.tri(bray_matrix)]), "\n")
cat("Mean distance:", mean(bray_matrix[upper.tri(bray_matrix)]), "\n")
cat("Median distance:", median(bray_matrix[upper.tri(bray_matrix)]), "\n")
cat("Distance standard deviation:", sd(bray_matrix[upper.tri(bray_matrix)]), "\n\n")

# 7. Calculate Bray-Curtis-based diversity metrics for each sample
sample_diversity <- data.frame(
  Sample = rownames(bray_matrix),
  # Mean distance: average dissimilarity between this sample and all other samples
  Mean_Distance = apply(bray_matrix, 1, function(x) mean(x[x > 0])),
  # Minimum distance: nearest neighbor
  Min_Distance = apply(bray_matrix, 1, function(x) min(x[x > 0])),
  # Maximum distance: farthest neighbor
  Max_Distance = apply(bray_matrix, 1, function(x) max(x[x > 0])),
  # Distance standard deviation: variation in dissimilarity between this sample and other samples
  SD_Distance = apply(bray_matrix, 1, function(x) sd(x[x > 0]))
)

# Add grouping information
sample_diversity <- sample_diversity %>%
  left_join(group_data, by = "Sample")

# 8. Save diversity metrics for each sample
write.csv(sample_diversity, "LG_Sample_BrayCurtis_Diversity.csv", row.names = FALSE)

# 9. Output sample-level diversity statistics
cat("===== Bray-Curtis diversity metrics for each sample =====\n")
print(sample_diversity)
cat("\n")

# 10. Visualization: mean Bray-Curtis distance for each sample
avg_dist_plot <- ggplot(sample_diversity, 
                        aes(x = reorder(Sample, Mean_Distance), 
                            y = Mean_Distance, 
                            fill = Group)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Mean_Distance - SD_Distance, 
                    ymax = Mean_Distance + SD_Distance), 
                width = 0.2, alpha = 0.7) +
  geom_text(aes(label = round(Mean_Distance, 3)), 
            vjust = -0.5, size = 3) +
  labs(
    title = "Mean Bray-Curtis Distance for Each Sample",
    subtitle = "Error bars indicate standard deviation",
    x = "Sample",
    y = "Mean Bray-Curtis distance",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  ) +
  scale_fill_brewer(palette = "Set1")

ggsave("LG_Sample_MeanBrayCurtis_Distance.pdf", avg_dist_plot, width = 12, height = 6)

# 11. Visualization: heatmap of the Bray-Curtis distance matrix
heatmap_plot <- ggplot(melt(bray_matrix), aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = median(bray_matrix[upper.tri(bray_matrix)]),
    name = "Bray-Curtis\ndistance"
  ) +
  labs(
    title = "Heatmap of the Bray-Curtis Distance Matrix",
    x = "Sample", 
    y = "Sample"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    legend.position = "right"
  )

ggsave("LG_BrayCurtis_DistanceHeatmap.pdf", heatmap_plot, width = 10, height = 8)

# 12. NMDS analysis based on Bray-Curtis distance
cat("===== NMDS analysis based on Bray-Curtis distance =====\n")
nmds_result <- metaMDS(t(ko_filtered), distance = "bray", k = 2)
cat("Stress value:", round(nmds_result$stress, 3), "\n\n")

# 13. Extract coordinates and merge with sample diversity data
nmds_points <- scores(nmds_result, "sites") %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  left_join(sample_diversity, by = "Sample")

# 14. Visualization: NMDS plot colored by mean Bray-Curtis distance
nmds_bray_plot <- ggplot(nmds_points, 
                         aes(x = NMDS1, y = NMDS2, 
                             color = Mean_Distance, 
                             shape = Group)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text_repel(
    aes(label = Sample), 
    box.padding = 0.5,
    show.legend = FALSE,
    size = 3
  ) +
  scale_color_gradient2(
    low = "blue", mid = "yellow", high = "red",
    midpoint = median(nmds_points$Mean_Distance),
    name = "Mean Bray-Curtis\ndistance"
  ) +
  scale_shape_manual(values = c(16, 17, 18, 15, 19)) +  # Different shapes indicate different groups
  labs(
    title = "NMDS Plot Colored by Mean Bray-Curtis Distance",
    subtitle = paste("Based on Bray-Curtis distance | Stress =", round(nmds_result$stress, 3))
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  )

ggsave("LG_NMDS_with_BrayCurtis_Color.pdf", nmds_bray_plot, width = 10, height = 8)

# 15. Correlation analysis between NMDS coordinates and Bray-Curtis distance
cat("===== Correlation analysis between NMDS coordinates and Bray-Curtis distance =====\n")

# Calculate NMDS distances
nmds_coords <- as.matrix(nmds_points[, c("NMDS1", "NMDS2")])
nmds_dist <- dist(nmds_coords)

# Mantel test: correlation between NMDS distance and Bray-Curtis distance
mantel_test <- mantel(bray_dist, nmds_dist, method = "pearson", permutations = 999)
cat("Mantel test result:\n")
cat("Correlation coefficient r =", round(mantel_test$statistic, 4), "\n")
cat("p-value =", mantel_test$signif, "\n\n")

# 16. Multidimensional scaling analysis to directly visualize Bray-Curtis distance
mds_result <- cmdscale(bray_dist, k = 2, eig = TRUE)
mds_points <- as.data.frame(mds_result$points)
colnames(mds_points) <- c("MDS1", "MDS2")
mds_points$Sample <- rownames(mds_points)
mds_points <- mds_points %>%
  left_join(sample_diversity, by = "Sample")

# Calculate percentage of explained variance
var_explained <- round(100 * mds_result$eig / sum(mds_result$eig), 2)

mds_plot <- ggplot(mds_points, 
                   aes(x = MDS1, y = MDS2, 
                       color = Mean_Distance, 
                       shape = Group)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text_repel(
    aes(label = Sample), 
    box.padding = 0.5,
    show.legend = FALSE,
    size = 3
  ) +
  scale_color_gradient2(
    low = "blue", mid = "yellow", high = "red",
    midpoint = median(mds_points$Mean_Distance),
    name = "Mean Bray-Curtis\ndistance"
  ) +
  scale_shape_manual(values = c(16, 17, 18, 15, 19)) +
  labs(
    title = "Classical MDS Analysis Based on Bray-Curtis Distance",
    subtitle = paste("Variance explained by Dimension 1:", var_explained[1], "% | Dimension 2:", var_explained[2], "%"),
    x = paste("MDS1 (", var_explained[1], "%)", sep = ""),
    y = paste("MDS2 (", var_explained[2], "%)", sep = "")
  ) +
  theme_minimal()

ggsave("LG_Classical_MDS_BrayCurtis.pdf", mds_plot, width = 10, height = 8)

# 17. PERMANOVA test based on Bray-Curtis distance
cat("===== PERMANOVA test based on Bray-Curtis distance =====\n")
permanova_result <- adonis2(bray_dist ~ Group, data = group_data, permutations = 999)
print(permanova_result)
cat("\n")

# 18. ANOSIM test based on Bray-Curtis distance
cat("===== ANOSIM test based on Bray-Curtis distance =====\n")
anosim_result <- anosim(t(ko_filtered), group_data$Group, distance = "bray")
print(anosim_result)
cat("\n")

# 19. Summary output
cat("===== Bray-Curtis diversity analysis completed =====\n")
cat("The following files have been generated:\n")
cat("1. LG_BrayCurtis_DistanceMatrix.csv - Complete Bray-Curtis distance matrix\n")
cat("2. LG_Sample_BrayCurtis_Diversity.csv - Bray-Curtis diversity metrics for each sample\n")
cat("3. LG_Sample_MeanBrayCurtis_Distance.pdf - Bar plot of the mean distance for each sample\n")
cat("4. LG_BrayCurtis_DistanceHeatmap.pdf - Bray-Curtis distance heatmap\n")
cat("5. LG_NMDS_with_BrayCurtis_Color.pdf - NMDS plot colored by Bray-Curtis distance\n")
cat("6. LG_Classical_MDS_BrayCurtis.pdf - Classical MDS analysis plot\n")
cat("\nMain analysis results:\n")
cat("1. Bray-Curtis distance matrix statistics: mean distance =", round(mean(bray_matrix[upper.tri(bray_matrix)]), 4), "\n")
cat("2. Mantel test: correlation coefficient between NMDS distance and original Bray-Curtis distance =", round(mantel_test$statistic, 4), "\n")
cat("3. PERMANOVA: significance of between-group differences (p =", permanova_result$`Pr(>F)`[1], ")\n")
cat("4. ANOSIM: significance of between-group differences (R =", anosim_result$statistic, ", p =", anosim_result$signif, ")\n")