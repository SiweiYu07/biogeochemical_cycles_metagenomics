
library(ggplot2)
library(reshape2)
library(grid)
library(tidyverse)
library(psych)
library(pheatmap)

#import data files
env <- read.csv("LG_env.csv", check.names = FALSE)
sp <- read.csv("LG_TPM.csv", check.names = FALSE)

#Initializes two matrices, one for storing correlation values 
# and the other for storing significance markers
correlation_matrix <- matrix(NA, nrow = ncol(env), ncol = ncol(sp))
significance_matrix <- matrix("", nrow = ncol(env), ncol = ncol(sp))
rownames(correlation_matrix) <- colnames(env)
colnames(correlation_matrix) <- colnames(sp)

#calculate correlations and significance, add custom significance markers
for (i in 1:ncol(env)) {
  for (j in 1:ncol(sp)) {
    result <- cor.test(env[,i], sp[,j], method = "spearman", use = "complete.obs")
    correlation_matrix[i, j] <- result$estimate
    p_value <- result$p.value
    if (p_value < 0.001) {
      significance_matrix[i, j] <- "***"
    } else if (p_value < 0.01) {
      significance_matrix[i, j] <- "**"
    } else if (p_value < 0.05) {
      significance_matrix[i, j] <- "*"
    }
    
  }
}

melted_correlation_matrix <- melt(correlation_matrix)
colnames(melted_correlation_matrix) <- c("Var1", "Var2", "Correlation")
melted_correlation_matrix$Correlation <- as.numeric(as.character(melted_correlation_matrix$Correlation))

#significance matrix convert into data frame incase the data type were correct
melted_significance_matrix <- melt(significance_matrix)
colnames(melted_significance_matrix) <- c("Var1", "Var2", "Significance")
melted_significance_matrix$Significance <- as.character(melted_significance_matrix$Significance)

#ensure the colnames were consistency
colnames(melted_correlation_matrix) <- c("Var1", "Var2", "Correlation")
colnames(melted_significance_matrix) <- c("Var1", "Var2", "Significance")

#transfer data type
melted_correlation_matrix[,1:2] <- lapply(melted_correlation_matrix[,1:2],as.character)
melted_significance_matrix[,1:2] <- lapply(melted_significance_matrix[,1:2],as.character)

#merge correlations and significance data frame
merged_data <- cbind(melted_correlation_matrix, melted_significance_matrix)
merged_data <- merged_data[,-c(4,5)]%>%na.omit()

#plot heatmap laving top space for custom text
heatmap_plot <- ggplot(merged_data, aes(Var1, Var2, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = Significance), color="black",size = 12) +
  scale_fill_gradient2(low = "#999999" , high = "#0066CC", mid = "white", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(hjust = 1, size = 20)
  ) 
#print heatmap
print(heatmap_plot)

#add the custom text
grid.text(" *: P < 0.05\n **: P < 0.01\n ***: P < 0.001",
          x = 0.88, y = 0.4, just = "left", gp = gpar(col = "black", fontsize = 10))
ggsave("LG_env_TPM.pdf", heatmap_plot, width = 10, height = 10)
ggsave("LG_env_TPM.png", heatmap_plot, width = 10, height = 10, dpi = 600)