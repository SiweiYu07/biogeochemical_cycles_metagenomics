# ==== 0. Load required packages ====
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(Polychrome)

# ==== 1. Read data ====
gene_data <- read.csv("LG_contribution.csv", header = TRUE, check.names = FALSE)
group_data <- read.csv("LG_group.csv", header = TRUE)

# ==== 2. Data preprocessing ====
sample_names <- colnames(gene_data)[5:ncol(gene_data)]

gene_long <- gene_data %>%
  select(Element.cycling, Function.category, KO, Gene, all_of(sample_names)) %>%
  pivot_longer(cols = all_of(sample_names), names_to = "Bins", values_to = "TPM")

merged_data <- merge(gene_long, group_data, by = "Bins")

# ==== 3. Calculate mean TPM table ====
function_phase_data <- merged_data %>%
  group_by(Element.cycling, Function.category, KO, Gene, Phase) %>%
  summarise(Mean_TPM = mean(TPM, na.rm = TRUE), .groups = "drop")

# ==== 4. Calculate log2FC ====
phases <- c("LG-4", "LG-3", "LG-2", "LG-1")

function_wide <- function_phase_data %>%
  pivot_wider(names_from = Phase, values_from = Mean_TPM, values_fill = 0)

log2fc_data <- function_wide %>%
  select(Element.cycling, Function.category, KO, Gene)

for (phase in phases[-1]) {
  col <- paste0(phase, "/LG-4")
  log2fc_data[[col]] <- log2((function_wide[[phase]] + 1) / (function_wide[["LG-4"]] + 1))
}

# ==== 5. Significance testing ====
pvalue_data <- log2fc_data[, 1:4]

for (phase in phases[-1]) {
  pval_col <- paste0(phase, "_pvalue")
  pvalue_data[[pval_col]] <- sapply(1:nrow(log2fc_data), function(i) {
    row <- log2fc_data[i, ]
    gene_expr <- merged_data %>%
      filter(Element.cycling == row$Element.cycling,
             Function.category == row$Function.category,
             KO == row$KO,
             Gene == row$Gene)
    g1 <- gene_expr %>% filter(Phase == "LG-4") %>% pull(TPM)
    g2 <- gene_expr %>% filter(Phase == phase) %>% pull(TPM)
    if (length(g1) < 2 || length(g2) < 2) return(NA)
    wilcox.test(g1, g2, exact = FALSE)$p.value
  })
}

# ==== 6. Build heatmap matrix ====
comparison_cols <- colnames(log2fc_data)[5:ncol(log2fc_data)]
heatmap_matrix <- as.matrix(log2fc_data[, comparison_cols])
rownames(heatmap_matrix) <- log2fc_data$Gene

star_matrix <- matrix("", nrow = nrow(heatmap_matrix), ncol = length(comparison_cols))
colnames(star_matrix) <- comparison_cols
rownames(star_matrix) <- log2fc_data$Gene

for (i in 1:length(comparison_cols)) {
  phase <- strsplit(comparison_cols[i], "/")[[1]][1]
  pval_col <- paste0(phase, "_pvalue")
  if (pval_col %in% colnames(pvalue_data)) {
    pval <- pvalue_data[[pval_col]]
    star_matrix[, i] <- ifelse(pval < 0.001, "***",
                               ifelse(pval < 0.01, "**",
                                      ifelse(pval < 0.05, "*", "")))
  }
}

# ==== 7. Select the top 5 functional genes in each category ====
fc_threshold <- 1
fc_mask <- apply(abs(heatmap_matrix), 1, function(x) any(x > fc_threshold, na.rm = TRUE))
sig_mask <- apply(star_matrix, 1, function(x) any(x != ""))
selected_mask <- fc_mask | sig_mask

heatmap_sig <- heatmap_matrix[selected_mask, , drop = FALSE]
star_sig <- star_matrix[selected_mask, , drop = FALSE]
max_fc <- apply(abs(heatmap_sig), 1, max)

annot <- log2fc_data %>%
  filter(Gene %in% rownames(heatmap_sig)) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  mutate(maxFC = max_fc[Gene])

top_genes_df <- annot %>%
  group_by(Function.category) %>%
  arrange(desc(maxFC)) %>%
  slice_head(n = 5) %>%
  ungroup()

selected_genes <- top_genes_df$Gene
heatmap_matrix_filtered <- heatmap_sig[selected_genes, , drop = FALSE]
star_matrix_filtered <- star_sig[selected_genes, , drop = FALSE]

# ==== 8. Build annotations ====
annotation_subset <- top_genes_df %>% arrange(Function.category)
sorted_rows <- annotation_subset$Gene
heatmap_matrix_filtered <- heatmap_matrix_filtered[sorted_rows, ]
star_matrix_filtered <- star_matrix_filtered[sorted_rows, ]

row_annotation <- data.frame(Function = annotation_subset$Function.category)
rownames(row_annotation) <- sorted_rows

function_list <- unique(row_annotation$Function)

high_contrast_colors <- createPalette(21, seedcolors = c("#CCCCFF", "#339933"))
names(high_contrast_colors) <- function_list
anno_colors <- list(
  Function = high_contrast_colors
)

# ==== 9. Build display matrix ====
display_matrix <- matrix(
  paste0(round(heatmap_matrix_filtered, 2), " ", star_matrix_filtered),
  nrow = nrow(heatmap_matrix_filtered),
  ncol = ncol(heatmap_matrix_filtered)
)

# ==== 10. Plot heatmap ====
pheatmap(heatmap_matrix_filtered,
         color = colorRampPalette(c("#339933", "white", "purple"))(100),
         breaks = seq(-max(abs(heatmap_matrix_filtered)), max(abs(heatmap_matrix_filtered)), length.out = 100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_row = row_annotation,
         annotation_colors = anno_colors,
         annotation_names_row = TRUE,
         annotation_legend = TRUE,
         legend_breaks = c(-4, -2, 0, 2, 4),
         main = "Functional Gene Expression log2FC (vs LG-4)",
         display_numbers = display_matrix,
         number_color = "black",
         fontsize_number = 28,
         fontsize = 28,
         fontsize_row = 38,
         fontsize_col = 38,
         cellheight = 35,
         cellwidth = 85,
         annotation_legend_side = "right",
         annotation_legend_cols = 1,
         angle_col = 45
)

# ==== 11. Export CSV results ====
write.csv(log2fc_data, "log2fc_results.csv", row.names = FALSE)
write.csv(pvalue_data, "pvalue_results.csv", row.names = FALSE)
