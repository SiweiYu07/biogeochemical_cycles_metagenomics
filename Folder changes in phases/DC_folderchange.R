# ==== 0. Load required packages ====
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(Polychrome)

# ==== 1. Read data ====
gene_data <- read.csv("DC_contribution.csv", header = TRUE, check.names = FALSE)
group_data <- read.csv("DC_group.csv", header = TRUE)

# Requirement: group_data must contain at least the columns Bins and Phase
# Bins must match the sample column names in gene_data
stopifnot(all(c("Bins","Phase") %in% colnames(group_data)))

# ==== 2. Select only numeric or numeric-like sample columns ====
meta_cols <- 4  # The first four columns are: Element.cycling, Function.category, KO, Gene
all_sample_cols <- colnames(gene_data)[(meta_cols + 1):ncol(gene_data)]

to_numeric_if_numeric_like <- function(x) {
  if (is.numeric(x)) return(x)
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    ok <- grepl("^\\s*-?\\d+(?:\\.\\d+)?\\s*$|^\\s*NA\\s*$", x)
    if (mean(ok | is.na(x)) > 0.9) return(suppressWarnings(as.numeric(x)))
  }
  x
}

gene_data[all_sample_cols] <- lapply(gene_data[all_sample_cols], to_numeric_if_numeric_like)
numeric_sample_cols <- names(Filter(is.numeric, gene_data[all_sample_cols]))
if (length(numeric_sample_cols) == 0) stop("No numeric sample columns were detected. Please check the data.")

# ==== 3. Convert to long format and merge with grouping information ====
gene_long <- gene_data %>%
  select(Element.cycling, Function.category, KO, Gene, all_of(numeric_sample_cols)) %>%
  pivot_longer(cols = all_of(numeric_sample_cols), names_to = "Bins", values_to = "TPM")

# Keep only bins with grouping information
merged_data <- gene_long %>% inner_join(group_data, by = "Bins")

# ==== 4. Calculate mean TPM for each phase ====
function_phase_data <- merged_data %>%
  group_by(Element.cycling, Function.category, KO, Gene, Phase) %>%
  summarise(Mean_TPM = mean(TPM, na.rm = TRUE), .groups = "drop")

# ==== 5. Determine phase order and reference phase, prioritizing DC-4 ====
phase_levels <- unique(function_phase_data$Phase)

# Try to sort by the trailing number, such as DC-1, DC-2, ...
phase_nums <- suppressWarnings(as.numeric(sub(".*-(\\d+)$", "\\1", phase_levels)))
if (all(!is.na(phase_nums))) {
  phases <- phase_levels[order(phase_nums, decreasing = TRUE)]  # Order as 4, 3, 2, 1
} else {
  phases <- sort(phase_levels, decreasing = TRUE)
}

# Use DC-4 as the reference phase; if absent, use phases[1]
ref_phase <- if ("DC-4" %in% phases) "DC-4" else phases[1]
comp_phases <- setdiff(phases, ref_phase)

# ==== 6. Calculate log2FC relative to ref_phase ====
function_wide <- function_phase_data %>%
  pivot_wider(names_from = Phase, values_from = Mean_TPM, values_fill = 0)

log2fc_data <- function_wide %>%
  select(Element.cycling, Function.category, KO, Gene)

for (phase in comp_phases) {
  col <- paste0(phase, "/", ref_phase)
  log2fc_data[[col]] <- log2((function_wide[[phase]] + 1) / (function_wide[[ref_phase]] + 1))
}

# ==== 7. Wilcoxon significance test based on bin-level distributions ====
pvalue_data <- log2fc_data[, 1:4]
for (phase in comp_phases) {
  pval_col <- paste0(phase, "_pvalue")
  pvalue_data[[pval_col]] <- vapply(1:nrow(log2fc_data), function(i) {
    row <- log2fc_data[i, ]
    gene_expr <- merged_data %>%
      filter(Element.cycling == row$Element.cycling,
             Function.category == row$Function.category,
             KO == row$KO,
             Gene == row$Gene)
    g_ref <- gene_expr %>% filter(Phase == ref_phase) %>% pull(TPM)
    g_cmp <- gene_expr %>% filter(Phase == phase) %>% pull(TPM)
    if (length(na.omit(g_ref)) < 2 || length(na.omit(g_cmp)) < 2) return(NA_real_)
    suppressWarnings(wilcox.test(g_ref, g_cmp, exact = FALSE)$p.value)
  }, numeric(1))
}

# ==== 8. Build heatmap and significance-star matrices ====
comparison_cols <- colnames(log2fc_data)[5:ncol(log2fc_data)]
heatmap_matrix <- as.matrix(log2fc_data[, comparison_cols])
rownames(heatmap_matrix) <- log2fc_data$Gene

star_matrix <- matrix("", nrow = nrow(heatmap_matrix), ncol = length(comparison_cols))
colnames(star_matrix) <- comparison_cols
rownames(star_matrix) <- rownames(heatmap_matrix)

for (i in seq_along(comparison_cols)) {
  phase_i <- strsplit(comparison_cols[i], "/")[[1]][1]
  pcol <- paste0(phase_i, "_pvalue")
  if (pcol %in% colnames(pvalue_data)) {
    pval <- pvalue_data[[pcol]]
    star_matrix[, i] <- ifelse(pval < 0.001, "***",
                               ifelse(pval < 0.01,  "**",
                                      ifelse(pval < 0.05,  "*", "")))
  }
}

# ==== 9. Select the top 5 genes in each Function.category by |maximum log2FC| ====
fc_threshold <- 1
fc_mask  <- apply(abs(heatmap_matrix), 1, function(x) any(x > fc_threshold, na.rm = TRUE))
sig_mask <- apply(star_matrix,   1, function(x) any(x != ""))
selected_mask <- fc_mask | sig_mask

heatmap_sig <- heatmap_matrix[selected_mask, , drop = FALSE]
star_sig    <- star_matrix[selected_mask, , drop = FALSE]
max_fc      <- apply(abs(heatmap_sig), 1, max)

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
star_matrix_filtered    <- star_sig[selected_genes,    , drop = FALSE]

# Row order: group by functional category while preserving the existing order within each group
annotation_subset <- top_genes_df %>% arrange(Function.category)
sorted_rows <- annotation_subset$Gene
heatmap_matrix_filtered <- heatmap_matrix_filtered[sorted_rows, ]
star_matrix_filtered    <- star_matrix_filtered[sorted_rows, ]

row_annotation <- data.frame(Function = annotation_subset$Function.category)
rownames(row_annotation) <- sorted_rows

# Colors
function_list <- unique(row_annotation$Function)
high_contrast_colors <- createPalette(max(21, length(function_list)), seedcolors = c("#CCCCFF", "#339933"))
high_contrast_colors <- setNames(high_contrast_colors[seq_along(function_list)], function_list)
anno_colors <- list(Function = high_contrast_colors)

# Numbers: value plus significance stars
display_matrix <- matrix(
  paste0(round(heatmap_matrix_filtered, 2), " ", star_matrix_filtered),
  nrow = nrow(heatmap_matrix_filtered),
  ncol = ncol(heatmap_matrix_filtered),
  dimnames = list(rownames(heatmap_matrix_filtered), colnames(heatmap_matrix_filtered))
)

# ==== 10. Plot heatmap ====
pheatmap(heatmap_matrix_filtered,
         color = colorRampPalette(c("#339933", "white", "#FF3366"))(100),
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
         main = paste0("Functional Gene Expression log2FC (vs ", ref_phase, ")"),
         display_numbers = display_matrix,
         number_color = "black",
         fontsize_number = 28,
         fontsize = 28,
         fontsize_row = 38,
         fontsize_col = 38,
         cellheight = 35,
         cellwidth  = 85,
         annotation_legend_side = "right",
         annotation_legend_cols = 1,
         angle_col = 45
)

# ==== 11. Export results ====
write.csv(log2fc_data, "DC_log2fc_results.csv", row.names = FALSE)
write.csv(pvalue_data, "DC_pvalue_results.csv", row.names = FALSE)

# ==== 12. Friendly reminder: bins missing from the grouping file, if any ====
missing_bins <- setdiff(numeric_sample_cols, unique(group_data$Bins))
if (length(missing_bins) > 0) {
  message("The following sample columns do not have matching Phase information in group_data and will be ignored:\n",
          paste(missing_bins, collapse = ", "))
}