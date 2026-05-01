# === 1. Load required packages ===
library(tidyverse)
library(igraph)
library(readr)
library(bipartite)  # For NODF calculation
library(vegan)      # For null model generation

# === 2. Read the list of data files matching the pattern, e.g., LG1902.csv ===
files <- list.files(pattern = "^LG\\d{4}\\.csv$")

# === 3. Initialize the output table ===
summary_result <- data.frame(
  File = character(),
  Period = character(),
  Nestedness = numeric(),
  Null_Nestedness = numeric(),
  Relative_Nestedness = numeric(),
  stringsAsFactors = FALSE
)

# === 4. Main loop ===
for (file in files) {
  # Extract the year
  period <- stringr::str_extract(file, "\\d{4}")
  
  # Read column names from the third row, then read data from the fourth row onward
  header <- read.csv(file, skip = 2, nrows = 1, header = FALSE)
  col_names <- as.character(header)
  df <- read.csv(file, skip = 3, header = FALSE)
  colnames(df) <- col_names
  
  # Extract the TPM matrix from the fourth column onward, using KOs as row names
  tpm_data <- df[, 4:ncol(df)]
  rownames(tpm_data) <- df[[1]]
  
  # Remove KOs with a total TPM of zero
  tpm_data <- tpm_data[rowSums(tpm_data) > 0, ]
  if (nrow(tpm_data) == 0) next
  
  # Build the edge list
  edge_df <- tpm_data %>%
    as.data.frame() %>%
    rownames_to_column("Function") %>%
    pivot_longer(-Function, names_to = "Bin", values_to = "TPM") %>%
    filter(TPM > 0)
  
  # === 5. Build the presence–absence bipartite matrix ===
  pa_matrix <- edge_df %>%
    distinct(Function, Bin) %>%           # Remove duplicated Function–Bin pairs
    mutate(value = 1) %>%
    pivot_wider(names_from = Bin, values_from = value, values_fill = 0) %>%
    column_to_rownames("Function") %>%
    as.matrix()
  
  if (nrow(pa_matrix) < 2 || ncol(pa_matrix) < 2) next  # Avoid errors
  
  # === 6. Calculate NODF nestedness ===
  nodf_obs <- nested(pa_matrix, method = "NODF") / 100  # Normalize to the 0–1 range
  
  # === 7. Generate null models 100 times and calculate mean nestedness Nr ===
  set.seed(123)
  null_vals <- replicate(100, {
    null_matrix <- permatfull(pa_matrix, mtype = "prab", fixedmar = "both", times = 1)$perm[[1]]
    nested(null_matrix, method = "NODF") / 100
  })
  nodf_null_mean <- mean(null_vals)
  
  # === 8. Calculate relative nestedness ===
  RN <- if (nodf_null_mean == 0) NA else (nodf_obs - nodf_null_mean) / nodf_null_mean
  
  # === 9. Save results ===
  summary_result <- rbind(summary_result, data.frame(
    File = file,
    Period = period,
    Nestedness = round(nodf_obs, 4),
    Null_Nestedness = round(nodf_null_mean, 4),
    Relative_Nestedness = round(RN, 4)
  ))
}

# === 10. Write the results to a CSV file ===
write.csv(summary_result, "KO_bin_nestedness_summary.csv", row.names = FALSE)

# === 11. Visualize the temporal trend of nestedness ===
ggplot(summary_result, aes(x = as.numeric(Period), y = Relative_Nestedness)) +
  geom_line(group = 1, color = "purple") +
  geom_point(size = 3, color = "darkblue") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Year",
    y = "Relative Nestedness (RN)",
    title = "Temporal Dynamics of Functional Nestedness (KO-bin Network)"
  )