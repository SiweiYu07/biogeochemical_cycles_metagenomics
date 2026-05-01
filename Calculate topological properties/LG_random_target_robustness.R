# === 1. Load required packages ===
library(tidyverse)
library(igraph)
library(pracma)

# === 2. Read all files ===
files <- list.files(pattern = "^LG\\d{4}\\.csv$")

# === 3. Initialize the summary table ===
auc_summary <- data.frame(
  File = character(),
  Year = character(),
  AUC_random_bin = numeric(),
  AUC_random_KO = numeric(),
  AUC_targeted_bin = numeric(),
  AUC_targeted_KO = numeric(),
  stringsAsFactors = FALSE
)

# === 4. Define a function for one-time random removal ===
rand.remov.once <- function(netRaw, rm.percent, sp.ra, node = "bin", abundance.weighted = TRUE) {
  if (node == "bin") {
    id.rm <- sample(1:ncol(netRaw), round(ncol(netRaw) * rm.percent))
    net.Raw <- netRaw
    net.Raw[, id.rm] <- 0
    net.strength <- if (abundance.weighted) {
      net.Raw * matrix(sp.ra, nrow = nrow(netRaw), ncol = ncol(netRaw), byrow = FALSE)
    } else net.Raw
    remain.count <- sum(rowSums(net.strength) > 0)
    remain.percent <- remain.count / nrow(netRaw)
  } else {
    id.rm <- sample(1:nrow(netRaw), round(nrow(netRaw) * rm.percent))
    net.Raw <- netRaw
    net.Raw[id.rm, ] <- 0
    net.strength <- if (abundance.weighted) {
      net.Raw * matrix(sp.ra, nrow = nrow(netRaw), ncol = ncol(netRaw), byrow = TRUE)
    } else net.Raw
    remain.count <- sum(colSums(net.strength) > 0)
    remain.percent <- remain.count / ncol(netRaw)
  }
  return(remain.percent)
}

# === 5. Random removal simulation function ===
rmsimu <- function(netRaw, rm.p.list, sp.ra, node = "bin", abundance.weighted = TRUE, nperm = 100) {
  t(sapply(rm.p.list, function(x) {
    remains <- sapply(1:nperm, function(i) {
      rand.remov.once(netRaw, x, sp.ra, node, abundance.weighted)
    })
    remain.mean <- mean(remains)
    remain.sd <- sd(remains)
    remain.se <- remain.sd / sqrt(nperm)
    c(remain.mean = remain.mean, remain.sd = remain.sd, remain.se = remain.se)
  }))
}

# === 6. Targeted removal function ===
targeted_removal <- function(netRaw, top_nodes, node = "bin", abundance.weighted = TRUE, sp.ra) {
  if (node == "bin") {
    net.Raw <- netRaw
    net.Raw[, top_nodes] <- 0
    net.strength <- if (abundance.weighted) {
      net.Raw * matrix(sp.ra, nrow = nrow(netRaw), ncol = ncol(netRaw), byrow = FALSE)
    } else net.Raw
    remain.count <- sum(rowSums(net.strength) > 0)
    remain.percent <- remain.count / nrow(netRaw)
  } else {
    net.Raw <- netRaw
    net.Raw[top_nodes, ] <- 0
    net.strength <- if (abundance.weighted) {
      net.Raw * matrix(sp.ra, nrow = nrow(netRaw), ncol = ncol(netRaw), byrow = TRUE)
    } else net.Raw
    remain.count <- sum(colSums(net.strength) > 0)
    remain.percent <- remain.count / ncol(netRaw)
  }
  return(remain.percent)
}

# === 7. Main loop for processing each file ===
for (file in files) {
  year <- str_extract(file, "\\d{4}")
  cat("Processing:", file, "\n")
  
  # === Data import ===
  header <- read.csv(file, skip = 2, nrows = 1, header = FALSE)
  col_names <- as.character(header)
  df <- read.csv(file, skip = 3, header = FALSE)
  colnames(df) <- col_names
  
  tpm_data <- df[, 4:ncol(df)]
  rownames(tpm_data) <- df[[1]]
  tpm_data <- tpm_data[rowSums(tpm_data) > 0, ]
  tpm_data <- tpm_data[, colSums(tpm_data) > 0]
  if (nrow(tpm_data) < 2 | ncol(tpm_data) < 2) next
  
  binary_mat <- (tpm_data > 0) * 1
  network.raw <- as.matrix(binary_mat)
  
  KO_abund <- rowSums(tpm_data) / sum(rowSums(tpm_data))
  Bin_abund <- colSums(tpm_data) / sum(colSums(tpm_data))
  
  # === Build an igraph network for targeted removal ===
  edge_df <- binary_mat %>%
    as.data.frame() %>%
    rownames_to_column("Function") %>%
    pivot_longer(-Function, names_to = "Bin", values_to = "Presence") %>%
    filter(Presence > 0)
  
  g <- graph_from_data_frame(edge_df, directed = FALSE)
  V(g)$type <- V(g)$name %in% edge_df$Function
  
  # Get the top five high-degree nodes
  bin_nodes <- V(g)[V(g)$type == FALSE]
  ko_nodes <- V(g)[V(g)$type == TRUE]
  bin_deg <- sort(degree(g, v = bin_nodes), decreasing = TRUE)
  ko_deg <- sort(degree(g, v = ko_nodes), decreasing = TRUE)
  top_bin_idx <- match(names(bin_deg)[1:5], colnames(network.raw))
  top_ko_idx <- match(names(ko_deg)[1:5], rownames(network.raw))
  
  # === Run random removal simulations ===
  rm.p.list <- seq(0.05, 1, by = 0.05)
  simu_bin <- rmsimu(network.raw, rm.p.list, KO_abund, node = "bin", abundance.weighted = TRUE, nperm = 100)
  simu_KO <- rmsimu(network.raw, rm.p.list, Bin_abund, node = "KO", abundance.weighted = TRUE, nperm = 100)
  
  df_bin <- data.frame(Proportion.removed = rm.p.list, simu_bin, Type = "KO remaining after bin removal")
  df_KO <- data.frame(Proportion.removed = rm.p.list, simu_KO, Type = "Bin remaining after KO removal")
  df_all <- bind_rows(df_bin, df_KO)
  
  # === AUC calculation ===
  auc_bin <- trapz(rm.p.list, df_bin$remain.mean)
  auc_KO <- trapz(rm.p.list, df_KO$remain.mean)
  
  # === Targeted removal ===
  tar_bin <- targeted_removal(network.raw, top_bin_idx, node = "bin", abundance.weighted = TRUE, sp.ra = KO_abund)
  tar_KO <- targeted_removal(network.raw, top_ko_idx, node = "KO", abundance.weighted = TRUE, sp.ra = Bin_abund)
  
  # === Summarize AUC and targeted removal results ===
  auc_summary <- rbind(auc_summary, data.frame(
    File = file,
    Year = year,
    AUC_random_bin = round(auc_bin, 4),
    AUC_random_KO = round(auc_KO, 4),
    AUC_targeted_bin = round(tar_bin, 4),
    AUC_targeted_KO = round(tar_KO, 4)
  ))
  
  # === Write the robustness results for each individual file ===
  write.csv(df_all, paste0("robustness_random_removal_", year, ".csv"), row.names = FALSE)
}

# === Save the summary results across all years ===
write.csv(auc_summary, "robustness_auc_summary_all.csv", row.names = FALSE)