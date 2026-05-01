# =========================================================
# Causal-ML: Partial Residual Plot Generation for LG Lake
# =========================================================

# -----------------------------
# 0. Packages
# -----------------------------
pkg_needed <- c("dplyr", "purrr", "stringr", "ranger", "ggplot2", "readr", "tibble")
pkg_to_install <- pkg_needed[!pkg_needed %in% rownames(installed.packages())]
if(length(pkg_to_install) > 0) install.packages(pkg_to_install)

library(dplyr)
library(purrr)
library(stringr)
library(ranger)
library(ggplot2)
library(readr)
library(tibble)

set.seed(123)

# -----------------------------
# 1. Parameters
# -----------------------------
N_FOLDS   <- 5
N_TREES   <- 200 
MIN_NODE  <- 3
N_REPEATS <- 20
B_BOOT    <- 100
B_PERM    <- 100

# -----------------------------
# 2. Data Cleaning
# -----------------------------
read_lake <- function(file){
  df <- read.csv(file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  names(df) <- trimws(names(df))
  clean_names <- names(df) %>%
    gsub("\\.", "", .) %>% gsub("/", "", .) %>% gsub("\\(", "", .) %>%
    gsub("\\)", "", .) %>% gsub("\\+", "", .) %>% gsub("-", "", .) %>% gsub("\\s+", "", .)
  names(df) <- clean_names
  
  keep_cols <- intersect(c("Nestedness", "Robustness", "Precip", "Temp", "TOCTN", "TOC", "TN", "TP", 
                           "C", "N", "FeMn", "NiCo", "VCr", "CuMoZn", "MoAl"), names(df))
  df <- df[, keep_cols, drop = FALSE]
  for(j in seq_along(df)){ df[[j]] <- suppressWarnings(as.numeric(trimws(as.character(df[[j]])))) }
  df <- df[complete.cases(df[, intersect(c("Nestedness", "Robustness"), names(df)), drop = FALSE]), ]
  return(df)
}

# -----------------------------
# 3. Helpers
# -----------------------------
make_folds <- function(n, K = N_FOLDS, seed = 123){
  set.seed(seed)
  sample(rep(1:K, length.out = n))
}

cf_predict <- function(df, yvar, xvars, K = N_FOLDS, seed = 123){
  n <- nrow(df)
  folds <- make_folds(n, K = K, seed = seed)
  pred <- rep(NA_real_, n)
  for(k in seq_len(K)){
    train_id <- which(folds != k); test_id  <- which(folds == k)
    train_df <- df[train_id, c(yvar, xvars), drop = FALSE]
    test_df  <- df[test_id, xvars, drop = FALSE]
    form <- as.formula(paste(yvar, "~", paste(xvars, collapse = " + ")))
    fit <- ranger(formula = form, data = train_df, num.trees = N_TREES, 
                  min.node.size = MIN_NODE, seed = seed + k, num.threads = 2)
    pred[test_id] <- predict(fit, data = test_df)$predictions
  }
  return(pred)
}

# ---------------------------------------------------------
# 4. NEW: Partial Residual Plot Function (LG Specific)
# ---------------------------------------------------------
plot_partial_residual_lg <- function(df, outcome, treatment, controls, lake_name = "LG", seed = 123) {
  cat("Plotting LG Partial Residual:", treatment, "->", outcome, "\n")
  
  # 1. Obtain the residuals of the outcome after adjusting for environmental factors
  y_hat <- cf_predict(df, outcome, controls, seed = seed + 50)
  y_tilde <- df[[outcome]] - y_hat
  
  # 2. Obtain the Treatment residuals after controlling for environmental factors
  d_hat <- cf_predict(df, treatment, controls, seed = seed + 60)
  d_tilde <- df[[treatment]] - d_hat
  
  plot_df <- data.frame(d_resid = d_tilde, y_resid = y_tilde)
  
  # 3. plotting
  p <- ggplot(plot_df, aes(x = d_resid, y = y_resid)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
    geom_point(alpha = 2, color = "#2F74BB", size = 4) + # LG 使用绿色系点
    geom_smooth(method = "lm", color = "#D73027", fill = "#D73027", alpha = 0.15, size = 2) + 
    labs(
      x = paste0("Partial Residuals: ", treatment),
      y = paste0("Partial Residuals: ", outcome),
      title = paste("LG Lake:", treatment, "Effect")
    ) +
    theme_bw(base_size = 20) +
    theme(panel.grid = element_blank(), axis.title = element_text(face = "bold"))
  

  file_name <- paste0(lake_name, "_partial_", treatment, "_", outcome, ".pdf")
  ggsave(file_name, p, width = 5, height = 4)
}

# -----------------------------
# 5. Modified Analysis Function
# -----------------------------
analyze_lg_with_plots <- function(df, lake_name = "LG"){
  

  climate_vars  <- intersect(c("Precip", "Temp","TOCTN"), names(df))
  nutrient_vars <- intersect(c("TN", "TP"), names(df))
  redox_vars    <- intersect(c("FeMn", "MoAl", "CuMoZn"), names(df))
  outcomes      <- intersect(c("Nestedness", "Robustness"), names(df))
  controls      <- unique(c(climate_vars, nutrient_vars))
  

  cat("\n--- Starting LG Residual Analysis ---\n")
  for(rv in redox_vars){
    for(yy in outcomes){
      # plotting partial residual
      try(plot_partial_residual_lg(df, outcome = yy, treatment = rv, 
                                   controls = controls, lake_name = lake_name))
    }
  }
}

# -----------------------------
# 6. Execution
# -----------------------------
if(file.exists("LG.csv")){
  LG_data <- read_lake("LG.csv")
  analyze_lg_with_plots(LG_data, "LG")
  cat("\nSuccess: LG partial residual plots generated.\n")
} else {
  cat("Error: LG.csv not found in work directory.\n")
}