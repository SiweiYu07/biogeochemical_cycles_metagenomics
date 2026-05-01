# =========================================================
# Causal-ML style analysis for ONE lake only: DC.csv
#
# Step 1: Climate vs Nutrient -> each redox proxy
# Step 2: each redox proxy -> Nestedness / Robustness
#
# Output:
#   DC_driver_to_redox.csv
#   DC_redox_to_network.csv
#   DC_permutation_validation.csv
#   DC_summary_driver.csv
#   DC_summary_network.csv
#   DC_driver_to_redox.png / pdf
#   DC_redox_to_network.png / pdf
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
# 2. Read and clean DC data
# -----------------------------
read_lake <- function(file){
  
  df <- read.csv(
    file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )
  
  names(df) <- trimws(names(df))
  
  clean_names <- names(df)
  clean_names <- gsub("\\.", "", clean_names)
  clean_names <- gsub("/", "", clean_names)
  clean_names <- gsub("\\(", "", clean_names)
  clean_names <- gsub("\\)", "", clean_names)
  clean_names <- gsub("\\+", "", clean_names)
  clean_names <- gsub("-", "", clean_names)
  clean_names <- gsub("\\s+", "", clean_names)
  
  names(df) <- clean_names
  
  keep_cols <- intersect(
    c("Nestedness", "Robustness",
      "Precip", "Temp",
      "TOCTN", "TOC", "TN", "TP", "C", "N",
      "FeMn", "NiCo", "VCr", "CuMoZn", "MoAl"),
    names(df)
  )
  
  df <- df[, keep_cols, drop = FALSE]
  
  for(j in seq_along(df)){
    df[[j]] <- suppressWarnings(as.numeric(trimws(as.character(df[[j]]))))
  }
  

  need_core <- intersect(c("Nestedness", "Robustness"), names(df))
  if(length(need_core) > 0){
    keep <- complete.cases(df[, need_core, drop = FALSE])
    df <- df[keep, , drop = FALSE]
  }
  
  rownames(df) <- NULL
  
  cat("\n========================\n")
  cat("Read file:", file, "\n")
  cat("Dimensions:", nrow(df), "x", ncol(df), "\n")
  cat("Columns:\n")
  print(names(df))
  cat("========================\n")
  
  return(df)
}

# -----------------------------
# 3. Helpers
# -----------------------------
make_folds <- function(n, K = N_FOLDS, seed = 123){
  set.seed(seed)
  sample(rep(1:K, length.out = n))
}

oof_r2 <- function(y, yhat){
  denom <- sum((y - mean(y))^2)
  if(denom == 0) return(NA_real_)
  1 - sum((y - yhat)^2) / denom
}

safe_intersect <- function(x, y){
  intersect(x, y)
}

# -----------------------------
# 4. Cross-fitted prediction
# -----------------------------
cf_predict <- function(df, yvar, xvars,
                       K = N_FOLDS, seed = 123,
                       num.trees = N_TREES,
                       min.node.size = MIN_NODE){
  
  stopifnot(yvar %in% names(df))
  stopifnot(all(xvars %in% names(df)))
  
  n <- nrow(df)
  folds <- make_folds(n, K = K, seed = seed)
  pred <- rep(NA_real_, n)
  
  for(k in seq_len(K)){
    train_id <- which(folds != k)
    test_id  <- which(folds == k)
    
    train_df <- df[train_id, c(yvar, xvars), drop = FALSE]
    test_df  <- df[test_id, xvars, drop = FALSE]
    
    form <- as.formula(paste(yvar, "~", paste(xvars, collapse = " + ")))
    
    fit <- ranger(
      formula = form,
      data = train_df,
      num.trees = num.trees,
      min.node.size = min.node.size,
      seed = seed + k,
      num.threads = 2
    )
    
    pred[test_id] <- predict(fit, data = test_df)$predictions
  }
  
  pred
}

# -----------------------------
# 5. Step 1: driver -> redox
# repeated cross-fitted delta-R2
# -----------------------------
cf_group_effect_once <- function(df, yvar, base_vars, add_vars, seed = 123){
  
  base_vars <- safe_intersect(base_vars, names(df))
  add_vars  <- safe_intersect(add_vars, names(df))
  
  if(length(base_vars) == 0 || length(add_vars) == 0){
    return(tibble(
      redox = yvar,
      r2_base = NA_real_,
      r2_full = NA_real_,
      delta_r2 = NA_real_
    ))
  }
  
  pred_base <- cf_predict(df, yvar, base_vars, seed = seed)
  pred_full <- cf_predict(df, yvar, unique(c(base_vars, add_vars)), seed = seed + 100)
  
  tibble(
    redox = yvar,
    r2_base = oof_r2(df[[yvar]], pred_base),
    r2_full = oof_r2(df[[yvar]], pred_full),
    delta_r2 = oof_r2(df[[yvar]], pred_full) - oof_r2(df[[yvar]], pred_base)
  )
}

cf_group_effect_repeat <- function(df, yvar, base_vars, add_vars,
                                   n_repeats = N_REPEATS, seed = 123){
  
  outs <- map_dfr(seq_len(n_repeats), function(i){
    cf_group_effect_once(df, yvar, base_vars, add_vars, seed = seed + i * 1000) %>%
      mutate(rep = i)
  })
  
  tibble(
    redox = yvar,
    r2_base = mean(outs$r2_base, na.rm = TRUE),
    r2_full = mean(outs$r2_full, na.rm = TRUE),
    delta_r2 = mean(outs$delta_r2, na.rm = TRUE),
    delta_r2_sd = sd(outs$delta_r2, na.rm = TRUE)
  )
}

boot_group_effect <- function(df, yvar, base_vars, add_vars,
                              B = B_BOOT, n_repeats = 10, seed = 123){
  
  set.seed(seed)
  vals <- numeric(B)
  
  for(b in seq_len(B)){
    if(b %% 50 == 0) cat("bootstrap driver:", yvar, b, "/", B, "\n")
    
    id <- sample(seq_len(nrow(df)), replace = TRUE)
    dsub <- df[id, , drop = FALSE]
    
    out <- cf_group_effect_repeat(
      dsub, yvar = yvar,
      base_vars = base_vars,
      add_vars = add_vars,
      n_repeats = n_repeats,
      seed = seed + b
    )
    
    vals[b] <- out$delta_r2
  }
  
  tibble(
    redox = yvar,
    ci_low = quantile(vals, 0.025, na.rm = TRUE),
    ci_high = quantile(vals, 0.975, na.rm = TRUE)
  )
}

# -----------------------------
# 6. Step 2: redox -> network
# continuous-treatment partial effect
# -----------------------------
dml_partial_effect_once <- function(df, outcome, treatment, controls, seed = 123){
  
  controls <- safe_intersect(controls, names(df))
  
  if(length(controls) == 0){
    fit <- lm(as.formula(paste(outcome, "~", treatment)), data = df)
    ss <- summary(fit)
    return(tibble(
      outcome = outcome,
      treatment = treatment,
      beta = coef(ss)[treatment, "Estimate"],
      se = coef(ss)[treatment, "Std. Error"],
      t = coef(ss)[treatment, "t value"],
      p = coef(ss)[treatment, "Pr(>|t|)"],
      r2 = ss$r.squared
    ))
  }
  
  y_hat <- cf_predict(df, outcome, controls, seed = seed + 1)
  d_hat <- cf_predict(df, treatment, controls, seed = seed + 2)
  
  y_tilde <- df[[outcome]] - y_hat
  d_tilde <- df[[treatment]] - d_hat
  
  fit <- lm(y_tilde ~ d_tilde)
  ss <- summary(fit)
  
  tibble(
    outcome = outcome,
    treatment = treatment,
    beta = coef(ss)["d_tilde", "Estimate"],
    se = coef(ss)["d_tilde", "Std. Error"],
    t = coef(ss)["d_tilde", "t value"],
    p = coef(ss)["d_tilde", "Pr(>|t|)"],
    r2 = ss$r.squared
  )
}

dml_partial_effect_repeat <- function(df, outcome, treatment, controls,
                                      n_repeats = N_REPEATS, seed = 123){
  
  outs <- map_dfr(seq_len(n_repeats), function(i){
    dml_partial_effect_once(
      df, outcome = outcome, treatment = treatment,
      controls = controls, seed = seed + i * 1000
    ) %>% mutate(rep = i)
  })
  
  tibble(
    outcome = outcome,
    treatment = treatment,
    beta = mean(outs$beta, na.rm = TRUE),
    beta_sd = sd(outs$beta, na.rm = TRUE),
    se = mean(outs$se, na.rm = TRUE),
    t = mean(outs$t, na.rm = TRUE),
    p = mean(outs$p, na.rm = TRUE),
    r2 = mean(outs$r2, na.rm = TRUE)
  )
}

boot_dml_effect <- function(df, outcome, treatment, controls,
                            B = B_BOOT, n_repeats = 10, seed = 123){
  
  set.seed(seed)
  vals <- numeric(B)
  
  for(b in seq_len(B)){
    if(b %% 50 == 0) cat("bootstrap DML:", treatment, "->", outcome, b, "/", B, "\n")
    
    id <- sample(seq_len(nrow(df)), replace = TRUE)
    dsub <- df[id, , drop = FALSE]
    
    out <- dml_partial_effect_repeat(
      dsub, outcome = outcome, treatment = treatment,
      controls = controls, n_repeats = n_repeats,
      seed = seed + b
    )
    vals[b] <- out$beta
  }
  
  tibble(
    outcome = outcome,
    treatment = treatment,
    ci_low = quantile(vals, 0.025, na.rm = TRUE),
    ci_high = quantile(vals, 0.975, na.rm = TRUE)
  )
}

# -----------------------------
# 7. Permutation validation
# -----------------------------
perm_test_dml <- function(df, outcome, treatment, controls,
                          B = B_PERM, n_repeats = 10, seed = 123){
  
  obs <- dml_partial_effect_repeat(
    df, outcome = outcome, treatment = treatment,
    controls = controls, n_repeats = n_repeats, seed = seed
  )$beta
  
  set.seed(seed)
  null_vals <- numeric(B)
  
  for(b in seq_len(B)){
    if(b %% 50 == 0) cat("permutation:", treatment, "->", outcome, b, "/", B, "\n")
    
    dperm <- df
    dperm[[treatment]] <- sample(dperm[[treatment]])
    
    out <- dml_partial_effect_repeat(
      dperm, outcome = outcome, treatment = treatment,
      controls = controls, n_repeats = 5,
      seed = seed + b
    )
    null_vals[b] <- out$beta
  }
  
  p_perm <- mean(abs(null_vals) >= abs(obs), na.rm = TRUE)
  
  tibble(
    outcome = outcome,
    treatment = treatment,
    beta_obs = obs,
    perm_p = p_perm,
    null_mean = mean(null_vals, na.rm = TRUE),
    null_sd = sd(null_vals, na.rm = TRUE)
  )
}

# -----------------------------
# 8. Analyze DC only
# -----------------------------
analyze_dc <- function(df, lake_name = "DC"){
  
  climate_vars  <- safe_intersect(c( "Precip","Temp"), names(df))
  nutrient_vars <- safe_intersect(c("TN", "TP"), names(df))
  redox_vars    <- safe_intersect(c("FeMn", "CuMoZn", "MoAl"), names(df))
  outcomes      <- safe_intersect(c("Nestedness", "Robustness"), names(df))
  
  cat("\n===================================\n")
  cat("Analyzing lake:", lake_name, "\n")
  cat("Climate vars :", paste(climate_vars, collapse = ", "), "\n")
  cat("Nutrient vars:", paste(nutrient_vars, collapse = ", "), "\n")
  cat("Redox vars   :", paste(redox_vars, collapse = ", "), "\n")
  cat("Outcomes     :", paste(outcomes, collapse = ", "), "\n")
  cat("===================================\n")
  
  # ---- Step 1: drivers -> redox
  driver_res <- list()
  
  for(rv in redox_vars){
    
    cat("\n[Step 1] Climate ->", rv, "\n")
    est_c <- cf_group_effect_repeat(
      df, yvar = rv,
      base_vars = nutrient_vars,
      add_vars = climate_vars
    )
    ci_c <- boot_group_effect(
      df, yvar = rv,
      base_vars = nutrient_vars,
      add_vars = climate_vars
    )
    
    out_c <- left_join(est_c, ci_c, by = "redox") %>%
      mutate(driver = "Climate", lake = lake_name)
    
    cat("\n[Step 1] Nutrient ->", rv, "\n")
    est_n <- cf_group_effect_repeat(
      df, yvar = rv,
      base_vars = climate_vars,
      add_vars = nutrient_vars
    )
    ci_n <- boot_group_effect(
      df, yvar = rv,
      base_vars = climate_vars,
      add_vars = nutrient_vars
    )
    
    out_n <- left_join(est_n, ci_n, by = "redox") %>%
      mutate(driver = "Nutrient", lake = lake_name)
    
    driver_res[[rv]] <- bind_rows(out_c, out_n)
  }
  
  driver_results <- bind_rows(driver_res) %>%
    relocate(lake, driver, redox)
  
  # ---- Step 2: redox -> network
  controls <- unique(c(climate_vars, nutrient_vars))
  network_res <- list()
  perm_res <- list()
  
  for(rv in redox_vars){
    for(yy in outcomes){
      
      cat("\n[Step 2]", rv, "->", yy, "\n")
      
      est <- dml_partial_effect_repeat(
        df, outcome = yy, treatment = rv, controls = controls
      )
      ci  <- boot_dml_effect(
        df, outcome = yy, treatment = rv, controls = controls
      )
      pm  <- perm_test_dml(
        df, outcome = yy, treatment = rv, controls = controls
      )
      
      network_res[[paste(rv, yy, sep = "_")]] <-
        left_join(est, ci, by = c("outcome", "treatment")) %>%
        mutate(lake = lake_name)
      
      perm_res[[paste(rv, yy, sep = "_")]] <-
        pm %>% mutate(lake = lake_name)
    }
  }
  
  network_results <- bind_rows(network_res) %>%
    relocate(lake, outcome, treatment)
  
  perm_results <- bind_rows(perm_res) %>%
    relocate(lake, outcome, treatment)
  
  # ---- Summary
  summary_driver <- driver_results %>%
    mutate(
      sig = ifelse(ci_low > 0 | ci_high < 0, "yes", "no"),
      stable = ifelse(!is.na(delta_r2_sd) & abs(delta_r2) > delta_r2_sd, "yes", "no")
    ) %>%
    arrange(redox, desc(delta_r2))
  
  summary_network <- network_results %>%
    left_join(
      perm_results %>% select(lake, outcome, treatment, perm_p),
      by = c("lake", "outcome", "treatment")
    ) %>%
    mutate(
      sig = ifelse(ci_low > 0 | ci_high < 0, "yes", "no"),
      stable = ifelse(!is.na(beta_sd) & abs(beta) > beta_sd, "yes", "no"),
      perm_sig = ifelse(perm_p < 0.05, "yes", "no")
    ) %>%
    arrange(outcome, desc(abs(beta)))
  
  # ---- Save results
  write_csv(driver_results,  "DC_driver_to_redox.csv")
  write_csv(network_results, "DC_redox_to_network.csv")
  write_csv(perm_results,    "DC_permutation_validation.csv")
  write_csv(summary_driver,  "DC_summary_driver.csv")
  write_csv(summary_network, "DC_summary_network.csv")
  
  # ---- Plot 1
  p1 <- ggplot(driver_results,
               aes(x = redox, y = delta_r2, fill = driver)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                  position = position_dodge(width = 0.7), width = 0.2) +
    labs(
      x = "Redox proxy",
      y = "Incremental cross-fitted R²",
      title = "DC: climate/nutrient effects on redox"
    ) +
    theme_classic(base_size = 12)
  
  ggsave("DC_driver_to_redox.png", p1, width = 7, height = 5, dpi = 300)
  ggsave("DC_driver_to_redox.pdf", p1, width = 7, height = 5)
  
  # ---- Plot 2
  p2 <- ggplot(network_results,
               aes(x = treatment, y = beta, fill = outcome)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                  position = position_dodge(width = 0.7), width = 0.2) +
    labs(
      x = "Redox proxy",
      y = "Cross-fitted partial effect",
      title = "DC: redox effects on network metrics"
    ) +
    theme_classic(base_size = 12)
  
  ggsave("DC_redox_to_network.png", p2, width = 7, height = 5, dpi = 300)
  ggsave("DC_redox_to_network.pdf", p2, width = 7, height = 5)
  
  return(list(
    driver_results = driver_results,
    network_results = network_results,
    perm_results = perm_results,
    summary_driver = summary_driver,
    summary_network = summary_network
  ))
}

# -----------------------------
# 9. Run DC only
# -----------------------------
DC <- read_lake("DC.csv")
res_DC <- analyze_dc(DC, "DC")

# -----------------------------
# 10. Print key tables
# -----------------------------
print(res_DC$summary_driver)
print(res_DC$summary_network)