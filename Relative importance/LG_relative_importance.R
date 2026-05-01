library(relaimpo)
library(car)
library(ggplot2)
library(dplyr)
library(readr)
library(gridExtra)

#=============================
# 2. Read data
#=============================
df <- read.csv("LG.csv", check.names = FALSE)

str(df)

#=============================
# 3. Select variables
#=============================
df2 <- df %>%
  dplyr::select(
    Nestedness,
    Robustness,
    `Precip.`,
    `Temp.`,
    `TN`,
    `C`,
    `Fe/Mn`,
    `Mo/Al`,
    `(Cu+Mo)/Zn`,
    `Ni/Co`
  )

# Convert all columns to numeric
df2[] <- lapply(df2, as.numeric)

# Remove missing values
df2 <- na.omit(df2)

#=============================
# 4. Standardize predictors
#=============================
predictor_names <- c("Precip.", "Temp.", "TN", "C", "Fe/Mn", "Mo/Al", "(Cu+Mo)/Zn", "Ni/Co")
df2[predictor_names] <- scale(df2[predictor_names])

#=============================
# 5. Check collinearity
#=============================
fit_vif_nested <- lm(
  Nestedness ~ `Precip.` + `Temp.` + `TN` + `C` + `Fe/Mn` + `Mo/Al` + `(Cu+Mo)/Zn` + `Ni/Co`,
  data = df2
)

fit_vif_robust <- lm(
  Robustness ~ `Precip.` + `Temp.` + `TN` + `C` + `Fe/Mn` + `Mo/Al` + `(Cu+Mo)/Zn` + `Ni/Co`,
  data = df2
)

cat("\n===== VIF for Nestedness model =====\n")
print(car::vif(fit_vif_nested))

cat("\n===== VIF for Robustness model =====\n")
print(car::vif(fit_vif_robust))

#=============================
# 6. Model fitting
#=============================
fit_nested <- lm(
  Nestedness ~ `Precip.` + `Temp.` + `TN` + `C` + `Fe/Mn` + `Mo/Al` + `(Cu+Mo)/Zn` + `Ni/Co`,
  data = df2
)

fit_robust <- lm(
  Robustness ~ `Precip.` + `Temp.` + `TN` + `C` + `Fe/Mn` + `Mo/Al` + `(Cu+Mo)/Zn` + `Ni/Co`,
  data = df2
)

cat("\n===== Summary: Nestedness model =====\n")
print(summary(fit_nested))

cat("\n===== Summary: Robustness model =====\n")
print(summary(fit_robust))

#=============================
# 7. Extract model statistics
#=============================
get_model_stats <- function(fit) {
  sm <- summary(fit)
  
  r2 <- sm$r.squared
  adj_r2 <- sm$adj.r.squared
  f_value <- sm$fstatistic[1]
  df1 <- sm$fstatistic[2]
  df2 <- sm$fstatistic[3]
  p_value <- pf(f_value, df1, df2, lower.tail = FALSE)
  
  stats_label <- paste0(
    "R² = ", sprintf("%.3f", r2), "\n",
    "Adj.R² = ", sprintf("%.3f", adj_r2), "\n",
    "F = ", sprintf("%.2f", f_value), "\n",
    "P = ", ifelse(p_value < 0.001, "< 0.001", sprintf("%.3f", p_value))
  )
  
  return(list(
    r2 = r2,
    adj_r2 = adj_r2,
    f_value = f_value,
    p_value = p_value,
    label = stats_label
  ))
}

stats_nested <- get_model_stats(fit_nested)
stats_robust <- get_model_stats(fit_robust)

cat("\n===== Model stats: Nestedness =====\n")
print(stats_nested)

cat("\n===== Model stats: Robustness =====\n")
print(stats_robust)

#=============================
# 8. Relative importance analysis
#=============================
rel_nested <- calc.relimp(fit_nested, type = "lmg", rela = TRUE)
rel_robust <- calc.relimp(fit_robust, type = "lmg", rela = TRUE)

cat("\n===== Relative importance: Nestedness =====\n")
print(rel_nested)

cat("\n===== Relative importance: Robustness =====\n")
print(rel_robust)

#=============================
# 9. Bootstrap confidence intervals
#=============================
set.seed(123)

boot_nested <- boot.relimp(
  fit_nested,
  b = 1000,
  type = "lmg",
  rank = TRUE,
  diff = TRUE,
  rela = TRUE
)

boot_robust <- boot.relimp(
  fit_robust,
  b = 1000,
  type = "lmg",
  rank = TRUE,
  diff = TRUE,
  rela = TRUE
)

cat("\n===== Bootstrap evaluation: Nestedness =====\n")
print(booteval.relimp(boot_nested))

cat("\n===== Bootstrap evaluation: Robustness =====\n")
print(booteval.relimp(boot_robust))

#=============================
# 10. Organize result tables
#=============================
nested_df <- data.frame(
  Variable = names(rel_nested$lmg),
  Importance = as.numeric(rel_nested$lmg) * 100
)

robust_df <- data.frame(
  Variable = names(rel_robust$lmg),
  Importance = as.numeric(rel_robust$lmg) * 100
)

nested_df <- nested_df %>% arrange(desc(Importance))
robust_df <- robust_df %>% arrange(desc(Importance))

write.csv(nested_df, "LG_Nestedness_relative_importance.csv", row.names = FALSE)
write.csv(robust_df, "LG_Robustness_relative_importance.csv", row.names = FALSE)

#=============================
# 11. Plotting function with model statistics annotation
#=============================
plot_relimp <- function(dat, title_text, fill_color, stats_obj) {
  
  xmax <- max(dat$Importance) * 1.28   # Leave space on the right for text
  ymax <- nrow(dat)
  
  ggplot(dat, aes(x = reorder(Variable, Importance), y = Importance)) +
    geom_col(fill = fill_color, width = 0.6) +
    geom_text(
      aes(label = sprintf("%.1f%%", Importance)),
      hjust = 0.0,
      size = 4
    ) +
    coord_flip() +
    annotate(
      "text",
      x = 3.0,                       # Vertical position after coordinate flipping
      y = xmax * 0.4,                # Horizontal position
      label = stats_obj$label,
      hjust = 0,
      vjust = 1,
      size = 4.5,
      fontface = "plain"
    ) +
    scale_y_continuous(
      limits = c(0, xmax),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = NULL,
      y = "Relative importance (%)",
      title = title_text
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 17, hjust = 0.5),
      axis.text.x = element_text(size = 17, color = "black"),
      axis.text.y = element_text(size = 17, color = "black"),
      axis.title = element_text(size = 17),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
}

p1 <- plot_relimp(nested_df, "Nestedness", "#897Cff", stats_nested)
p2 <- plot_relimp(robust_df, "Robustness", "#897Cff", stats_robust)

#=============================
# 12. Save individual plots
#=============================
ggsave("LG_Nestedness_relative_importance.jpg", p1, width = 6, height = 3.5, dpi = 600)
ggsave("LG_Nestedness_relative_importance.pdf", p1, width = 6, height = 3.5)

ggsave("LG_Robustness_relative_importance.jpg", p2, width = 6, height = 3.5, dpi = 600)
ggsave("LG_Robustness_relative_importance.pdf", p2, width = 6, height = 3.5)

#=============================
# 13. Save combined plot
#=============================
p_combined <- gridExtra::grid.arrange(p1, p2, ncol = 2)

ggsave("LG_relative_importance_combined.jpg", p_combined, width = 6, height = 3.6, dpi = 600)
ggsave("LG_relative_importance_combined.pdf", p_combined, width = 6, height = 3.6)

#=============================
# 14. Output message
#=============================
cat("\nAnalysis completed. Generated files:\n")
cat("LG_Nestedness_relative_importance.csv\n")
cat("LG_Robustness_relative_importance.csv\n")
cat("LG_Nestedness_relative_importance.jpg / pdf\n")
cat("LG_Robustness_relative_importance.jpg / pdf\n")
cat("LG_relative_importance_combined.jpg / pdf\n")