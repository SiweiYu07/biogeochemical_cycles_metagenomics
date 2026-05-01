library(tidyverse)
library(mgcv)

#=========================
# 1. Read data
#=========================
df <- read.csv("LG_loess.csv", check.names = FALSE)

colnames(df)[1:3] <- c("Year_mid", "Nestedness", "Robustness")

df$Year_mid <- as.numeric(df$Year_mid)
df$Nestedness <- as.numeric(df$Nestedness)
df$Robustness <- as.numeric(df$Robustness)

df <- df %>% arrange(Year_mid)

#=========================
# 2. Convert to long format
#=========================
df_long <- df %>%
  pivot_longer(
    cols = c(Nestedness, Robustness),
    names_to = "Metric",
    values_to = "Value"
  )

df_long$Metric <- factor(
  df_long$Metric,
  levels = c("Nestedness", "Robustness")
)

#=========================
# 3. Fit GAM models separately
#=========================
metric_list <- levels(df_long$Metric)

fit_list <- lapply(metric_list, function(m) {
  sub <- df_long %>% filter(Metric == m)
  
  gam_fit <- gam(Value ~ s(Year_mid, k = 5), data = sub, method = "REML")
  sm <- summary(gam_fit)
  
  r2 <- sm$r.sq
  pval <- sm$s.table[1, "p-value"]
  
  newdat <- data.frame(
    Year_mid = seq(min(sub$Year_mid), max(sub$Year_mid), length.out = 200)
  )
  
  pred <- predict(gam_fit, newdata = newdat, se.fit = TRUE)
  
  pred_df <- data.frame(
    Year_mid = newdat$Year_mid,
    Fit = pred$fit,
    Lower = pred$fit - 1.96 * pred$se.fit,
    Upper = pred$fit + 1.96 * pred$se.fit,
    Metric = m
  )
  
  y_max <- max(sub$Value, na.rm = TRUE)
  y_min <- min(sub$Value, na.rm = TRUE)
  y_pos <- y_max - 0.03 * (y_max - y_min)
  
  p_text <- if (pval < 0.001) {
    "'<'*0.001"
  } else {
    sprintf("'='*%.3f", pval)
  }
  
  anno_df <- data.frame(
    Metric = m,
    x = min(sub$Year_mid) + 0.04 * (max(sub$Year_mid) - min(sub$Year_mid)),
    y = y_pos,
    label = paste0(
      "italic(R)^2~'='~", sprintf("%.3f", r2),
      "*','~~italic(P)~", p_text
    )
  )
  
  list(model = gam_fit, pred = pred_df, anno = anno_df)
})

pred_all <- bind_rows(lapply(fit_list, `[[`, "pred"))
anno_all <- bind_rows(lapply(fit_list, `[[`, "anno"))

# Optional: check model results
for (i in seq_along(metric_list)) {
  cat("\n============================\n")
  cat("Metric:", metric_list[i], "\n")
  print(summary(fit_list[[i]]$model))
}

#=========================
#=========================
# 4. Plot
#=========================
x_breaks <- seq(
  floor(min(df$Year_mid) / 20) * 20,
  ceiling(max(df$Year_mid) / 20) * 20,
  by = 20
)

p <- ggplot(df_long, aes(x = Year_mid, y = Value)) +
  geom_point(size = 6, color = "#897cff") +
  geom_line(linewidth = 2, color = "#897cff", alpha = 0.75) +
  geom_ribbon(
    data = pred_all,
    aes(x = Year_mid, ymin = Lower, ymax = Upper),
    inherit.aes = FALSE,
    fill = "#999999",
    alpha = 0.35
  ) +
  geom_line(
    data = pred_all,
    aes(x = Year_mid, y = Fit),
    inherit.aes = FALSE,
    color = "#CC3333",
    linewidth = 3
  ) +
  geom_text(
    data = anno_all,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    size = 8
  ) +
  facet_wrap(~Metric, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = x_breaks) +
  labs(
    x = "Midpoint year of sliding window",
    y = NULL
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.text = element_text(size = 30),
    axis.text.x = element_text(color = "black", size = 30),
    axis.text.y = element_text(color = "black", size = 30),
    axis.title.x = element_text(face = "bold", size = 30),
    axis.title.y = element_text(face = "bold", size = 30),
    axis.ticks = element_line(linewidth = 1),
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(20, 20, 20, 20)
  )

print(p)
ggsave("LG_Nestedness_Robustness_GAM_panels.pdf", p, width = 8, height = 7.5)
ggsave("LG_Nestedness_Robustness_GAM_panels.png", p, width = 8, height = 7.5, dpi = 600)