library(SiZer)
library(ggplot2)
library(segmented)

# Read data
seg <- read.csv("DC_segment.csv", check.names = FALSE)

# Make column names compatible
colnames(seg)[colnames(seg) == "Fe/Mn"] <- "Fe_Mn"

# Convert to numeric
seg$Fe_Mn <- as.numeric(seg$Fe_Mn)
seg$Nestedness <- as.numeric(seg$Nestedness)

# Remove missing values
seg <- seg[!is.na(seg$Fe_Mn) & !is.na(seg$Nestedness), ]

# First inspect the scatter plot and smoothed trend
ggplot(seg, aes(Fe_Mn, Nestedness)) +
  geom_point(size = 5, color = "#33CCCC") +
  geom_smooth(method = "loess", linewidth = 2, color = "blue", fill = "#FF9999") +
  xlab("Fe/Mn") +
  ylab("Nestedness") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 28),
    axis.text.y = element_text(size = 28),
    axis.title.x = element_text(size = 28),
    axis.title.y = element_text(size = 28),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# Linear model
fit_lm <- lm(Nestedness ~ Fe_Mn, data = seg)
summary(fit_lm)

# Davies test
dav <- davies.test(fit_lm, seg.Z = ~ Fe_Mn)
print(dav)

# Fit only one breakpoint
lm_seg1 <- segmented(fit_lm, seg.Z = ~ Fe_Mn, npsi = 1)
sum_seg <- summary(lm_seg1)
ci_bp <- confint(lm_seg1)

print(sum_seg)
print(ci_bp)

# Extract results
bp <- lm_seg1$psi[1, "Est."]
bp_se <- lm_seg1$psi[1, "St.Err"]

# Breakpoint ± standard error range
bp_low <- bp - bp_se
bp_high <- bp + bp_se

R2_lm <- summary(fit_lm)$r.squared
R2_seg <- summary(lm_seg1)$r.squared
p_dav <- dav$p.value

p_label <- ifelse(p_dav < 0.001,
                  "Davies test P < 0.001",
                  paste0("Davies test P = ", signif(p_dav, 3)))

# Plot
plot(lm_seg1,
     xlab = "Fe/Mn",
     ylab = "Nestedness",
     col = "blue",
     lwd = 7,
     cex.axis = 2,
     cex.lab = 2,
     axes = FALSE,
     frame.plot = FALSE)

axis(1, lwd = 2, lwd.ticks = 2, cex.axis = 2)
axis(2, lwd = 2, lwd.ticks = 2, cex.axis = 2)

box(lwd = 2, col = "black")

usr <- par("usr")

# Shaded area for breakpoint ± standard error
rect(xleft = bp_low,
     ybottom = usr[3],
     xright = bp_high,
     ytop = usr[4],
     col = rgb(1, 0, 0, alpha = 0.15),
     border = NA)

# Redraw the fitted line
plot(lm_seg1, add = TRUE, col = "#cc3333", lwd = 7)

# Add scatter points
points(Nestedness ~ Fe_Mn,
       data = seg,
       pch = 19,
       col = "#6eb5a1",
       cex = 2.7)

# Add vertical breakpoint line
abline(v = bp, lty = 2, lwd = 3, col = "red")

# Add text to the plot
text(x = usr[1] + 0.05 * (usr[2] - usr[1]),
     y = usr[4] - 0.08 * (usr[4] - usr[3]),
     labels = paste0("Breakpoint = ", round(bp, 2)),
     adj = c(0, 1),
     cex = 1.5)

text(x = usr[1] + 0.05 * (usr[2] - usr[1]),
     y = usr[4] - 0.16 * (usr[4] - usr[3]),
     labels = paste0("SE range = ", round(bp_low, 2), " to ", round(bp_high, 2)),
     adj = c(0, 1),
     cex = 1.5)

text(x = usr[1] + 0.05 * (usr[2] - usr[1]),
     y = usr[4] - 0.24 * (usr[4] - usr[3]),
     labels = paste0("Linear R2 = ", round(R2_lm, 3)),
     adj = c(0, 1),
     cex = 1.5)

text(x = usr[1] + 0.05 * (usr[2] - usr[1]),
     y = usr[4] - 0.32 * (usr[4] - usr[3]),
     labels = paste0("Segmented R2 = ", round(R2_seg, 3)),
     adj = c(0, 1),
     cex = 1.5)

text(x = usr[1] + 0.05 * (usr[2] - usr[1]),
     y = usr[4] - 0.40 * (usr[4] - usr[3]),
     labels = p_label,
     adj = c(0, 1),
     cex = 1.5)

dev.copy2pdf(file = "DC_Nestedness_FeMn_breakpoint.pdf", width = 7, height = 7)
dev.copy(png, filename = "DC_Nestedness_FeMn_breakpoint.png", width = 2000, height = 2000, res = 600)
dev.off()