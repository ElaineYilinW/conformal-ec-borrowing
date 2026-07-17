#!/usr/bin/env Rscript

# Reproducible figures for the four longitudinal-simulation slides.

source(file.path("sim", "engine_longitudinal.R"))

out_dir <- file.path("CC_ppt", "longitudinal_assets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

navy <- "#17365D"
blue <- "#2F75B5"
light_blue <- "#9DC3E6"
gray <- "#7F8C8D"
red <- "#C74440"
green <- "#2E8B57"

# Seed 2 is a representative draw: detection 94.1%, wrong exclusion 9.1%,
# close to the 200-replication averages 93.4% and 8.3%.
set.seed(2)
dat <- long_generate_data()
screen <- long_screen_cv(dat, score_configs = long_score_configs(FALSE))

fit <- long_fit_outcome(dat$Yc, dat$Xc, longitudinal = TRUE)
rr <- dat$Yc - fit$predict_matrix(dat$Xc)
re <- dat$Ye - fit$predict_matrix(dat$Xe)
clean <- dat$B == 0
bad <- dat$B == 1

mean_se <- function(x) {
  cbind(mean = colMeans(x), se = apply(x, 2, sd) / sqrt(nrow(x)))
}
m_r <- mean_se(rr)
m_c <- mean_se(re[clean, , drop = FALSE])
m_b <- mean_se(re[bad, , drop = FALSE])

png(file.path(out_dir, "one_draw_two_sided.png"), width = 1600, height = 1450,
    res = 180, bg = "white")
par(mfrow = c(2, 1), mar = c(4.6, 5.2, 3.7, 1.4), oma = c(0, 0, 1.0, 0),
    family = "HersheySans", las = 1)

ylim <- range(c(m_r[, "mean"] + 1.96 * c(-m_r[, "se"], m_r[, "se"]),
                m_c[, "mean"] + 1.96 * c(-m_c[, "se"], m_c[, "se"]),
                m_b[, "mean"] + 1.96 * c(-m_b[, "se"], m_b[, "se"])))
ylim <- c(min(-8.5, ylim[1]), max(1.5, ylim[2]))
plot(LONG_TIMES, m_r[, "mean"], type = "n", xlim = c(1, 12), ylim = ylim,
     xlab = "Visit", ylab = "Mean residual after RCT-control GLS",
     main = "All 12 visits expose the incompatible trajectory", axes = FALSE,
     cex.main = 1.12, col.main = navy)
axis(1, at = LONG_TIMES)
axis(2)
box(col = "#D9E2F3")
abline(h = 0, lty = 2, col = "#B7B7B7")
band <- function(m, col) {
  polygon(c(LONG_TIMES, rev(LONG_TIMES)),
          c(m[, "mean"] - 1.96 * m[, "se"], rev(m[, "mean"] + 1.96 * m[, "se"])),
          col = adjustcolor(col, alpha.f = 0.16), border = NA)
  lines(LONG_TIMES, m[, "mean"], col = col, lwd = 3)
  points(LONG_TIMES, m[, "mean"], col = col, pch = 16, cex = 0.65)
}
band(m_r, gray)
band(m_c, blue)
band(m_b, red)
legend(x = 7.2, y = -1.1, c("RCT controls", "Compatible EC", "Incompatible EC"),
       col = c(gray, blue, red), lwd = 3, pch = 16, bty = "n", cex = 0.93)

p_clean <- screen$p$avg[clean]
p_bad <- screen$p$avg[bad]
br <- seq(0, 1, by = 0.1)
h_clean <- hist(p_clean, breaks = br, plot = FALSE)
h_bad <- hist(p_bad, breaks = br, plot = FALSE)
ymax <- 1.12 * max(h_clean$counts, h_bad$counts)
plot(NA, xlim = c(0, 1), ylim = c(0, ymax), xlab = "Weighted conformal p-value",
     ylab = "External controls", main = "Two-sided p-values separate bad ECs without one-tail trimming",
     axes = FALSE, cex.main = 1.12, col.main = navy)
axis(1, at = seq(0, 1, 0.1))
axis(2)
box(col = "#D9E2F3")
rect(br[-length(br)] + 0.005, 0, br[-1] - 0.005, h_clean$counts,
     col = adjustcolor(blue, alpha.f = 0.58), border = blue)
rect(br[-length(br)] + 0.025, 0, br[-1] - 0.025, h_bad$counts,
     col = adjustcolor(red, alpha.f = 0.62), border = red)
abline(v = dat$config$gamma, lwd = 3, lty = 2, col = navy)
text(dat$config$gamma + 0.015, ymax * 0.93, "gamma = 0.10",
     pos = 4, col = navy, font = 2)
text(0.045, ymax * 0.80, "drop", col = red, font = 2)
text(0.57, ymax * 0.80, "retain", col = green, font = 2)
legend("topright", c(sprintf("Compatible EC (n=%d)", sum(clean)),
                      sprintf("Incompatible EC (n=%d)", sum(bad))),
       fill = c(adjustcolor(blue, alpha.f = 0.58), adjustcolor(red, alpha.f = 0.62)),
       border = c(blue, red), bty = "n", cex = 0.90)
dev.off()

# Score comparison uses the actual 200-replication summaries.
scr <- read.csv(file.path("sim", "results_longitudinal_main", "longitudinal_main_screening.csv"))
want <- c("avg", "final", "c1", "max", "quad")
labels <- c(avg = "Average / internal CS", final = "Visit 12 only",
            c1 = "Visit 12 - visit 1", max = "Maximum / diagonal",
            quad = "Quadratic / diagonal")
scr <- scr[match(want, scr$score), ]
vals <- rbind(scr$detection, scr$wrong_detection)

png(file.path(out_dir, "score_comparison.png"), width = 1900, height = 570,
    res = 180, bg = "white")
par(mar = c(3.8, 12.0, 2.8, 1.2), family = "HersheySans", las = 1)
bp <- barplot(vals, beside = TRUE, horiz = TRUE, names.arg = unname(labels[want]),
              xlim = c(0, 1), col = c(blue, "#D9E2F3"), border = c(blue, "#8EA9DB"),
              xlab = "Probability across 200 Monte Carlo replications",
              main = "Whole-trajectory scores preserve the level signal\nBlue = detection; pale = wrong exclusion",
              cex.names = 0.92, cex.axis = 0.92, cex.lab = 0.96, col.main = navy)
axis(1, at = seq(0, 1, 0.1))
abline(v = 0.10, lty = 3, col = gray)
text(as.vector(vals) + 0.018, as.vector(bp), sprintf("%.1f%%", 100 * as.vector(vals)),
     cex = 0.80, col = navy, xpd = NA)
dev.off()

cat(sprintf("Representative draw: detection %.3f; wrong detection %.3f; clean=%d; bad=%d\n",
            mean(!screen$selected$avg[bad]), mean(!screen$selected$avg[clean]),
            sum(clean), sum(bad)))
