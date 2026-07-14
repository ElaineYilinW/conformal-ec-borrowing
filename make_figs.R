# make_figs.R  ->  fig_outcome.pdf , fig_pval.pdf   (Figure 8 of the blueprint)
# Built on the REFERENCE-CELL / Table-1 setup: DGP1 hetero, n_r0=60, n_ec=200,
# piB=0.30, Delta=8, correct model, c2 contrast, weighted CV+ calibration.
# Panel 1: density of the visit-average outcome Z=c2'Y (why the raw outcome confounds
#          covariate shift with contamination). Panel 2: the weighted residual-rank
#          conformal p-values and the screen's truncation at gamma.
source("engine2.R")
set.seed(20260713L)
n_r1 <- 120L; n_r0 <- 60L; n_ec <- 200L; piB <- 0.30; Delta <- 8
dat <- gen_data(1L, n_r1, n_r0, n_ec, "A", Delta, piB, FALSE)   # controls only (effect not used by the screen)
Yc <- dat$Yc; Xc <- dat$Xc; Ye <- dat$Ye; Xe <- dat$Xe; B <- dat$B; cc <- C2
Zr <- as.numeric(Yc %*% cc); Ze <- as.numeric(Ye %*% cc)
gray <- "#888780"; blue <- "#378ADD"; amber <- "#EF9F27"; red <- "#B23B3B"; green <- "#639922"
tr <- function(c, a = 0.32) adjustcolor(c, alpha.f = a)
fc <- mean(B == 0); fk <- mean(B == 1)
dc <- function(v, sc = 1) { d <- density(v, n = 220, adjust = 1.1); d$y <- d$y * sc; d }

# ---- Panel 1: outcome Y density ----
pdf("fig_outcome.pdf", width = 6.8, height = 3.0, pointsize = 10)
par(mar = c(3.2, 0.8, 0.6, 0.6), mgp = c(1.9, 0.6, 0))
d1 <- dc(Zr); d2 <- dc(Ze[B == 0], fc); d3 <- dc(Ze[B == 1], fk)
plot(NA, xlim = range(d1$x, d2$x, d3$x), ylim = c(0, max(d1$y)),
     xlab = "outcome  Z = visit-average control level", ylab = "", yaxt = "n", bty = "n")
polygon(d3, col = tr(amber), border = amber, lwd = 1.6)
polygon(d2, col = tr(blue),  border = blue,  lwd = 1.6)
polygon(d1, col = tr(gray),  border = gray,  lwd = 1.6)
legend("topright", bty = "n", cex = 0.85,
       legend = c("RCT control", "compatible EC  (borrow)", "non-compatible EC  (drop)"),
       fill = c(tr(gray), tr(blue), tr(amber)), border = c(gray, blue, amber))
dev.off()

# ---- Panel 2: weighted residual-rank p-values + truncation at gamma ----
P  <- tails(1L, dat, "c2", "internal", TRUE, "correct", calib = "cv")
pv <- P$abs                       # symmetric weighted conformal p-value, one per EC
gam <- 0.10                       # illustrative screen threshold (adaptive gamma lands ~0.05-0.10)
br <- seq(0, 1, by = 0.05)
hc <- hist(pv[B == 0], breaks = br, plot = FALSE)   # compatible
hi <- hist(pv[B == 1], breaks = br, plot = FALSE)   # incompatible
ym <- max(hc$counts, hi$counts)
pdf("fig_pval.pdf", width = 6.8, height = 3.0, pointsize = 10)
par(mar = c(3.2, 3.0, 0.6, 0.6), mgp = c(1.9, 0.6, 0))
plot(NA, xlim = c(0, 1), ylim = c(0, ym * 1.12),
     xlab = expression("weighted residual-rank conformal " * italic(p) * "-value"),
     ylab = "count of external controls", bty = "l")
rect(0, 0, gam, ym * 1.12, col = tr(red, 0.08), border = NA)   # drop zone p <= gamma
for (i in seq_along(hi$counts)) rect(br[i], 0, br[i+1], hi$counts[i], col = tr(amber, 0.55), border = amber)
for (i in seq_along(hc$counts)) rect(br[i], 0, br[i+1], hc$counts[i], col = tr(blue, 0.40), border = blue)
abline(v = gam, lty = 2, lwd = 1.6, col = red)
text(gam, ym * 1.10, expression(gamma), col = red, pos = 4, cex = 1.0)
text(gam/2, ym * 0.60, "drop", col = red, srt = 90, cex = 0.9)
text(0.58, ym * 0.85, "keep  (p > gamma)", col = green, cex = 0.95, font = 2)
legend("right", bty = "n", cex = 0.82, inset = 0.02,
       legend = c("compatible EC", "non-compatible EC"),
       fill = c(tr(blue, 0.40), tr(amber, 0.55)), border = c(blue, amber))
dev.off()
cat(sprintf("wrote fig_outcome.pdf, fig_pval.pdf | keep-frac compatible=%.2f, incompat dropped=%.2f (gamma=%.2f)\n",
            mean(pv[B==0] > gam), mean(pv[B==1] <= gam), gam))
