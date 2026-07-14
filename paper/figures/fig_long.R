# fig_long.R  ->  fig_long.pdf   (Figure 1 of the blueprint)
# Why screening must be longitudinal, and why the score must match the drift shape.
# Clean version: the RCT controls are shown as a shaded "compatible band" (not spaghetti).
# MOST external controls are compatible (blue, tracking the band); a FEW are incompatible
# (amber, drifting out) -- and the SHAPE of the drift decides which estimand catches it.
set.seed(7L)
gray <- "#8A8880"; blue <- "#2F7DC8"; amber <- "#E0891E"; ink <- "#2B2B2B"
tr <- function(c, a) adjustcolor(c, alpha.f = a)
tt <- 1:4
mc  <- function(t) 2 + 0.7 * t          # control natural progression
sdc <- 0.42
gen_comp <- function(n) t(sapply(seq_len(n), function(i) mc(tt) + rnorm(1, 0, 0.22) + rnorm(4, 0, 0.11)))
gen_inc  <- function(mi, n) t(sapply(seq_len(n), function(i) mi(tt) + rnorm(1, 0, 0.16) + rnorm(4, 0, 0.11)))
panels <- list(
  list(mi = function(t) 2 + 0.7 * t + 1.35 * pmax(t - 2, 0)^1.7,
       ttl = "Late deterioration", note = "flagged by visit-average & slope (not baseline)"),
  list(mi = function(t) 2 + 0.7 * t + 2.0,
       ttl = "Global level shift", note = "flagged by baseline & visit-average (not slope)"),
  list(mi = function(t) 2 + 0.7 * t - 2.0 * exp(-((t - 2.5)^2) / 0.5),
       ttl = "Transient dip, then recovery", note = "flagged by trajectory-shape (Mahalanobis) only"))

pdf("fig_long.pdf", width = 10.4, height = 3.75, pointsize = 10)
par(mfrow = c(1, 3), oma = c(2.5, 0, 0.3, 0), mar = c(3.4, 3.0, 2.2, 0.8), mgp = c(1.9, 0.6, 0))
for (p in panels) {
  comp <- gen_comp(4); inc <- gen_inc(p$mi, 2)
  up <- mc(tt) + 1.7 * sdc; lo <- mc(tt) - 1.7 * sdc
  yl <- range(lo, up, comp, inc) + c(-0.2, 0.4)
  plot(NA, xlim = c(0.85, 4.15), ylim = yl, xlab = "visit", ylab = "outcome Y",
       xaxt = "n", bty = "l", main = p$ttl, font.main = 1, cex.main = 1.1)
  axis(1, at = tt)
  polygon(c(tt, rev(tt)), c(up, rev(lo)), col = tr(gray, 0.16), border = NA)  # compatible band
  lines(tt, mc(tt), col = gray, lwd = 2.6)                                    # control mean
  for (i in seq_len(nrow(comp))) lines(tt, comp[i, ], col = tr(blue, 0.7), lwd = 1.2)  # compatible ECs (most)
  for (i in seq_len(nrow(inc)))  lines(tt, inc[i, ],  col = amber, lwd = 2.3)          # incompatible ECs (few)
  mtext(p$note, side = 1, line = 2.3, cex = 0.72, col = ink)
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
legend("bottom", horiz = TRUE, bty = "n", cex = 0.9, seg.len = 1.7,
       legend = c("RCT control range", "compatible EC  (most)", "incompatible EC  (few)"),
       col = c(tr(gray, 0.5), tr(blue, 0.9), amber), lwd = c(7, 2, 2.3))
dev.off()
cat("wrote fig_long.pdf\n")
