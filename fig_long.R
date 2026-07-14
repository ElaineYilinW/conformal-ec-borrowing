# fig_long.R  ->  fig_long.pdf
# Why screening must be longitudinal, and why the score must match the drift shape.
# Three incompatibility patterns; RCT controls (gray) follow a linear natural progression,
# a few external controls (amber) depart from it in different ways. Each pattern is caught
# by a DIFFERENT estimand/contrast -- motivating the estimand-aligned longitudinal score.
set.seed(7L)
gray <- "#8A8880"; amber <- "#E0891E"; ink <- "#2B2B2B"
tr <- function(c, a) adjustcolor(c, alpha.f = a)
tt <- 1:4
mc  <- function(t) 2 + 0.7 * t                       # control natural progression (linear)
gen <- function(mfun, n, sdI = 0.30, sdE = 0.18) {   # subjects = random intercept + visit noise
  t(sapply(seq_len(n), function(i) mfun(tt) + rnorm(1, 0, sdI) + rnorm(4, 0, sdE)))
}
panels <- list(
  list(mi = function(t) 2 + 0.7 * t + 1.3 * pmax(t - 2, 0)^1.7,
       ttl = "Late deterioration", note = "flagged by: visit-average, slope  (not baseline)"),
  list(mi = function(t) 2 + 0.7 * t + 2.1,
       ttl = "Global level shift", note = "flagged by: baseline, visit-average  (not slope)"),
  list(mi = function(t) 2 + 0.7 * t - 2.0 * exp(-((t - 2.5)^2) / 0.5),
       ttl = "Transient dip, then recovery", note = "flagged by: trajectory-shape (Mahalanobis) only"))

pdf("fig_long.pdf", width = 10.4, height = 3.5, pointsize = 10)
par(mfrow = c(1, 3), oma = c(2.3, 0, 0.3, 0), mar = c(3.6, 3.0, 2.2, 0.8), mgp = c(1.9, 0.6, 0))
for (p in panels) {
  Yc <- gen(mc, 12); Yi <- gen(p$mi, 6)
  yl <- range(Yc, Yi) + c(-0.2, 0.4)
  plot(NA, xlim = c(0.85, 4.15), ylim = yl, xlab = "visit", ylab = "outcome Y",
       xaxt = "n", bty = "l", main = p$ttl, font.main = 1, cex.main = 1.1)
  axis(1, at = tt)
  for (i in 1:nrow(Yc)) lines(tt, Yc[i, ], col = tr(gray, 0.35), lwd = 1)
  for (i in 1:nrow(Yi)) lines(tt, Yi[i, ], col = tr(amber, 0.45), lwd = 1)
  lines(tt, mc(tt), col = gray, lwd = 3.2)
  lines(tt, p$mi(tt), col = amber, lwd = 3.2)
  mtext(p$note, side = 1, line = 2.35, cex = 0.72, col = ink)
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
legend("bottom", horiz = TRUE, bty = "n", cex = 0.9, lwd = 3,
       col = c(gray, amber), legend = c("RCT controls (natural progression)", "incompatible external controls"))
dev.off()
cat("wrote fig_long.pdf\n")
