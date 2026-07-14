# Isolate the RCT-only type-I: is the 0.077 a borrowing effect or a bootstrap-CI floor?
# Pure randomized difference of means on the c2 projection, null DGP (effect=FALSE).
# Compare three CIs at the SAME estimator: percentile boot B=100, B=999, and normal-theory (Welch).
source("engine2.R")
set.seed(20260713L)
n_r1 <- 120L; n_r0 <- 60L; n_ec <- 200L; piB <- 0.30; Delta <- 8; pattern <- "A"; cc <- C2
nsim <- 2000L
one <- function(r) {
  set.seed(20260713L + r)
  d <- gen_data(1, n_r1, n_r0, n_ec, pattern, Delta, piB, FALSE)   # NULL: true tau=0
  Zt <- as.numeric(d$Yt %*% cc); Zr <- as.numeric(d$Yc %*% cc)
  est <- mean(Zt) - mean(Zr)
  # normal-theory Welch CI
  se <- sqrt(var(Zt)/length(Zt) + var(Zr)/length(Zr))
  rej_norm <- abs(est) > 1.96 * se
  # percentile bootstrap, B=100 and B=999 (same draws, nested)
  bd <- numeric(999)
  for (b in seq_len(999)) {
    it <- sample.int(n_r1, n_r1, TRUE); ic <- sample.int(n_r0, n_r0, TRUE)
    bd[b] <- mean(Zt[it]) - mean(Zr[ic])
  }
  q100 <- quantile(bd[1:100], c(.025, .975)); q999 <- quantile(bd, c(.025, .975))
  c(rej_b100 = !(q100[1] <= 0 & 0 <= q100[2]),
    rej_b999 = !(q999[1] <= 0 & 0 <= q999[2]),
    rej_norm = rej_norm)
}
R <- t(vapply(seq_len(nsim), one, numeric(3)))
cat(sprintf("RCT-only difference-of-means, NULL DGP (true tau=0), %d reps\n", nsim))
cat(sprintf("  percentile bootstrap B=100 : type-I = %.3f  (SE %.3f)\n", mean(R[,1]), sqrt(mean(R[,1])*(1-mean(R[,1]))/nsim)))
cat(sprintf("  percentile bootstrap B=999 : type-I = %.3f  (SE %.3f)\n", mean(R[,2]), sqrt(mean(R[,2])*(1-mean(R[,2]))/nsim)))
cat(sprintf("  normal-theory Welch CI     : type-I = %.3f  (SE %.3f)\n", mean(R[,3]), sqrt(mean(R[,3])*(1-mean(R[,3]))/nsim)))
