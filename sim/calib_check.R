# is the weighted conformal p-value calibrated? (core correctness, ref-independent)-
# under the null (CLEAN EC, correct model, weighted), the EC p-value must be
# ~super-uniform: P(p <= gamma) <= gamma (approx =gamma). test the two-sided |p|,
# upper p_hi, lower p_lo for DGP1 (binary, w_of weights) and DGP3 (continuous,
# logistic weights, post-fix). pool clean-EC p-values across reps; compare to gamma.

source("engine2.R")
args <- commandArgs(trailingOnly = TRUE)
nsim  <- if (length(args) >= 1) as.integer(args[1]) else 1500L
ncore <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores() - 2L)
set.seed(20260713L)
n_r0 <- 60L; n_ec <- 200L; Delta <- 5; pattern <- "A"
gam <- c(0.05, 0.10, 0.20)

collect <- function(dgp, dratio) {
  one <- function(r) tryCatch({
    set.seed(20260713L + r + dgp * 7919L)
    dat <- gen_data(dgp, 2L, n_r0, n_ec, pattern, Delta, 0.0, FALSE)   # piB=0: ALL clean
    P <- screen_tails(dgp, dat, "c2", "internal", TRUE, "correct", dratio = dratio)
    cbind(P$abs, P$hi, P$lo)                                            # all ECs are clean
  }, error = function(e) NULL)
  R <- parallel::mclapply(seq_len(nsim), one, mc.cores = ncore)
  R <- R[vapply(R, function(x) is.matrix(x) && ncol(x) == 3, TRUE)]
  M <- do.call(rbind, R)
  t(vapply(gam, function(g) c(abs = mean(M[, 1] <= g), hi = mean(M[, 2] <= g), lo = mean(M[, 3] <= g)), numeric(3)))
}

for (cfg in list(list(1L, "oracle", "DGP1 binary (w_of exact weights)"),
                 list(3L, "oracle", "DGP3 continuous (oracle density ratio)"),
                 list(3L, "logistic", "DGP3 continuous (logistic density ratio, post-fix)"))) {
  tab <- collect(cfg[[1]], cfg[[2]])
  cat(sprintf("\n== %s ==\n", cfg[[3]]))
  out <- data.frame(gamma = gam, P_abs_le = round(tab[, 1], 3), P_hi_le = round(tab[, 2], 3), P_lo_le = round(tab[, 3], 3))
  cat("P(p <= gamma) should be ~<= gamma (calibrated/conservative):\n"); print(out, row.names = FALSE)
}
