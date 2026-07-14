# robustness.R  ->  fills the robustness table (tab:robust) of borrowing_summary.tex
# ---------------------------------------------------------------------------
# One-axis-at-a-time stress test off the sym-ada / CV+ baseline. Every row keeps
# the base cell (hetero base DGP, n_r0=60, n_ec=200, piB=0.30, Delta=8, pattern A,
# contrast-c2 score, internal covariance, weighted, correct model, CV+ calibration,
# symmetric adaptive screen) and changes exactly ONE thing:
#   1 baseline                       (nothing changed)
#   2 severe model misspecification  model = severe  (drop covariate from working mean)
#   3 unweighted                     weighted = FALSE (drop covariate-shift weighting)
#   4 split instead of CV+           calib = split    (one-time train/calibrate split)
#   5 C-slope drift, c2              pattern = C      (ramp drift; visit-average score)
#   6 C-slope drift, Mahalanobis     pattern = C, score = quad
#   7 C-slope drift, max (trap)      pattern = C, score = max
# Reported per row = the SCREEN (sym-ada) under that setting; rct/oracle also written
# to the CSV for reference. Metrics: bias, rmse, type-I, power, detection, false-excl.
# Consolidates the former newsetup6/7/8 exploratory scripts into one reproducible run.
#
# Usage:  Rscript robustness.R [nsim] [Bboot] [ncore] [Delta]
#   defaults: nsim=300  Bboot=500  ncore=detectCores()-2  Delta=8
source("engine2.R")
args  <- commandArgs(trailingOnly = TRUE)
nsim  <- if (length(args) >= 1) as.integer(args[1]) else 300L
Bboot <- if (length(args) >= 2) as.integer(args[2]) else 500L
ncore <- if (length(args) >= 3) as.integer(args[3]) else max(1L, parallel::detectCores() - 2L)
Delta <- if (length(args) >= 4) as.integer(args[4]) else 8L
set.seed(20260713L)
n_r1 <- 120L; n_r0 <- 60L; n_ec <- 200L; piB <- 0.30; cc <- C2; tau0 <- sum(cc * tau_vec)
gone <- c(0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30); Ks <- paste0("s", seq_along(gone))

pick <- function(th0, bm, cols) {           # intFRT CSB de-biased MSE threshold rule
  d2 <- (th0[cols] - th0["rct"])^2
  v  <- apply(bm[, cols, drop = FALSE], 2, var, na.rm = TRUE)
  vp <- apply(bm[, cols, drop = FALSE], 2, function(col) var(col - bm[, "rct"], na.rm = TRUE))
  w  <- which.min(pmax(d2 - vp, 0) + v); cols[if (length(w)) w else 1] }
ci_rej <- function(v) { q <- tryCatch(quantile(v, c(.025, .975), na.rm = TRUE), error = function(e) c(NA, NA))
  if (any(is.na(q))) return(FALSE); !(q[1] <= 0 & 0 <= q[2]) }

one_rep <- function(r, shape, model, weighted, pattern, calib) tryCatch({
  set.seed(20260713L + r)
  cand <- function(dat) {
    Zt <- as.numeric(dat$Yt %*% cc); Zr <- as.numeric(dat$Yc %*% cc); Ze <- as.numeric(dat$Ye %*% cc)
    P <- tails(1L, dat, shape, "internal", weighted, model, calib = calib)
    e <- function(sel) est_tau_c(Zt, Zr, dat$Xc, Ze, dat$Xe, sel)
    list(v = c(rct = e(rep(FALSE, n_ec)), oracle = e(dat$B == 0),
               setNames(vapply(gone, function(g) e(P$abs > g), 0), Ks)), P = P) }
  run <- function(effect) {
    dat <- gen_data(1L, n_r1, n_r0, n_ec, pattern, Delta, piB, effect, vmode = "hetero")
    C0 <- cand(dat); th0 <- C0$v
    bm <- matrix(NA, Bboot, length(th0), dimnames = list(NULL, names(th0)))
    for (b in seq_len(Bboot)) {
      it <- sample.int(n_r1, n_r1, TRUE); ic <- sample.int(n_r0, n_r0, TRUE); ie <- sample.int(n_ec, n_ec, TRUE)
      db <- list(Yt = dat$Yt[it, , drop = FALSE], Xt = dat$Xt[it], Yc = dat$Yc[ic, , drop = FALSE],
                 Xc = dat$Xc[ic], Ye = dat$Ye[ie, , drop = FALSE], Xe = dat$Xe[ie], B = dat$B[ie])
      bm[b, ] <- cand(db)$v
    }
    sy <- pick(th0, bm, Ks); sel <- C0$P$abs > gone[as.integer(sub("s", "", sy))]
    key <- c(rct = "rct", oracle = "oracle", scr = sy)
    list(th = setNames(as.numeric(th0[key]), names(key)),
         rej = setNames(vapply(key, function(k) ci_rej(bm[, k]), TRUE), names(key)),
         det = sum(!sel & dat$B == 1), nd = sum(dat$B == 1),
         fer = sum(!sel & dat$B == 0), nf = sum(dat$B == 0))
  }
  e <- run(TRUE); n0 <- run(FALSE)
  list(th = e$th, rej_e = e$rej, rej_n = n0$rej, det = e$det, nd = e$nd, fer = e$fer, nf = e$nf)
}, error = function(e) NULL)

# grid: shape, model, weighted, pattern, calib, label  (one axis off baseline each)
grid <- list(
  list("c2",   "correct", TRUE,  "A", "cv",    "baseline (sym-ada, CV+)"),
  list("c2",   "severe",  TRUE,  "A", "cv",    "severe model misspec."),
  list("c2",   "correct", FALSE, "A", "cv",    "unweighted (correct)"),
  list("c2",   "correct", TRUE,  "A", "split", "split instead of CV+"),
  list("c2",   "correct", TRUE,  "C", "cv",    "C-slope drift, c2"),
  list("quad", "correct", TRUE,  "C", "cv",    "C-slope drift, Mahalanobis"),
  list("max",  "correct", TRUE,  "C", "cv",    "C-slope drift, max (trap)"))

t0 <- Sys.time(); rows <- list()
for (cf in grid) {
  shape <- cf[[1]]; model <- cf[[2]]; weighted <- cf[[3]]; pattern <- cf[[4]]; calib <- cf[[5]]; lab <- cf[[6]]
  res <- parallel::mclapply(seq_len(nsim), function(r) one_rep(r, shape, model, weighted, pattern, calib), mc.cores = ncore)
  res <- res[vapply(res, function(x) is.list(x) && length(x$th) == 3, TRUE)]
  TH <- t(vapply(res, function(x) x$th, numeric(3))); RE <- t(vapply(res, function(x) x$rej_e, logical(3)))
  RN <- t(vapply(res, function(x) x$rej_n, logical(3))); colnames(TH) <- colnames(RE) <- colnames(RN) <- c("rct","oracle","scr")
  det <- sum(vapply(res, function(x) x$det, 0)) / max(sum(vapply(res, function(x) x$nd, 0)), 1)
  fer <- sum(vapply(res, function(x) x$fer, 0)) / sum(vapply(res, function(x) x$nf, 0))
  for (m in c("rct","oracle","scr")) rows[[length(rows)+1]] <- data.frame(
    setting = lab, method = m,
    bias = round(mean(TH[,m]) - tau0, 3), rmse = round(sqrt(mean((TH[,m]-tau0)^2)), 3),
    typeI = round(mean(RN[,m]), 3), power = round(mean(RE[,m]), 3),
    detection = if (m=="scr") round(det,3) else if (m=="oracle") 1 else NA,
    false_excl = if (m=="scr") round(fer,3) else NA)
}
out <- do.call(rbind, rows); write.csv(out, "robustness_results.csv", row.names = FALSE)
scr <- out[out$method == "scr", c("setting","bias","rmse","typeI","power","detection","false_excl")]
cat(sprintf("completed %d reps/setting in %.0f s (Bboot=%d, Delta=%d, tau=%.4f)\n",
            nsim, as.numeric(difftime(Sys.time(), t0, "secs")), Bboot, Delta, tau0))
cat("\n== Robustness table (tab:robust) -- SCREEN (sym-ada) per setting ==\n")
print(scr, row.names = FALSE)
cat("\n(full rct/oracle/screen breakdown in robustness_results.csv)\n")
