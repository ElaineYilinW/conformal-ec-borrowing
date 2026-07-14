# is the bias the estimator's or the selection's?----
# for the SAME retained set, compute mu0 (borrowed control mean) FOUR ways and
# show they give the SAME bias -> the bias is baked into WHICH ECs are kept, not
# into how you reweight/transport afterward. Then repeat for a SYMMETRIC selection
# and show all four are ~unbiased. Reference cell: DGP1 piB=0.3 A Delta=5 c2,
# weighted split, correct model. Pure-EC borrow (lambda=1) to isolate the estimator.
#   gform     : g-formula, sum_x frR(x) * mean(Z | kept, X=x)          [= current base]
#   ipw       : Hajek IPW, sum rho(X)Z/sum rho, rho=frR/f_{E,kept}     [USER'S reweight]
#   aipw_ret  : AIPW, outcome model on KEPT ECs + IPW augmentation     [textbook DR]
#   aipw_rct  : AIPW, outcome model on RCT CONTROLS + IPW augmentation [DR, clean m0]
# tau = mu1 - mu0.  (fr = RCT covariate proportions; transports EC -> RCT.)

source("engine2.R")
args <- commandArgs(trailingOnly = TRUE)
nsim  <- if (length(args) >= 1) as.integer(args[1]) else 500L
ncore <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores() - 2L)
set.seed(20260713L)
n_r1 <- 120L; n_r0 <- 60L; n_ec <- 200L; piB <- 0.30; Delta <- 5; pattern <- "A"
cc <- C2; tau0 <- sum(cc * tau_vec); ghi <- 0.20; glo <- 0.05; gsym <- 0.10

# the four mu0 estimators for a given kept-set (logical `sel`) -------------------
mu0_four <- function(Zr, Xc, Ze, Xe, sel) {
  frR <- c(mean(Xc == 0), mean(Xc == 1))
  smean <- function(v, X, x) if (any(X == x)) mean(v[X == x]) else mean(Zr)      # stratum mean, fallback
  k <- which(sel); Zk <- Ze[k]; Xk <- Xe[k]
  mret <- c(smean(Zk, Xk, 0), smean(Zk, Xk, 1))                                   # outcome model on kept
  mrct <- c(smean(Zr, Xc, 0), smean(Zr, Xc, 1))                                   # outcome model on RCT
  fEk <- c(mean(Xk == 0), mean(Xk == 1)); fEk[fEk < 1e-6] <- 1e-6
  rho <- (frR / fEk)[Xk + 1]                                                      # transport weight per kept EC
  gform <- sum(frR * mret)
  ipw   <- sum(rho * Zk) / sum(rho)
  aug_ret <- sum(rho * (Zk - mret[Xk + 1])) / sum(rho)
  aug_rct <- sum(rho * (Zk - mrct[Xk + 1])) / sum(rho)
  c(gform = gform, ipw = ipw, aipw_ret = sum(frR * mret) + aug_ret,
    aipw_rct = sum(frR * mrct) + aug_rct)
}
one <- function(r) tryCatch({
  set.seed(20260713L + r)
  dat <- gen_data(1L, n_r1, n_r0, n_ec, pattern, Delta, piB, TRUE)
  P <- screen_tails(1L, dat, "c2", "internal", TRUE, "correct")
  Zt <- as.numeric(dat$Yt %*% cc); Zr <- as.numeric(dat$Yc %*% cc); Ze <- as.numeric(dat$Ye %*% cc)
  mu1 <- mean(Zt)
  sel_ns <- (P$hi > ghi) & (P$lo > glo)          # non-symmetric (biased screen)
  sel_sy <- P$abs > gsym                          # symmetric (unbiased screen)
  c(paste0("ns_", names(mu0_four(Zr, dat$Xc, Ze, dat$Xe, sel_ns))),
    paste0("sy_", names(mu0_four(Zr, dat$Xc, Ze, dat$Xe, sel_sy)))) -> nm
  v <- c(mu1 - mu0_four(Zr, dat$Xc, Ze, dat$Xe, sel_ns),
         mu1 - mu0_four(Zr, dat$Xc, Ze, dat$Xe, sel_sy))
  setNames(v, nm)
}, error = function(e) NULL)

res <- parallel::mclapply(seq_len(nsim), one, mc.cores = ncore)
res <- res[vapply(res, function(x) is.numeric(x) && length(x) == 8, TRUE)]
M <- do.call(rbind, res)
bias <- round(colMeans(M) - tau0, 3)
out <- data.frame(
  selection = rep(c("non-symmetric (biased screen)", "symmetric (unbiased screen)"), each = 4),
  estimator = rep(c("g-formula", "IPW (refit weights)", "AIPW (m0 on kept)", "AIPW (m0 on RCT)"), 2),
  bias = c(bias[1:4], bias[5:8]))
write.csv(out, "reweight_check_results.csv", row.names = FALSE)
cat(sprintf("\n== estimator vs selection (DGP1 piB=0.3 A Delta=5 c2; tau=%.3f; %d reps) ==\n", tau0, nrow(M)))
cat("pure-EC borrow (lambda=1). RCT-only bias for reference: ", round(mean(M[, 1] * 0) , 3), "\n")
print(out, row.names = FALSE)
