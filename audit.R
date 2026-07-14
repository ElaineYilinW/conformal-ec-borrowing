# THOROUGH AUDIT: is the non-symmetric screen bias a code bug or the selection?--
# three independent checks, each falsifiable:
#  CHECK 1  estimand+estimator sanity: huge CLEAN sample, no screen -> rct & full
#           must both equal tau0 (else DGP/estimand/est_tau_c is broken).
#  CHECK 2  ZERO contamination (piB=0, EVERY EC clean): screen bias.
#           if non-symmetric is STILL biased with NOTHING to delete correctly,
#           the bias IS the asymmetric residual truncation, not a bug and not
#           "removing contamination wrong". symmetric must stay ~0.
#  CHECK 3  piB=0.3: the full picture. full<0 (keeps contamination), oracle~0,
#           sym~0, nonsym>0 (over-cuts high tail). signs must be coherent.
# DGP1 A Delta=5 c2 weighted split correct.  (falsification, not confirmation.)

source("engine2.R")
args <- commandArgs(trailingOnly = TRUE)
nsim  <- if (length(args) >= 1) as.integer(args[1]) else 1000L
ncore <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores() - 2L)
set.seed(20260713L)
n_r1 <- 120L; n_r0 <- 60L; n_ec <- 200L; Delta <- 5; pattern <- "A"; cc <- C2
tau0 <- sum(cc * tau_vec); ghi <- 0.20; glo <- 0.05; gsym <- 0.10

# ---- CHECK 1: estimand + estimator on a huge CLEAN sample (no screen) ---------
set.seed(1L)
big <- gen_data(1L, 20000L, 20000L, 20000L, pattern, Delta, 0.0, TRUE)   # piB=0, all clean
Ztb <- as.numeric(big$Yt %*% cc); Zrb <- as.numeric(big$Yc %*% cc); Zeb <- as.numeric(big$Ye %*% cc)
tau_full <- est_tau_c(Ztb, Zrb, big$Xc, Zeb, big$Xe, rep(TRUE,  length(Zeb)))
tau_rct  <- est_tau_c(Ztb, Zrb, big$Xc, Zeb, big$Xe, rep(FALSE, length(Zeb)))
cat(sprintf("\nCHECK 1  huge clean n=20000:  tau0=%.4f   rct=%.4f   full(borrow all clean)=%.4f\n",
            tau0, tau_rct, tau_full))
cat("         -> both must be ~= tau0; if so DGP + estimand + est_tau_c are correct.\n")

# ---- CHECKS 2 & 3: screen bias at piB=0 (all clean) and piB=0.3 ---------------
one <- function(r, piB) tryCatch({
  set.seed(20260713L + r + as.integer(piB * 1e4))
  dat <- gen_data(1L, n_r1, n_r0, n_ec, pattern, Delta, piB, TRUE)
  P <- screen_tails(1L, dat, "c2", "internal", TRUE, "correct")
  Zt <- as.numeric(dat$Yt %*% cc); Zr <- as.numeric(dat$Yc %*% cc); Ze <- as.numeric(dat$Ye %*% cc)
  e <- function(sel) est_tau_c(Zt, Zr, dat$Xc, Ze, dat$Xe, sel)
  ns <- (P$hi > ghi) & (P$lo > glo); sy <- P$abs > gsym; cl <- dat$B == 0
  c(rct = e(rep(FALSE, n_ec)), full = e(rep(TRUE, n_ec)), oracle = e(cl),
    sym = e(sy), nonsym = e(ns),
    fer_sym = if (any(cl)) mean(!sy[cl]) else NA, fer_ns = if (any(cl)) mean(!ns[cl]) else NA)
}, error = function(e) NULL)

for (piB in c(0.0, 0.3)) {
  res <- parallel::mclapply(seq_len(nsim), function(r) one(r, piB), mc.cores = ncore)
  res <- res[vapply(res, function(x) is.numeric(x) && length(x) == 7, TRUE)]
  M <- do.call(rbind, res)
  b <- colMeans(M[, 1:5]) - tau0
  tab <- data.frame(method = c("rct","full","oracle","sym","nonsym"),
                    bias = round(b, 3),
                    false_excl_clean = c(NA, NA, NA, round(mean(M[,"fer_sym"], na.rm=TRUE), 3),
                                         round(mean(M[,"fer_ns"], na.rm=TRUE), 3)))
  cat(sprintf("\nCHECK piB=%.1f  (%d reps):%s\n", piB, nrow(M),
              if (piB == 0) "  NO contamination -> any nonsym bias is pure asymmetric truncation" else ""))
  print(tab, row.names = FALSE)
}
cat("\nverdict: if CHECK1 ~0, CHECK2 nonsym>0 with sym~0, CHECK3 signs coherent -> no bug; bias is the selection.\n")
