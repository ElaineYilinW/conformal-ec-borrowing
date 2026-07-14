# reference cell -> fills Table 1----
# DGP 1, piB=0.3, pattern A, Delta=8, estimand c2, weighted CV+ calibration,
# correct model. Compares: rct / full / oracle / (quad,max at internal, one-sided
# adaptive) / (c2 at raw,source; nonsym adaptive) / (c2 internal: sym-ada, nonsym-fix,
# nonsym-ada). Metrics: detection, false-excl, bias, RMSE, type-I, power.
# lambda = aware inverse-variance (kept); adaptive gamma = intFRT CSB de-biased MSE
# (max((tau(g)-tau_rct)^2 - Var_b(paired diff), 0) + bootstrap-var). calib = CV+ (K=5).
#
# Usage:  Rscript reference_cell.R [nsim] [Bboot] [ncore]
#   defaults: nsim=300  Bboot=500  ncore=detectCores()-2
# Type-I is a full-pipeline percentile-bootstrap CI: it carries a shared ~0.05-0.07
# floor set by Bboot (see typeI_check.R); Bboot=500 keeps that floor near nominal.

source("engine2.R")
args <- commandArgs(trailingOnly = TRUE)
nsim  <- if (length(args) >= 1) as.integer(args[1]) else 300L
Bboot <- if (length(args) >= 2) as.integer(args[2]) else 500L
ncore <- if (length(args) >= 3) as.integer(args[3]) else max(1L, parallel::detectCores() - 2L)
set.seed(20260713L)

n_r1 <- 120L; n_r0 <- 60L; n_ec <- 200L; piB <- 0.30; Delta <- 8; pattern <- "A"
cc <- C2; tau0 <- sum(cc * tau_vec)

gone <- c(0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30)     # one-sided / sym grid (finer low end)
gh <- c(0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.20) # non-symmetric hi grid, down to 0.01
gl <- c(0.00, 0.02, 0.05)                               # so adaptive can pick a GENTLE cut
g2 <- expand.grid(gh = gh, gl = gl)
Kq <- paste0("quad@", gone); Km <- paste0("max@", gone); Ksym <- paste0("c2sym@", gone)
Kraw <- paste0("c2raw@", g2$gh, "_", g2$gl); Ksrc <- paste0("c2src@", g2$gh, "_", g2$gl)
Kada <- paste0("c2ada@", g2$gh, "_", g2$gl)

# 5 unique (shape,cov) screens shared across methods
screens5 <- function(dat) list(
  qi = tails(1, dat, "quad", "internal", TRUE, "correct", calib = "cv"),
  mi = tails(1, dat, "max",  "internal", TRUE, "correct", calib = "cv"),
  cr = tails(1, dat, "c2",   "raw",      TRUE, "correct", calib = "cv"),
  ci = tails(1, dat, "c2",   "internal", TRUE, "correct", calib = "cv"),
  cs = tails(1, dat, "c2",   "source",   TRUE, "correct", calib = "cv"))

cand <- function(dat, P) {
  e <- function(sel) tau_of(1, dat, cc, sel)
  v <- c(rct = e(rep(FALSE, n_ec)), full = e(rep(TRUE, n_ec)), oracle = e(dat$B == 0),
         setNames(vapply(gone, function(g) e(P$qi$hi > g), 0), Kq),
         setNames(vapply(gone, function(g) e(P$mi$hi > g), 0), Km),
         setNames(vapply(gone, function(g) e(P$ci$abs > g), 0), Ksym),
         setNames(vapply(seq_len(nrow(g2)), function(i) e((P$cr$hi > g2$gh[i]) & (P$cr$lo > g2$gl[i])), 0), Kraw),
         setNames(vapply(seq_len(nrow(g2)), function(i) e((P$cs$hi > g2$gh[i]) & (P$cs$lo > g2$gl[i])), 0), Ksrc),
         setNames(vapply(seq_len(nrow(g2)), function(i) e((P$ci$hi > g2$gh[i]) & (P$ci$lo > g2$gl[i])), 0), Kada),
         c2fix = e((P$ci$hi > 0.30) & (P$ci$lo > 0.05)))
  v
}
pick <- function(th0, bm, cols) { d2 <- (th0[cols] - th0["rct"])^2   # intFRT CSB de-biased MSE
  v <- apply(bm[, cols, drop = FALSE], 2, var, na.rm = TRUE)
  vp <- apply(bm[, cols, drop = FALSE], 2, function(col) var(col - bm[, "rct"], na.rm = TRUE))
  w <- which.min(pmax(d2 - vp, 0) + v); cols[if (length(w)) w else 1] }

rep_meth <- c("rct", "full", "oracle", "quad", "max", "c2raw", "c2src", "c2sym", "c2fix", "c2ada")
screen_meth <- c("oracle", "quad", "max", "c2raw", "c2src", "c2sym", "c2fix", "c2ada")

one_rep <- function(r) {
  set.seed(20260713L + r)
  de <- gen_data(1, n_r1, n_r0, n_ec, pattern, Delta, piB, TRUE)
  dn <- gen_data(1, n_r1, n_r0, n_ec, pattern, Delta, piB, FALSE)
  run <- function(dat) {
    P <- screens5(dat); th0 <- cand(dat, P)
    bm <- matrix(NA, Bboot, length(th0), dimnames = list(NULL, names(th0)))
    for (b in seq_len(Bboot)) {
      it <- sample.int(n_r1, n_r1, TRUE); ic <- sample.int(n_r0, n_r0, TRUE); ie <- sample.int(n_ec, n_ec, TRUE)
      db <- list(Yt = dat$Yt[it, , drop = FALSE], Xt = dat$Xt[it], Yc = dat$Yc[ic, , drop = FALSE],
                 Xc = dat$Xc[ic], Ye = dat$Ye[ie, , drop = FALSE], Xe = dat$Xe[ie], B = dat$B[ie])
      bm[b, ] <- cand(db, screens5(db))
    }
    qs <- pick(th0, bm, Kq); ms <- pick(th0, bm, Km); sy <- pick(th0, bm, Ksym)
    rw <- pick(th0, bm, Kraw); sr <- pick(th0, bm, Ksrc); ad <- pick(th0, bm, Kada)
    key <- c(rct = "rct", full = "full", oracle = "oracle", quad = qs, max = ms,
             c2raw = rw, c2src = sr, c2sym = sy, c2fix = "c2fix", c2ada = ad)
    th <- setNames(as.numeric(th0[key]), names(key))
    rej <- setNames(vapply(key, function(k) { q <- quantile(bm[, k], c(.025, .975), na.rm = TRUE); !(q[1] <= 0 & 0 <= q[2]) }, TRUE), names(key))
    # selections for detection/FER (recompute on point p-values at chosen gamma)
    gg <- function(k) as.numeric(strsplit(sub(".*@", "", k), "_")[[1]])
    sel <- list(oracle = dat$B == 0,
                quad = P$qi$hi > gg(qs), max = P$mi$hi > gg(ms),
                c2raw = { g <- gg(rw); (P$cr$hi > g[1]) & (P$cr$lo > g[2]) },
                c2src = { g <- gg(sr); (P$cs$hi > g[1]) & (P$cs$lo > g[2]) },
                c2sym = P$ci$abs > gg(sy),
                c2fix = (P$ci$hi > 0.30) & (P$ci$lo > 0.05),
                c2ada = { g <- gg(ad); (P$ci$hi > g[1]) & (P$ci$lo > g[2]) })
    list(th = th, rej = rej, sel = sel)
  }
  re <- run(de); rn <- run(dn)
  det <- fer <- setNames(numeric(length(screen_meth)), screen_meth)
  for (m in screen_meth) { ex <- !re$sel[[m]]; det[m] <- sum(ex & de$B == 1); fer[m] <- sum(ex & de$B == 0) }
  list(th = re$th, rej_e = re$rej, rej_n = rn$rej, det = det, fer = fer,
       nd = sum(de$B == 1), nf = sum(de$B == 0))
}

t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nsim), one_rep, mc.cores = ncore)
ok <- !vapply(res, is.null, TRUE) & vapply(res, function(x) !inherits(x, "try-error"), TRUE)
res <- res[ok]; cat(sprintf("completed %d/%d reps in %.0f s\n", length(res), nsim, as.numeric(difftime(Sys.time(), t0, "secs"))))

TH <- t(vapply(res, function(x) x$th, numeric(length(rep_meth))))
RE <- t(vapply(res, function(x) x$rej_e, logical(length(rep_meth))))
RN <- t(vapply(res, function(x) x$rej_n, logical(length(rep_meth))))
colnames(TH) <- colnames(RE) <- colnames(RN) <- rep_meth
det <- rowSums(vapply(res, function(x) x$det, numeric(length(screen_meth)))); names(det) <- screen_meth
fer <- rowSums(vapply(res, function(x) x$fer, numeric(length(screen_meth)))); names(fer) <- screen_meth
nd <- sum(vapply(res, function(x) x$nd, 0)); nf <- sum(vapply(res, function(x) x$nf, 0))

row <- function(m) data.frame(method = m,
  detection = if (m %in% screen_meth) round(det[m] / nd, 3) else NA,
  false_excl = if (m %in% screen_meth) round(fer[m] / nf, 3) else NA,
  bias = round(mean(TH[, m]) - tau0, 3), rmse = round(sqrt(mean((TH[, m] - tau0)^2)), 3),
  typeI = round(mean(RN[, m]), 3), power = round(mean(RE[, m]), 3))
out <- do.call(rbind, lapply(rep_meth, row))
write.csv(out, "reference_cell_results.csv", row.names = FALSE)
cat("\n== Reference cell (DGP1, piB=0.3, A, Delta=8, c2); tau=", round(tau0, 4), " ==\n", sep = "")
print(out, row.names = FALSE)
