# full simulation harness (all blueprint setups)----
# Block DIAG : screening quality on DGP1 -- shape x cov x weight x calib x model
#              x (null + pattern x Delta). metrics: false-excl, detection.
# Block EST  : method comparison (rct/full/oracle/quad/max/aligned{one,sym,fix,nonsym})
#              x estimand x piB x (null + pattern x Delta), on DGP1(hetero),
#              DGP2(homo,source), DGP3(two-cov; weighted x dratio). metrics:
#              detection, false-excl, bias, rmse, type-I, power, coverage.
# Adaptive gamma = argmin over grid of (tau(g)-tau_rct)^2 + bootstrap variance.
# parallel over configs via mclapply (mc.cores = detectCores()-2 by default).
# smoke test: Rscript run_all.R 5 10 6

source("engine2.R")
args  <- commandArgs(trailingOnly = TRUE)
nsim  <- if (length(args) >= 1) as.integer(args[1]) else 5L
Bboot <- if (length(args) >= 2) as.integer(args[2]) else 10L
ncore <- if (length(args) >= 3) as.integer(args[3]) else max(1L, parallel::detectCores() - 2L)
set.seed(20260713L)
n_r1 <- 120L; n_r0 <- 60L; n_ec <- 200L
gone <- c(0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30)
gh <- c(0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.20); gl <- c(0.00, 0.02, 0.05); g2 <- expand.grid(gh = gh, gl = gl)
pick <- function(th0, bm, cols) { d2 <- (th0[cols] - th0["rct"])^2   # intFRT CSB de-biased MSE
  v <- apply(bm[, cols, drop = FALSE], 2, var, na.rm = TRUE)
  vp <- apply(bm[, cols, drop = FALSE], 2, function(col) var(col - bm[, "rct"], na.rm = TRUE))
  w <- which.min(pmax(d2 - vp, 0) + v); cols[if (length(w)) w else 1] }

## ==== BLOCK DIAG (DGP1 screening quality) ====================================
scen <- rbind(data.frame(pattern = "A", Delta = 0),
              expand.grid(pattern = c("A", "B", "C", "D"), Delta = c(1, 3, 5, 8), stringsAsFactors = FALSE))
diag_grid <- expand.grid(shape = c("quad", "max", "c1", "c2"), cov = c("raw", "internal", "source"),
                         weighted = c(TRUE, FALSE), calib = c("split", "cv"),
                         model = c("correct", "moderate", "severe"), si = seq_len(nrow(scen)),
                         stringsAsFactors = FALSE)
diag_one <- function(j) tryCatch({
  cfg <- diag_grid[j, ]; pat <- scen$pattern[cfg$si]; Del <- scen$Delta[cfg$si]
  dc <- da <- nc <- na <- 0
  for (r in seq_len(nsim)) {
    dat <- gen_data(1L, 2L, n_r0, n_ec, pat, Del, 0.5, FALSE)
    P <- tails(1L, dat, cfg$shape, cfg$cov, cfg$weighted, cfg$model, cfg$calib)
    del <- if (cfg$shape %in% c("quad", "max")) (P$hi <= 0.10) else (P$abs <= 0.10)
    active <- if (pat == "D") (dat$B == 1 & dat$Xe == 1) else (dat$B == 1)
    dc <- dc + sum(del & !active); nc <- nc + sum(!active)
    da <- da + sum(del & active);  na <- na + sum(active)
  }
  data.frame(block = "diag", dgp = 1L, vmode = "hetero", estimand = NA, piB = 0.5,
             shape = cfg$shape, cov = cfg$cov, weighted = cfg$weighted, calib = cfg$calib,
             model = cfg$model, dratio = NA, pattern = pat, Delta = Del, method = "screen",
             detection = if (na > 0) round(da / na, 4) else NA, false_excl = round(dc / max(nc, 1), 4),
             bias = NA, rmse = NA, typeI = NA, power = NA, coverage = NA)
}, error = function(e) NULL)

## ==== BLOCK EST (method comparison) =========================================
est_cells <- rbind(
  expand.grid(dgp = 1L, vmode = "hetero", weighted = TRUE, dratio = "none",
              estimand = c("c1", "c2"), piB = c(0.15, 0.30, 0.50), si = seq_len(nrow(scen)),
              model = "correct", stringsAsFactors = FALSE),
  expand.grid(dgp = 2L, vmode = c("homo", "source"), weighted = TRUE, dratio = "none",
              estimand = c("c1", "c2"), piB = c(0.15, 0.30, 0.50), si = seq_len(nrow(scen)),
              model = "correct", stringsAsFactors = FALSE),
  expand.grid(dgp = 3L, vmode = "homo", weighted = c(TRUE, FALSE), dratio = c("oracle", "logistic"),
              estimand = c("c1", "c2"), piB = c(0.30), si = seq_len(nrow(scen)),
              model = c("correct", "severe"), stringsAsFactors = FALSE))
reps_est <- c("rct", "full", "oracle", "quad", "max", "one", "sym", "fix", "nonsym")
scr_est  <- c("oracle", "quad", "max", "one", "sym", "fix", "nonsym")

est_one <- function(j) tryCatch({
  cfg <- est_cells[j, ]; dgp <- cfg$dgp; vmode <- cfg$vmode; level <- cfg$model
  wt <- cfg$weighted; pat <- scen$pattern[cfg$si]; Del <- scen$Delta[cfg$si]; piB <- cfg$piB
  ash <- if (cfg$estimand == "c1") "c1" else "c2"; cc <- if (cfg$estimand == "c1") C1 else C2
  tau0 <- sum(cc * tau_vec)
  Kq <- paste0("quad@", gone); Km <- paste0("max@", gone); Ko <- paste0("one@", gone)
  Ks <- paste0("sym@", gone); Ka <- paste0("nonsym@", g2$gh, "_", g2$gl)
  dr <- cfg$dratio                                       # DGP3 transport ratio: refit on retained inside est_tau_dgp3
  screens <- function(dat) list(q = tails(dgp, dat, "quad", "internal", wt, level, dratio = dr),
                                m = tails(dgp, dat, "max", "internal", wt, level, dratio = dr),
                                a = tails(dgp, dat, ash, "internal", wt, level, dratio = dr))
  cand <- function(dat, P) {
    e <- function(sel) tau_of(dgp, dat, cc, sel, dr)
    ne <- if (dgp == 3) length(dat$Xbe) else length(dat$Xe); B <- dat$B
    c(rct = e(rep(FALSE, ne)), full = e(rep(TRUE, ne)), oracle = e(B == 0),
      setNames(vapply(gone, function(g) e(P$q$hi > g), 0), Kq),
      setNames(vapply(gone, function(g) e(P$m$hi > g), 0), Km),
      setNames(vapply(gone, function(g) e(P$a$hi > g), 0), Ko),
      setNames(vapply(gone, function(g) e(P$a$abs > g), 0), Ks),
      setNames(vapply(seq_len(nrow(g2)), function(i) e((P$a$hi > g2$gh[i]) & (P$a$lo > g2$gl[i])), 0), Ka),
      fix = e((P$a$hi > 0.30) & (P$a$lo > 0.05)))
  }
  run <- function(effect) {
    dat <- gen_data(dgp, n_r1, n_r0, n_ec, pat, Del, piB, effect, vmode = vmode)
    P <- screens(dat); th0 <- cand(dat, P)
    bm <- matrix(NA, Bboot, length(th0), dimnames = list(NULL, names(th0)))
    for (b in seq_len(Bboot)) {
      ie <- sample.int(n_ec, n_ec, TRUE); ic <- sample.int(n_r0, n_r0, TRUE); it <- sample.int(n_r1, n_r1, TRUE)
      db <- if (dgp == 3)
        list(Yt = dat$Yt[it, , drop = FALSE], Xbt = dat$Xbt[it], Xct = dat$Xct[it],
             Yc = dat$Yc[ic, , drop = FALSE], Xbc = dat$Xbc[ic], Xcc = dat$Xcc[ic],
             Ye = dat$Ye[ie, , drop = FALSE], Xbe = dat$Xbe[ie], Xce = dat$Xce[ie], B = dat$B[ie])
      else list(Yt = dat$Yt[it, , drop = FALSE], Xt = dat$Xt[it], Yc = dat$Yc[ic, , drop = FALSE],
                Xc = dat$Xc[ic], Ye = dat$Ye[ie, , drop = FALSE], Xe = dat$Xe[ie], B = dat$B[ie])
      bm[b, ] <- cand(db, screens(db))
    }
    qs <- pick(th0, bm, Kq); ms <- pick(th0, bm, Km); os <- pick(th0, bm, Ko)
    ss <- pick(th0, bm, Ks); as_ <- pick(th0, bm, Ka)
    key <- c(rct = "rct", full = "full", oracle = "oracle", quad = qs, max = ms,
             one = os, sym = ss, fix = "fix", nonsym = as_)
    th <- setNames(as.numeric(th0[key]), names(key))
    rej <- setNames(vapply(key, function(k) { q <- quantile(bm[, k], c(.025, .975), na.rm = TRUE); !(q[1] <= 0 & 0 <= q[2]) }, TRUE), names(key))
    cov <- setNames(vapply(key, function(k) { q <- quantile(bm[, k], c(.025, .975), na.rm = TRUE); q[1] <= tau0 & tau0 <= q[2] }, TRUE), names(key))
    gg <- function(k) as.numeric(strsplit(sub(".*@", "", k), "_")[[1]])
    Be <- dat$B; Xe_ <- if (dgp == 3) dat$Xbe else dat$Xe
    active <- if (pat == "D") (Be == 1 & Xe_ == 1) else (Be == 1)
    sell <- list(oracle = Be == 0, quad = P$q$hi > gg(qs), max = P$m$hi > gg(ms),
                 one = P$a$hi > gg(os), sym = P$a$abs > gg(ss),
                 fix = (P$a$hi > 0.30) & (P$a$lo > 0.05),
                 nonsym = { g <- gg(as_); (P$a$hi > g[1]) & (P$a$lo > g[2]) })
    dd <- ff <- setNames(numeric(length(scr_est)), scr_est)
    for (m in scr_est) { ex <- !sell[[m]]; dd[m] <- sum(ex & active); ff[m] <- sum(ex & !active) }
    list(th = th, rej = rej, cov = cov, dd = dd, ff = ff, na = sum(active), nf = sum(!active))
  }
  acc_th <- matrix(NA, nsim, length(reps_est), dimnames = list(NULL, reps_est))
  acc_re <- acc_rn <- acc_cv <- matrix(0, nsim, length(reps_est), dimnames = list(NULL, reps_est))
  dd <- ff <- setNames(numeric(length(scr_est)), scr_est); na <- nf <- 0
  for (r in seq_len(nsim)) {
    e <- run(TRUE); n0 <- run(FALSE)
    acc_th[r, ] <- e$th; acc_re[r, ] <- e$rej; acc_cv[r, ] <- e$cov; acc_rn[r, ] <- n0$rej
    dd <- dd + e$dd; ff <- ff + e$ff; na <- na + e$na; nf <- nf + e$nf
  }
  do.call(rbind, lapply(reps_est, function(m) data.frame(
    block = "est", dgp = dgp, vmode = vmode, estimand = cfg$estimand, piB = piB,
    shape = NA, cov = "internal", weighted = wt, calib = "split", model = level, dratio = cfg$dratio,
    pattern = pat, Delta = Del, method = m,
    detection = if (m %in% scr_est && na > 0) round(dd[m] / na, 4) else NA,
    false_excl = if (m %in% scr_est) round(ff[m] / max(nf, 1), 4) else NA,
    bias = round(mean(acc_th[, m]) - tau0, 4), rmse = round(sqrt(mean((acc_th[, m] - tau0)^2)), 4),
    typeI = round(mean(acc_rn[, m]), 4), power = round(mean(acc_re[, m]), 4),
    coverage = round(mean(acc_cv[, m]), 4))))
}, error = function(e) NULL)

## ==== run ====================================================================
t0 <- Sys.time()
cat(sprintf("DIAG cells: %d  |  EST cells: %d  | nsim=%d B=%d cores=%d\n",
            nrow(diag_grid), nrow(est_cells), nsim, Bboot, ncore))
dR <- parallel::mclapply(seq_len(nrow(diag_grid)), diag_one, mc.cores = ncore)
cat(sprintf("  diag done (%.0f s)\n", as.numeric(difftime(Sys.time(), t0, "secs"))))
eR <- parallel::mclapply(seq_len(nrow(est_cells)), est_one, mc.cores = ncore)
cat(sprintf("  est done (%.0f s)\n", as.numeric(difftime(Sys.time(), t0, "secs"))))
out <- rbind(do.call(rbind, dR[!vapply(dR, is.null, TRUE)]),
             do.call(rbind, eR[!vapply(eR, is.null, TRUE)]))
write.csv(out, "run_all_results.csv", row.names = FALSE)
cat(sprintf("\nTOTAL rows: %d  (diag %d + est %d method-rows). NA metrics: %d. %.0f s\n",
            nrow(out), sum(out$block == "diag"), sum(out$block == "est"),
            sum(is.na(out$bias) & out$block == "est"), as.numeric(difftime(Sys.time(), t0, "secs"))))
cat("blocks x dgp:\n"); print(table(out$block, out$dgp))
