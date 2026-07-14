# comprehensive simulation engine (v3) for the full blueprint grid----
# builds on engine.R (constants, score_vec, cov_stratum, fit_mean, wrank_p, wcv_p).
# adds: DGP 2 (homoscedastic) & DGP 3 (two covariates), non-symmetric two-sided
# screens, adaptive gamma (MSE), general transported AIPW + continuous density
# ratio (DGP 3), and a one-config metric driver. Base-R + closed-form fits only.

source("engine.R")

xcterm <- c(0.7, 1.0, 1.3, 1.6)                      # X_c effect over visits = 0.7+0.3t

## ---- DGP generators ---------------------------------------------------------
# per-subject sd given binary X and variance mode + source tag
sig_of <- function(X, vmode, src) {
  if (vmode == "hetero") ifelse(X == 1, sigma1, sigma0)
  else if (vmode == "source") rep(if (src == "E") 1.8 else 1.0, length(X))
  else rep(1.0, length(X))                            # homo
}
traj_from <- function(mu, sig) {
  n <- nrow(mu); b <- rnorm(n, 0, sd_b)
  eps <- matrix(rnorm(n * 4), n, 4) %*% R_chol
  mu + matrix(b, n, 4) + matrix(sig, n, 4) * eps
}
mu_bin <- function(X) matrix(mean_base, length(X), 4, byrow = TRUE) + outer(X, mean_xterm)
mu_two <- function(Xb, Xc) mu_bin(Xb) + outer(Xc, xcterm)

# generate a full trial. dgp in {1,2,3}. vmode: hetero (dgp1) / homo (dgp2) /
# source (dgp2 source-different). returns treated/control/ec pieces (+ Xc for dgp3).
gen_data <- function(dgp, n_r1, n_r0, n_ec, pattern, Delta, piB, effect,
                     vmode = "hetero") {
  if (dgp == 3) {
    Xbt <- rbinom(n_r1, 1, px_rct); Xct <- rnorm(n_r1, 0, 1)
    Xbc <- rbinom(n_r0, 1, px_rct); Xcc <- rnorm(n_r0, 0, 1)
    Xbe <- rbinom(n_ec, 1, px_ec);  Xce <- rnorm(n_ec, 0.6, 1)
    B <- rbinom(n_ec, 1, piB)
    tt <- if (effect) matrix(tau_vec, n_r1, 4, byrow = TRUE) else 0
    Yt <- traj_from(mu_two(Xbt, Xct), rep(1, n_r1)) + tt
    Yc <- traj_from(mu_two(Xbc, Xcc), rep(1, n_r0))
    M <- matrix(Delta * drift_vec(pattern), n_ec, 4, byrow = TRUE); if (pattern == "D") M <- M * Xbe
    Ye <- traj_from(mu_two(Xbe, Xce), rep(1, n_ec)) + M * B
    return(list(Yt = Yt, Xbt = Xbt, Xct = Xct, Yc = Yc, Xbc = Xbc, Xcc = Xcc,
                Ye = Ye, Xbe = Xbe, Xce = Xce, B = B))
  }
  # dgp 1/2 : single binary X ; vmode hetero/homo/source
  Xt <- rbinom(n_r1, 1, px_rct); Xc <- rbinom(n_r0, 1, px_rct); Xe <- rbinom(n_ec, 1, px_ec)
  B <- rbinom(n_ec, 1, piB)
  tt <- if (effect) matrix(tau_vec, n_r1, 4, byrow = TRUE) else 0
  Yt <- traj_from(mu_bin(Xt), sig_of(Xt, vmode, "R")) + tt
  Yc <- traj_from(mu_bin(Xc), sig_of(Xc, vmode, "R"))
  M <- matrix(Delta * drift_vec(pattern), n_ec, 4, byrow = TRUE); if (pattern == "D") M <- M * Xe
  Ye <- traj_from(mu_bin(Xe), sig_of(Xe, vmode, "E")) + M * B
  list(Yt = Yt, Xt = Xt, Yc = Yc, Xc = Xc, Ye = Ye, Xe = Xe, B = B)
}

## ---- DGP 3 helpers: two-covariate mean model, single cov, density ratio ------
# working mean model on (Xb,Xc): correct = per-visit lm(Y~Xb+Xc); moderate =
# common covariate slopes across visits; severe = omit Xc (per-visit lm(Y~Xb)).
fit_mean_two <- function(Y, Xb, Xc, level) {
  if (level == "moderate") {                          # common slopes, visit intercepts
    n <- nrow(Y); Yl <- as.vector(Y); vis <- factor(rep(1:4, each = n))
    Xbl <- rep(Xb, 4); Xcl <- rep(Xc, 4)
    fit <- lm(Yl ~ vis + Xbl + Xcl)
    co <- coef(fit); ai <- c(co[1], co[1] + co[2:4]); bb <- co["Xbl"]; cc <- co["Xcl"]
    return(function(nb, nc) outer(rep(1, length(nb)), ai) + outer(nb * bb + nc * cc, rep(1, 4)))
  }
  cols <- if (level == "severe") lapply(1:4, function(t) lm(Y[, t] ~ Xb))
          else lapply(1:4, function(t) lm(Y[, t] ~ Xb + Xc))
  function(nb, nc) {
    nd <- data.frame(Xb = nb, Xc = nc)
    sapply(1:4, function(t) as.numeric(predict(cols[[t]], nd)))
  }
}
# single 4x4 residual covariance (homoscedastic) with ridge
cov_one <- function(resid) cov(resid) + 1e-8 * diag(4)
# density ratio omega = f_E/f_R on (Xb,Xc) [screen weight]; rho = f_R/f_E [transport]
dr_dgp3 <- function(Xb, Xc, method, src = NULL, nR = NULL, nE = NULL) {
  if (method == "oracle") {
    fR <- dbinom(Xb, 1, px_rct) * dnorm(Xc, 0, 1)
    fE <- dbinom(Xb, 1, px_ec)  * dnorm(Xc, 0.6, 1)
    return(fE / fR)                                   # omega = f_E/f_R
  }
  # src = coefs of glm.fit(S ~ 1+Xb+Xc), S=1 RCT / S=0 EC, so eta = logit P(RCT|x).
  # omega = f_E/f_R = [P(EC|x)/P(RCT|x)]*(nR/nE) = ((1-p_rct)/p_rct)*(nR/nE) = exp(-eta)*(nR/nE).
  eta <- src[1] + src[2] * Xb + src[3] * Xc
  exp(-eta) * (nR / nE)                                # f_E/f_R via Bayes
}

## ---- weighted-rank p-value with explicit weights (continuous or binary) ------
wrank_pw <- function(v_cal, w_cal, v_test, w_test) {
  cmp <- outer(v_test, v_cal, function(a, b) b >= a)
  (w_test + as.numeric(cmp %*% w_cal)) / (w_test + sum(w_cal))
}
# single-covariance score (dgp3, no strata)
score_one <- function(resid, S, shape) {
  cvec <- if (shape %in% c("c1", "c1a")) C1 else if (shape %in% c("c2", "c2a")) C2 else NULL
  if (is.null(S)) S <- diag(4)
  switch(shape,
    quad = rowSums((resid %*% solve(S)) * resid),
    max  = apply(sweep(resid, 2, sqrt(diag(S)), "/"), 1, max),
    c1a = , c2a = abs(as.numeric(resid %*% cvec)) / sqrt(as.numeric(t(cvec) %*% S %*% cvec)),
    as.numeric(resid %*% cvec) / sqrt(as.numeric(t(cvec) %*% S %*% cvec)))
}

## ---- screen: EC tail p-values (hi/lo/abs) for one shape+cov+weight -----------
# dgp1/2 path uses binary strata (engine.R score_vec + cov_stratum + w_of).
# dgp3 path uses single cov + continuous density-ratio weights.
# base shape is the signed contrast ("c1"/"c2") or "quad"/"max"; hi/lo from signed.
screen_tails <- function(dgp, dat, shape, cov_type, weighted, level, prop = 0.5,
                         dratio = "oracle") {
  if (dgp == 3) {
    Yc <- dat$Yc; Xbc <- dat$Xbc; Xcc <- dat$Xcc
    Ye <- dat$Ye; Xbe <- dat$Xbe; Xce <- dat$Xce
    n <- nrow(Yc); tr <- sample.int(n, floor(prop * n)); cal <- setdiff(seq_len(n), tr)
    mf <- fit_mean_two(Yc[tr, , drop = FALSE], Xbc[tr], Xcc[tr], level)
    rc <- Yc[cal, , drop = FALSE] - mf(Xbc[cal], Xcc[cal]); re <- Ye - mf(Xbe, Xce)
    S <- if (cov_type == "raw") NULL else if (cov_type == "source") cov_one(re) else cov_one(Yc[tr, , drop = FALSE] - mf(Xbc[tr], Xcc[tr]))
    if (weighted) {
      src <- if (dratio == "logistic") fit_src_dgp3(Xbc, Xcc, Xbe, Xce) else NULL
      wc <- omega_dgp3(Xbc[cal], Xcc[cal], dratio, src, length(Xbc), length(Xbe))
      we <- omega_dgp3(Xbe, Xce, dratio, src, length(Xbc), length(Xbe))
    } else { wc <- rep(1, length(cal)); we <- rep(1, length(Xbe)) }
    sc <- score_one(rc, S, shape); se <- score_one(re, S, shape)
    absc <- if (shape %in% c("c1", "c2")) score_one(rc, S, paste0(shape, "a")) else abs(sc)
    abse <- if (shape %in% c("c1", "c2")) score_one(re, S, paste0(shape, "a")) else abs(se)
    return(list(hi = wrank_pw(sc, wc, se, we), lo = wrank_pw(-sc, wc, -se, we),
                abs = wrank_pw(absc, wc, abse, we)))
  }
  # dgp 1/2
  Yc <- dat$Yc; Xc <- dat$Xc; Ye <- dat$Ye; Xe <- dat$Xe
  n <- nrow(Yc); tr <- sample.int(n, floor(prop * n)); cal <- setdiff(seq_len(n), tr)
  mf <- fit_mean(Yc[tr, , drop = FALSE], Xc[tr], level)
  covR <- if (cov_type == "raw") NULL else cov_stratum(Yc[tr, , drop = FALSE] - mf(Xc[tr]), Xc[tr])
  covE <- if (cov_type == "source") cov_stratum(Ye - mf(Xe), Xe) else covR
  rc <- Yc[cal, , drop = FALSE] - mf(Xc[cal]); Xcal <- Xc[cal]; re <- Ye - mf(Xe)
  wc <- w_of(Xcal, weighted); we <- w_of(Xe, weighted)
  sc <- score_vec(rc, Xcal, covR, shape); se <- score_vec(re, Xe, covE, shape)
  absh <- if (shape %in% c("c1", "c2")) paste0(shape, "a") else shape
  absc <- if (shape %in% c("c1", "c2")) score_vec(rc, Xcal, covR, absh) else abs(sc)
  abse <- if (shape %in% c("c1", "c2")) score_vec(re, Xe, covE, absh) else abs(se)
  list(hi = wrank_pw(sc, wc, se, we), lo = wrank_pw(-sc, wc, -se, we),
       abs = wrank_pw(absc, wc, abse, we))
}

## ---- selection given orientation + threshold(s) -----------------------------
sel_screen <- function(p, orient, gh, gl = 0) {
  if (orient == "sym") p$abs > gh
  else if (orient == "one") p$hi > gh
  else (p$hi > gh) & (p$lo > gl)                       # non-symmetric two-sided
}

## ---- estimators -------------------------------------------------------------
# dgp1/2: per-stratum transported estimator (engine.R est_tau_c). dgp3: general AIPW.
est_tau_dgp3 <- function(Zt, Zr, Xbc, Xcc, Ze, Xbe, Xce, sel, dratio = "logistic") {
  mu1 <- mean(Zt); mu0r <- mean(Zr); vR <- var(Zr) / length(Zr)
  ret <- which(sel)
  if (length(ret) < 3) return(mu1 - mu0r)
  Xbr <- Xbe[ret]; Xcr <- Xce[ret]
  df <- data.frame(Z = Ze[ret], Xb = Xbr, Xc = Xcr)
  fit <- lm(Z ~ Xb + Xc, data = df)
  pr <- predict(fit, data.frame(Xb = Xbc, Xc = Xcc)); gform <- mean(pr)
  # transport ratio rho_gamma = f_R / f_{E,gamma}, RE-ESTIMATED ON THE RETAINED sample
  # (blueprint eq 15: f_{E,gamma} is the covariate law among retained ECs). Refit the
  # source model S~X on RCT controls (S=1) vs RETAINED ECs (S=0), not the full EC set.
  if (dratio == "oracle") {
    w <- 1 / dr_dgp3(Xbr, Xcr, "oracle"); w <- pmin(pmax(w, 1e-3), quantile(w, 0.99, na.rm = TRUE))
  } else {
    src <- fit_src_dgp3(Xbc, Xcc, Xbr, Xcr)
    w <- pmin(1 / omega_dgp3(Xbr, Xcr, "logistic", src, length(Xbc), length(ret)), 50)
  }
  resid_e <- df$Z - predict(fit)
  aug <- sum(w * resid_e) / sum(w); mu0e <- gform + aug
  vE <- var(pr) / length(Zr) + sum((w * resid_e)^2) / (sum(w)^2)
  lam <- if (!is.finite(vE) || vE <= 0) 0 else vR / (vR + vE)
  mu1 - ((1 - lam) * mu0r + lam * mu0e)
}
# unified tau: given data, contrast c, selection, dgp. dratio picks the DGP3 ratio.
tau_of <- function(dgp, dat, cc, sel, dratio = "logistic") {
  Zt <- as.numeric(dat$Yt %*% cc); Zr <- as.numeric(dat$Yc %*% cc); Ze <- as.numeric(dat$Ye %*% cc)
  if (dgp == 3) est_tau_dgp3(Zt, Zr, dat$Xbc, dat$Xcc, Ze, dat$Xbe, dat$Xce, sel, dratio)
  else est_tau_c(Zt, Zr, dat$Xc, Ze, dat$Xe, sel)
}

## ---- logistic density ratio for DGP 3 (estimated weights) --------------------
# fit source model S~Xb+Xc (S=1 RCT control, S=0 EC) and return omega=f_E/f_R fn
fit_src_dgp3 <- function(Xbc, Xcc, Xbe, Xce) {
  X <- cbind(1, c(Xbc, Xbe), c(Xcc, Xce))               # design: intercept, Xb, Xc
  y <- c(rep(1, length(Xbc)), rep(0, length(Xbe)))      # S=1 RCT control, S=0 EC
  suppressWarnings(glm.fit(X, y, family = binomial()))$coefficients   # fast; returns (b0,bXb,bXc)
}
omega_dgp3 <- function(Xb, Xc, dratio, src = NULL, nR = NULL, nE = NULL) {
  w <- if (dratio == "oracle") dr_dgp3(Xb, Xc, "oracle")
       else dr_dgp3(Xb, Xc, "logistic", src = src, nR = nR, nE = nE)
  cap <- quantile(w, 0.99, na.rm = TRUE); pmin(pmax(w, 1e-3), cap)
}

## ---- screening-only driver: FER + detection for one config (no bootstrap) ----
# returns per-replicate deleted-clean / deleted-contaminated counts.
screen_metrics <- function(dgp, vmode, shape, cov_type, weighted, level, calib,
                           pattern, Delta, piB, n_r0, n_ec) {
  dat <- gen_data(dgp, 2L, n_r0, n_ec, pattern, Delta, piB, FALSE, vmode = vmode)
  P <- screen_tails(dgp, dat, shape, cov_type, weighted, level)
  del <- if (shape %in% c("quad", "max")) (P$hi <= 0.10) else (P$abs <= 0.10)  # sym/one at gamma .10
  active <- if (pattern == "D" && dgp != 3) (dat$B == 1 & dat$Xe == 1)
            else if (pattern == "D" && dgp == 3) (dat$B == 1 & dat$Xbe == 1)
            else (dat$B == 1)
  clean <- !active
  c(del_clean = sum(del & clean), n_clean = sum(clean),
    del_active = sum(del & active), n_active = sum(active))
}

## ---- CV+ calibration tails (dgp1/2) for the diagnostic block -----------------
screen_tails_cv <- function(dgp, dat, shape, cov_type, weighted, level, K = 5) {
  Yc <- dat$Yc; Xc <- dat$Xc; Ye <- dat$Ye; Xe <- dat$Xe
  n <- nrow(Yc); fold <- sample(rep(seq_len(K), length.out = n))
  mfull <- fit_mean(Yc, Xc, level)
  covR <- if (cov_type == "raw") NULL else cov_stratum(Yc - mfull(Xc), Xc)
  covE <- if (cov_type == "source") cov_stratum(Ye - mfull(Xe), Xe) else covR
  mods <- lapply(seq_len(K), function(k) fit_mean(Yc[fold != k, , drop = FALSE], Xc[fold != k], level))
  residC <- matrix(0, n, 4)
  for (k in seq_len(K)) { idx <- which(fold == k); residC[idx, ] <- Yc[idx, , drop = FALSE] - mods[[k]](Xc[idx]) }
  absh <- if (shape %in% c("c1", "c2")) paste0(shape, "a") else shape
  scC <- score_vec(residC, Xc, covR, shape)
  scCa <- if (shape %in% c("c1", "c2")) score_vec(residC, Xc, covR, absh) else abs(scC)
  vhi <- vlo <- vab <- matrix(0, nrow(Ye), K)
  for (k in seq_len(K)) {
    re <- Ye - mods[[k]](Xe); s <- score_vec(re, Xe, covE, shape)
    vhi[, k] <- s; vlo[, k] <- -s
    vab[, k] <- if (shape %in% c("c1", "c2")) score_vec(re, Xe, covE, absh) else abs(s)
  }
  list(hi = wcv_p(scC, fold, vhi, Xc, Xe, weighted),
       lo = wcv_p(-scC, fold, vlo, Xc, Xe, weighted),
       abs = wcv_p(scCa, fold, vab, Xc, Xe, weighted))
}
# dispatch split/cv
tails <- function(dgp, dat, shape, cov_type, weighted, level, calib = "split", dratio = "oracle") {
  if (calib == "cv" && dgp != 3) screen_tails_cv(dgp, dat, shape, cov_type, weighted, level)
  else screen_tails(dgp, dat, shape, cov_type, weighted, level, dratio = dratio)
}
