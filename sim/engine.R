# consolidated simulation engine (v2)----
# toy: 4 visits, binary severity X, RCT vs EC covariate shift.
# two estimands, four drift patterns x Delta, piB=0.5. twelve nonconformity
# scores: shape {quad, max, c1, c2} x cov {raw, internal, source};
# quad unsigned, max and contrast signed. weighted/unweighted x split/CV.

# constants (blueprint eq 17-22)----
R_mat <- 0.5^abs(outer(0:3, 0:3, "-"))
R_chol <- chol(R_mat)
sigma0 <- 1
sigma1 <- 1.8
sd_b <- 0.5
px_rct <- 0.30
px_ec <- 0.70
mean_base <- c(2, 2.5, 3, 3.5)
mean_xterm <- c(1, 1.4, 1.8, 2.2)
tau_vec <- c(0, 0.25, 0.60, 1.00)
C1 <- c(-1, 0, 0, 1)
C2 <- c(1, 1, 1, 1) / 4
estimands <- list(c1 = C1, c2 = C2)
tau_true <- c(c1 = sum(C1 * tau_vec), c2 = sum(C2 * tau_vec))
omega0 <- c(mild = (1 - px_ec) / (1 - px_rct), severe = px_ec / px_rct)

# data-generating process----
gen_traj <- function(X) {
  n <- length(X)
  b <- rnorm(n, 0, sd_b)
  eps <- matrix(rnorm(n * 4), n, 4) %*% R_chol
  sig <- ifelse(X == 1, sigma1, sigma0)
  mu <- matrix(mean_base, n, 4, byrow = TRUE) + outer(X, mean_xterm)
  mu + matrix(b, n, 4) + matrix(sig, n, 4) * eps
}

gen_controls <- function(n, px) {
  X <- rbinom(n, 1, px)
  list(Y = gen_traj(X), X = X)
}

gen_treated <- function(n, px) {
  X <- rbinom(n, 1, px)
  list(Y = gen_traj(X) + matrix(tau_vec, n, 4, byrow = TRUE), X = X)
}

drift_vec <- function(pattern) {
  switch(pattern, A = c(1, 1, 1, 1), B = c(0, 0, 0, 1),
         C = c(0, 1/3, 2/3, 1), D = c(0, 1/3, 2/3, 1))
}

gen_ec <- function(n, pattern, Delta, piB) {
  X <- rbinom(n, 1, px_ec)
  B <- rbinom(n, 1, piB)
  M <- matrix(Delta * drift_vec(pattern), n, 4, byrow = TRUE)
  if (pattern == "D") M <- M * X
  list(Y = gen_traj(X) + M * B, X = X, B = B)
}

# working mean model: correct / moderate / severe----
fit_mean <- function(Y, X, level) {
  ybar <- colMeans(Y)
  if (level == "severe") {
    return(function(nx) matrix(ybar, length(nx), 4, byrow = TRUE))
  }
  if (level == "moderate") {
    xbar <- mean(X)
    xc <- X - xbar
    den <- 4 * sum(xc^2)
    gam <- if (den > 0) sum(vapply(1:4, function(t) sum(xc * Y[, t]), numeric(1))) / den else 0
    return(function(nx) {
      matrix(ybar, length(nx), 4, byrow = TRUE) + outer((nx - xbar) * gam, rep(1, 4))
    })
  }
  m0 <- colMeans(Y[X == 0, , drop = FALSE])
  m1 <- colMeans(Y[X == 1, , drop = FALSE])
  function(nx) {
    out <- matrix(0, length(nx), 4)
    if (any(nx == 0)) out[nx == 0, ] <- matrix(m0, sum(nx == 0), 4, byrow = TRUE)
    if (any(nx == 1)) out[nx == 1, ] <- matrix(m1, sum(nx == 1), 4, byrow = TRUE)
    out
  }
}

# per-stratum residual covariance (used for internal and source scaling)----
cov_stratum <- function(resid, X) {
  ridge <- 1e-8 * diag(4)
  pooled <- cov(resid) + ridge
  c0 <- if (sum(X == 0) >= 4) cov(resid[X == 0, , drop = FALSE]) + ridge else pooled
  c1 <- if (sum(X == 1) >= 4) cov(resid[X == 1, , drop = FALSE]) + ridge else pooled
  list(`0` = c0, `1` = c1)
}

# scalar nonconformity score from residual matrix----
# covlist = NULL for raw (identity); otherwise per-stratum covariance.
score_vec <- function(resid, X, covlist, shape) {
  v <- numeric(nrow(resid))
  cvec <- if (shape %in% c("c1", "c1a")) C1 else if (shape %in% c("c2", "c2a")) C2 else NULL
  for (x in c(0, 1)) {
    idx <- which(X == x)
    if (!length(idx)) next
    r <- resid[idx, , drop = FALSE]
    S <- if (is.null(covlist)) diag(4) else covlist[[as.character(x)]]
    v[idx] <- switch(shape,
      quad = rowSums((r %*% solve(S)) * r),
      max = apply(sweep(r, 2, sqrt(diag(S)), "/"), 1, max),
      c1a = , c2a = abs(as.numeric(r %*% cvec)) / sqrt(as.numeric(t(cvec) %*% S %*% cvec)),
      as.numeric(r %*% cvec) / sqrt(as.numeric(t(cvec) %*% S %*% cvec)))
  }
  as.numeric(v)
}

# weights----
w_of <- function(X, weighted) if (weighted) ifelse(X == 1, omega0["severe"], omega0["mild"]) else rep(1, length(X))

# weighted upper-tail conformal p-value from calibration and test scores----
wrank_p <- function(v_cal, X_cal, v_test, X_test, weighted) {
  wc <- w_of(X_cal, weighted)
  wt <- w_of(X_test, weighted)
  cmp <- outer(v_test, v_cal, FUN = function(a, b) b >= a)
  (wt + as.numeric(cmp %*% wc)) / (wt + sum(wc))
}

# weighted cv+ p-value, eq (13)----
wcv_p <- function(v_cal, fold, v_ec_by, Xc, Xe, weighted) {
  wc <- w_of(Xc, weighted)
  we <- w_of(Xe, weighted)
  num <- we
  for (k in sort(unique(fold))) {
    ik <- which(fold == k)
    cmp <- outer(v_ec_by[, k], v_cal[ik], FUN = function(a, b) b >= a)
    num <- num + as.numeric(cmp %*% wc[ik])
  }
  num / (we + sum(wc))
}

# three covariance treatments for scoring (cal cov, ec cov)----
cov_choices <- function(covR, covE) {
  list(raw = list(NULL, NULL), internal = list(covR, covR), source = list(covR, covE))
}

# cv+ screen: returns EC p-value vector per shape.cov.weight key----
screen_cv <- function(Yc, Xc, Ye, Xe, level, K = 5) {
  n <- nrow(Yc)
  fold <- sample(rep(seq_len(K), length.out = n))
  mfull <- fit_mean(Yc, Xc, level)
  covsets <- cov_choices(cov_stratum(Yc - mfull(Xc), Xc), cov_stratum(Ye - mfull(Xe), Xe))
  mods <- lapply(seq_len(K), function(k) fit_mean(Yc[fold != k, , drop = FALSE], Xc[fold != k], level))
  residC <- matrix(0, n, 4)
  for (k in seq_len(K)) {
    idx <- which(fold == k)
    residC[idx, ] <- Yc[idx, , drop = FALSE] - mods[[k]](Xc[idx])
  }
  residE <- lapply(seq_len(K), function(k) Ye - mods[[k]](Xe))
  out <- list()
  for (shape in c("quad", "max", "c1", "c2")) for (cn in names(covsets)) {
    cc <- covsets[[cn]]
    vcal <- score_vec(residC, Xc, cc[[1]], shape)
    vec <- vapply(seq_len(K), function(k) score_vec(residE[[k]], Xe, cc[[2]], shape), numeric(nrow(Ye)))
    for (wt in c(FALSE, TRUE)) {
      out[[paste(shape, cn, if (wt) "w" else "u", sep = ".")]] <- wcv_p(vcal, fold, vec, Xc, Xe, wt)
    }
  }
  out
}

# one-time split screen: same keys----
screen_split <- function(Yc, Xc, Ye, Xe, level, prop = 0.5) {
  n <- nrow(Yc)
  tr <- sample.int(n, floor(prop * n))
  cal <- setdiff(seq_len(n), tr)
  mf <- fit_mean(Yc[tr, , drop = FALSE], Xc[tr], level)
  covsets <- cov_choices(cov_stratum(Yc[tr, , drop = FALSE] - mf(Xc[tr]), Xc[tr]),
                         cov_stratum(Ye - mf(Xe), Xe))
  residC <- Yc[cal, , drop = FALSE] - mf(Xc[cal])
  Xcal <- Xc[cal]
  residE <- Ye - mf(Xe)
  out <- list()
  for (shape in c("quad", "max", "c1", "c2")) for (cn in names(covsets)) {
    cc <- covsets[[cn]]
    vcal <- score_vec(residC, Xcal, cc[[1]], shape)
    vec <- score_vec(residE, Xe, cc[[2]], shape)
    for (wt in c(FALSE, TRUE)) {
      out[[paste(shape, cn, if (wt) "w" else "u", sep = ".")]] <- wrank_p(vcal, Xcal, vec, Xe, wt)
    }
  }
  out
}

# lean cv+ screen for a single (shape, cov, weight), for the bootstrap loop----
screen_one <- function(Yc, Xc, Ye, Xe, level, shape, cov_type, weighted, K = 5) {
  n <- nrow(Yc)
  fold <- sample(rep(seq_len(K), length.out = n))
  mfull <- fit_mean(Yc, Xc, level)
  covR <- if (cov_type == "raw") NULL else cov_stratum(Yc - mfull(Xc), Xc)
  covE <- if (cov_type == "source") cov_stratum(Ye - mfull(Xe), Xe) else covR
  mods <- lapply(seq_len(K), function(k) fit_mean(Yc[fold != k, , drop = FALSE], Xc[fold != k], level))
  residC <- matrix(0, n, 4)
  for (k in seq_len(K)) {
    idx <- which(fold == k)
    residC[idx, ] <- Yc[idx, , drop = FALSE] - mods[[k]](Xc[idx])
  }
  vcal <- score_vec(residC, Xc, covR, shape)
  vec <- vapply(seq_len(K), function(k) score_vec(Ye - mods[[k]](Xe), Xe, covE, shape), numeric(nrow(Ye)))
  wcv_p(vcal, fold, vec, Xc, Xe, weighted)
}

# transported contrast aipw for one estimand, inverse-variance borrowing----
est_tau_c <- function(Zt, Zr, Xc, Ze, Xe, sel) {
  mu1 <- mean(Zt)
  mu0r <- mean(Zr)
  fr <- c(mean(Xc == 0), mean(Xc == 1))
  zrx <- function(x) if (any(Xc == x)) mean(Zr[Xc == x]) else mu0r
  zr <- c(zrx(0), zrx(1))
  vR <- var(Zr) / length(Zr)
  m0e <- numeric(2)
  vt <- numeric(2)
  for (x in 0:1) {
    zx <- Ze[sel & Xe == x]
    if (length(zx) >= 2) {
      m0e[x + 1] <- mean(zx)
      vt[x + 1] <- fr[x + 1]^2 * var(zx) / length(zx)
    } else {
      m0e[x + 1] <- zr[x + 1]
      vt[x + 1] <- fr[x + 1]^2 * vR
    }
  }
  vE <- sum(vt)
  lam <- if (sum(sel) == 0 || !is.finite(vE)) 0 else vR / (vR + vE)
  mu1 - ((1 - lam) * mu0r + lam * sum(fr * m0e))
}
