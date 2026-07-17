# sunfish / jrss-a data-generating engine----
# follows zhou, zhu, drake & pang (2025) jrss-a 188:791-818, section 5, and the
# rdborrow vignette "primary_simulation_workflow.Rmd" for the numerical values.
#
# covariates      : (x1..x4) via a gaussian ar1 copula; x5 continuous, from x1..x4.
# selection       : logit pi_S(W) = alpha'W, sampled to fixed n (paper's formulation).
#                   settings 3,4 add an x4^2 term to the TRUE selection model.
# outcome         : Y_t = A*effect_t + beta_t'[1,x1..x5] + b_i + s_t*eps_it.
#                   settings 2,4 add an x4^2 term to the TRUE outcome model.
# unmeasured conf : setting 5 hides x5 (it becomes U); only x1..x4 are observable.
#
# deviation from jrss-a, on purpose: their noise is drawn independently at each
# visit, so the residual covariance is diagonal and every longitudinal score
# (mahalanobis / contrast) collapses to a per-visit one. we split each visit's
# variance into a subject random intercept plus an ar1 within-subject term,
# preserving their marginal sd.

# constants----
rho_cop <- 0.8
rho_ar <- 0.6
icc <- 0.5
sd_visit <- c(4, 4, 4, 5)
sd_x5 <- 10

# x4 is truncated: the untruncated exp(1/10) tail reaches ~110, at which point a
# quadratic term swamps the sd-4 noise by two orders of magnitude. the rdborrow
# synthetic data tops out at 32, and sunfish part 2 enrolled ages 2-25.
x4_cap <- 32

# jrss-a outcome coefficients, ordered (intercept, x1, x2, x3, x4, x5)
beta_visit <- list(
  c(10.0, 0.05, -1.5, -1.0, -0.2, -0.10),
  c(6.0, 0.50, -0.5, -1.0, -0.3, -0.06),
  c(5.0, 1.90, 1.4, -1.3, -0.4, -0.15),
  c(1.2, 1.00, 2.0, -0.5, -0.4, -0.10)
)

# x5 conditional mean in the source population (rdborrow vignette, internal arm)
x5_coef <- c(30, 10, 7, -6, -0.5)

# treatment effect per visit; grows to the vignette's alt_effect = 2.0
tau_sunfish <- c(0, 0.8, 1.4, 2.0)

# selection model coefficients (intercept, x1, x2, x3, x4, x5); alpha_q is the
# x4^2 term used only by the TRUE model in settings 3 and 4. alpha5 drives the
# x5 shift, reproducing the vignette's ~20-point gap in E[x5] between sources.
alpha_sel <- c(2.2, 0.30, 0.40, -0.30, -0.02, -0.10)
alpha_q <- 0.003
beta_q <- c(0.010, 0.010, 0.015, 0.015)

# contrasts, reused from engine.R conventions
C1_sf <- c(-1, 0, 0, 1)
C2_sf <- c(1, 1, 1, 1) / 4

# which settings carry a quadratic term in which true model----
q_outcome <- function(setting) setting %in% c(2, 4)
q_select <- function(setting) setting %in% c(3, 4)
# setting 5 hides x5: it is the unmeasured confounder U
obs_vars <- function(setting) if (setting == 5) c("x1", "x2", "x3", "x4") else c("x1", "x2", "x3", "x4", "x5")

expit <- function(z) 1 / (1 + exp(-z))

# covariates: gaussian ar1 copula over binom/binom/binom/exp margins----
gen_X_copula <- function(n) {
  R <- rho_cop^abs(outer(1:4, 1:4, "-"))
  U <- pnorm(matrix(rnorm(n * 4), n, 4) %*% chol(R))
  x1 <- qbinom(U[, 1], 1, 0.7)
  x2 <- qbinom(U[, 2], 1, 0.9)
  x3 <- qbinom(U[, 3], 1, 0.3)
  x4 <- round(qexp(U[, 4] * pexp(x4_cap - 1, rate = 1 / 10), rate = 1 / 10)) + 1
  x5 <- x5_coef[1] + x5_coef[2] * x1 + x5_coef[3] * x2 + x5_coef[4] * x3 +
    x5_coef[5] * x4 + rnorm(n, 0, sd_x5)
  data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)
}

# correlated visit noise: random intercept + ar1, preserving the marginal sd----
gen_noise <- function(n) {
  R <- rho_ar^abs(outer(1:4, 1:4, "-"))
  b <- rnorm(n, 0, sqrt(icc) * sd_visit[1])
  eps <- matrix(rnorm(n * 4), n, 4) %*% chol(R)
  s_t <- sqrt(1 - icc) * sd_visit
  matrix(b, n, 4) + sweep(eps, 2, s_t, "*")
}

# true conditional mean of the four visits given covariates and treatment----
mu_sunfish <- function(X, A, setting, effect) {
  n <- nrow(X)
  D <- cbind(1, X$x1, X$x2, X$x3, X$x4, X$x5)
  out <- vapply(1:4, function(t) {
    m <- as.numeric(D %*% beta_visit[[t]]) + A * effect[t]
    if (q_outcome(setting)) m <- m + beta_q[t] * X$x4^2
    m
  }, numeric(n))
  matrix(out, n, 4)
}

# drift applied to contaminated ecs (settings 1-4 only)----
drift_sunfish <- function(pattern) {
  switch(pattern,
    A = c(1, 1, 1, 1), B = c(0, 0, 0, 1),
    C = c(0, 1 / 3, 2 / 3, 1), D = c(0, 1 / 3, 2 / 3, 1)
  )
}

# generate one trial: sample rct / ec from a pool by the true selection model----
# effect is the per-visit treatment effect vector (zeros under the null).
# piB / Delta / pattern add binary contamination; setting 5 ignores them and
# relies on the unmeasured confounder x5 for incompatibility instead.
gen_data_sunfish <- function(setting, n_r1 = 120, n_r0 = 60, n_ec = 200,
                             effect = rep(0, 4), piB = 0.30, Delta = 3,
                             pattern = "A", pool_mult = 40) {
  n_r <- n_r1 + n_r0
  pool <- gen_X_copula(pool_mult * (n_r + n_ec))
  lp <- alpha_sel[1] + alpha_sel[2] * pool$x1 + alpha_sel[3] * pool$x2 +
    alpha_sel[4] * pool$x3 + alpha_sel[5] * pool$x4 + alpha_sel[6] * pool$x5
  if (q_select(setting)) lp <- lp + alpha_q * pool$x4^2
  pi_S <- expit(lp)

  idx_r <- sample.int(nrow(pool), n_r, prob = pi_S)
  idx_e <- sample.int(nrow(pool), n_ec, prob = 1 - pi_S)
  Xr <- pool[idx_r, ]
  Xe <- pool[idx_e, ]

  A <- rep(0L, n_r)
  A[sample.int(n_r, n_r1)] <- 1L

  Yr <- mu_sunfish(Xr, A, setting, effect) + gen_noise(n_r)
  Ye <- mu_sunfish(Xe, rep(0L, n_ec), setting, effect) + gen_noise(n_ec)

  B <- rep(0L, n_ec)
  if (setting != 5 && piB > 0) {
    B <- rbinom(n_ec, 1, piB)
    M <- matrix(Delta * drift_sunfish(pattern), n_ec, 4, byrow = TRUE)
    if (pattern == "D") M <- M * Xe$x1
    Ye <- Ye + M * B
  }

  list(
    Yt = Yr[A == 1, , drop = FALSE], Xt = Xr[A == 1, ],
    Yc = Yr[A == 0, , drop = FALSE], Xc = Xr[A == 0, ],
    Ye = Ye, Xe = Xe, B = B, setting = setting
  )
}

# working outcome model: per-visit lm on the OBSERVED covariates----
# returns a predictor function; coefficients are genuinely estimated here.
fit_mean_sunfish <- function(Y, X, setting) {
  vars <- obs_vars(setting)
  fm <- reformulate(vars, response = "y")
  fits <- lapply(1:4, function(t) lm(fm, data = cbind(X[vars], y = Y[, t])))
  function(newX) {
    vapply(1:4, function(t) as.numeric(predict(fits[[t]], newX)), numeric(nrow(newX)))
  }
}

# return the fitted coefficient matrix (4 visits x p), for inspection----
coef_sunfish <- function(Y, X, setting) {
  vars <- obs_vars(setting)
  fm <- reformulate(vars, response = "y")
  t(vapply(1:4, function(t) coef(lm(fm, data = cbind(X[vars], y = Y[, t]))),
    numeric(length(vars) + 1)
  ))
}

# working selection model: logistic on the OBSERVED covariates----
# returns omega = f_E/f_R evaluated at newX, via bayes on the fitted propensity.
fit_omega_sunfish <- function(Xc, Xe, setting) {
  vars <- obs_vars(setting)
  D <- as.matrix(cbind(1, rbind(Xc[vars], Xe[vars])))
  y <- c(rep(1, nrow(Xc)), rep(0, nrow(Xe)))
  cf <- suppressWarnings(glm.fit(D, y, family = binomial()))$coefficients
  nR <- nrow(Xc)
  nE <- nrow(Xe)
  function(newX) {
    eta <- as.numeric(as.matrix(cbind(1, newX[vars])) %*% cf)
    w <- exp(-eta) * (nR / nE)
    pmin(pmax(w, 1e-3), quantile(w, 0.99, na.rm = TRUE))
  }
}
