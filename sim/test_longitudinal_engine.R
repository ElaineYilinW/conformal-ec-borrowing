source("engine_longitudinal.R")

set.seed(20260717)
cfg <- long_default_config()
dat_alt <- long_generate_data(cfg)
dat <- long_set_null(dat_alt)

stopifnot(
  identical(dim(dat$Yt), c(cfg$n_r1, cfg$T)),
  identical(dim(dat$Yc), c(cfg$n_r0, cfg$T)),
  identical(dim(dat$Ye), c(cfg$n_ec, cfg$T)),
  all(abs((dat_alt$Yt - dat$Yt) - cfg$tau) < 1e-10),
  length(dat$B) == cfg$n_ec
)

fit <- long_fit_outcome(dat$Yc, dat$Xc, LONG_XVARS, longitudinal = TRUE)
stopifnot(
  !fit$fallback,
  is.finite(fit$sigma2), fit$sigma2 > 0,
  is.finite(fit$rho), fit$rho > -1 / (cfg$T - 1), fit$rho < 1,
  identical(dim(fit$predict_matrix(dat$Xc)), dim(dat$Yc))
)

# Every implemented symmetric score must be invariant to a sign flip.
r <- matrix(rnorm(20 * cfg$T), 20, cfg$T)
for (shape in c("avg", "final", "c1", "max", "quad")) {
  stopifnot(isTRUE(all.equal(
    long_score(r, long_cs_matrix(), shape),
    long_score(-r, long_cs_matrix(), shape),
    tolerance = 1e-12
  )))
}

spec <- long_working_spec("CC")
screen <- long_screen_cv(dat, spec, long_score_configs(TRUE), cfg$gamma, cfg$K)
stopifnot(
  all(vapply(screen$p, length, integer(1)) == cfg$n_ec),
  all(vapply(screen$p, function(x) all(is.finite(x) & x >= 0 & x <= 1), logical(1))),
  all(vapply(screen$selected, length, integer(1)) == cfg$n_ec)
)

# Density-ratio orientation: rho transports external covariates to the RCT,
# while omega transports RCT covariates to the external population.
set.seed(771)
xr_shift <- data.frame(X1 = rnorm(4000, 0.4, 1))
xe_shift <- data.frame(X1 = rnorm(4000, -0.4, 1))
src_shift <- long_fit_source(xr_shift, xe_shift, "X1")
wr <- src_shift$rho(xe_shift)
we <- src_shift$omega(xr_shift)
stopifnot(
  abs(weighted.mean(xe_shift$X1, wr) - mean(xr_shift$X1)) < 0.08,
  abs(weighted.mean(xr_shift$X1, we) - mean(xe_shift$X1)) < 0.08
)

est0 <- long_estimate_methods(dat, screen, spec, c("avg", "max", "quad"))$estimates
est1_direct <- long_estimate_methods(dat_alt, screen, spec, c("avg", "max", "quad"))$estimates
stopifnot(isTRUE(all.equal(unname(est1_direct - est0), rep(cfg$tau, length(est0)), tolerance = 1e-7)))

db <- long_resample_data(dat)
stopifnot(
  identical(dim(db$Yt), dim(dat$Yt)),
  identical(dim(db$Yc), dim(dat$Yc)),
  identical(dim(db$Ye), dim(dat$Ye))
)

# Large-sample GLS sanity check: the fitted covariance and shared coefficients
# should recover the DGP, not merely converge numerically.
set.seed(20260717)
cfg_big <- cfg
cfg_big$n_r1 <- 400L
cfg_big$n_r0 <- 400L
cfg_big$n_ec <- 100L
cfg_big$pi_bad <- 0
big <- long_set_null(long_generate_data(cfg_big, pool_mult = 10L))
fit_big <- long_fit_outcome(big$Yc, big$Xc, LONG_XVARS, longitudinal = TRUE)
cf <- fit_big$coefficients
stopifnot(
  abs(fit_big$sigma2 - cfg$sigma2) < 0.5,
  abs(fit_big$rho - cfg$rho) < 0.08,
  abs(cf["time"] - 0.05) < 0.03,
  abs(cf["Y0"] + 0.50) < 0.12,
  abs(cf["X1"] - cfg$outcome_x1) < 0.15
)

cat("All longitudinal engine checks passed.\n")
