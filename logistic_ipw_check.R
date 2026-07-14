# is my empirical-proportion IPW == literal glm(S~X) logistic-source IPW?----
# binary X (DGP1). for the SAME non-symmetric kept-set, compute mu0 THREE ways:
#   gform    : g-formula  sum_x frR(x)*mean(Z|kept,X=x)             [my est_tau_c]
#   ipw_emp  : Hajek IPW, rho = frR/f_{E,kept} from EMPIRICAL props [reweight_check]
#   ipw_glm  : Hajek IPW, rho from LITERAL glm(S ~ factor(X)),      [TEXTBOOK IPW]
#              S=1 RCT control, S=0 kept EC; w = p/(1-p), p=P(S=1|X)
# claim: for binary X, glm(S~X) is SATURATED -> ipw_glm == ipw_emp == gform to
# machine precision. print max|diff| and the shared bias. (DGP1 A Delta=5 c2.)

source("engine2.R")
args <- commandArgs(trailingOnly = TRUE)
nsim  <- if (length(args) >= 1) as.integer(args[1]) else 500L
ncore <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores() - 2L)
set.seed(20260713L)
n_r1 <- 120L; n_r0 <- 60L; n_ec <- 200L; piB <- 0.30; Delta <- 5; pattern <- "A"
cc <- C2; tau0 <- sum(cc * tau_vec); ghi <- 0.20; glo <- 0.05

mu0_three <- function(Zr, Xc, Ze, Xe, sel) {
  frR <- c(mean(Xc == 0), mean(Xc == 1))
  k <- which(sel); Zk <- Ze[k]; Xk <- Xe[k]
  mret <- c(if (any(Xk == 0)) mean(Zk[Xk == 0]) else mean(Zr),
            if (any(Xk == 1)) mean(Zk[Xk == 1]) else mean(Zr))
  gform <- sum(frR * mret)
  # empirical-proportion IPW (Hajek)
  fEk <- c(mean(Xk == 0), mean(Xk == 1)); fEk[fEk < 1e-9] <- 1e-9
  rho_emp <- (frR / fEk)[Xk + 1]
  ipw_emp <- sum(rho_emp * Zk) / sum(rho_emp)
  # LITERAL glm(S ~ factor(X)) logistic source model, S=1 RCT, S=0 kept EC
  Sdf <- data.frame(S = c(rep(1, length(Zr)), rep(0, length(Zk))),
                    X = factor(c(Xc, Xk), levels = c(0, 1)))
  fit <- suppressWarnings(glm(S ~ X, data = Sdf, family = binomial()))
  p <- predict(fit, data.frame(X = factor(Xk, levels = c(0, 1))), type = "response")
  p <- pmin(pmax(p, 1e-9), 1 - 1e-9)
  rho_glm <- p / (1 - p)                     # w(x) proportional to f_R/f_E
  ipw_glm <- sum(rho_glm * Zk) / sum(rho_glm)
  c(gform = gform, ipw_emp = ipw_emp, ipw_glm = ipw_glm)
}
one <- function(r) tryCatch({
  set.seed(20260713L + r)
  dat <- gen_data(1L, n_r1, n_r0, n_ec, pattern, Delta, piB, TRUE)
  P <- screen_tails(1L, dat, "c2", "internal", TRUE, "correct")
  Zt <- as.numeric(dat$Yt %*% cc); Zr <- as.numeric(dat$Yc %*% cc); Ze <- as.numeric(dat$Ye %*% cc)
  sel <- (P$hi > ghi) & (P$lo > glo)
  mean(Zt) - mu0_three(Zr, dat$Xc, Ze, dat$Xe, sel)
}, error = function(e) NULL)

res <- parallel::mclapply(seq_len(nsim), one, mc.cores = ncore)
res <- res[vapply(res, function(x) is.numeric(x) && length(x) == 3, TRUE)]
M <- do.call(rbind, res)
maxdiff <- max(abs(M[, "ipw_glm"] - M[, "ipw_emp"]), abs(M[, "ipw_emp"] - M[, "gform"]))
cat(sprintf("\n== binary X: g-formula vs empirical-IPW vs LITERAL glm(S~X) IPW (DGP1 A Delta=5 c2; %d reps) ==\n", nrow(M)))
cat(sprintf("per-rep max |ipw_glm - ipw_emp - gform| = %.2e   (0 => algebraically identical)\n\n", maxdiff))
print(data.frame(estimator = c("g-formula (my est_tau_c)", "IPW empirical props", "IPW literal glm(S~X)"),
                 bias = round(colMeans(M) - tau0, 4), rmse = round(sqrt(colMeans((M - tau0)^2)), 4)),
      row.names = FALSE)
