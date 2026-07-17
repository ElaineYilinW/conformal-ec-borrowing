args <- commandArgs(trailingOnly = TRUE)
result_dir <- if (length(args) >= 1) args[1] else "results_longitudinal_main"
mode <- if (length(args) >= 2) args[2] else "main"

read_result <- function(suffix) {
  read.csv(file.path(result_dir, paste0("longitudinal_", mode, "_", suffix, ".csv")))
}

est <- read_result("estimation")
scr <- read_result("screening")
bor <- read_result("borrowing")
nui <- read_result("nuisance")
raw <- readRDS(file.path(result_dir, paste0("longitudinal_", mode, "_raw.rds")))

pick <- function(d, method) d[d$spec == "CC" & d$method == method, , drop = FALSE]
avg <- pick(est, "avg")
rct <- pick(est, "rct_long")
avg_scr <- scr[scr$spec == "CC" & scr$score == "avg", , drop = FALSE]
if (nrow(avg) != 1L || nrow(rct) != 1L || nrow(avg_scr) != 1L) {
  stop("Expected exactly one CC row for avg, rct_long, and avg screening.")
}

n <- avg$n_rep
mcse_type1 <- sqrt(0.05 * 0.95 / n)
mcse_power <- sqrt(avg$power * (1 - avg$power) / n)
mcse_bias <- avg$empirical_sd / sqrt(n)
mcse_coverage <- sqrt(0.95 * 0.05 / n)

valid_results <- raw$results[vapply(raw$results, is.list, logical(1))]
boot_valid <- do.call(rbind, lapply(valid_results, function(x) x$boot_valid))
min_boot_valid <- min(boot_valid, na.rm = TRUE)
boot_quality <- data.frame(
  method = colnames(boot_valid),
  mean_valid_fraction = colMeans(boot_valid, na.rm = TRUE),
  p05_valid_fraction = apply(boot_valid, 2, quantile, 0.05, na.rm = TRUE),
  min_valid_fraction = apply(boot_valid, 2, min, na.rm = TRUE),
  row.names = NULL
)

audit <- data.frame(
  criterion = c(
    "all_replications_valid", "bootstrap_valid_fraction", "bad_EC_detection",
    "clean_EC_wrong_detection", "absolute_bias", "RMSE_below_RCT_long",
    "coverage_compatible_95pct", "type1_compatible_5pct", "power_above_RCT_long"
  ),
  value = c(
    length(valid_results), min_boot_valid, avg_scr$detection,
    avg_scr$wrong_detection, abs(avg$bias), avg$rmse / rct$rmse,
    avg$coverage, avg$type1, avg$power - rct$power
  ),
  target = c(
    paste0(nrow(raw$jobs), "/", nrow(raw$jobs)), ">=0.99", ">=0.90", "<=0.10",
    "<=0.05 and within 2 MCSE", "<1", "within 2 MCSE of 0.95",
    "within 2 MCSE of 0.05", ">0"
  ),
  pass = c(
    length(valid_results) == nrow(raw$jobs), min_boot_valid >= 0.99,
    avg_scr$detection >= 0.90, avg_scr$wrong_detection <= 0.10,
    abs(avg$bias) <= 0.05 && abs(avg$bias) <= 2 * mcse_bias,
    avg$rmse < rct$rmse,
    abs(avg$coverage - 0.95) <= 2 * mcse_coverage,
    abs(avg$type1 - 0.05) <= 2 * mcse_type1,
    avg$power > rct$power
  )
)

out_path <- file.path(result_dir, paste0("longitudinal_", mode, "_gate_audit.csv"))
boot_path <- file.path(result_dir, paste0("longitudinal_", mode, "_bootstrap_quality.csv"))
write.csv(audit, out_path, row.names = FALSE)
write.csv(boot_quality, boot_path, row.names = FALSE)

cat(sprintf("CC avg audit: n=%d; elapsed=%.1f min\n", n, raw$elapsed / 60))
cat(sprintf(
  "bias=%+.4f (MCSE %.4f), RMSE ratio=%.3f, type-I=%.3f (MCSE %.3f), power=%.3f (MCSE %.3f)\n",
  avg$bias, mcse_bias, avg$rmse / rct$rmse, avg$type1, mcse_type1,
  avg$power, mcse_power
))
print(audit, row.names = FALSE)
cat(sprintf("Gate audit written to %s\n", out_path))
cat(sprintf("Method-level bootstrap quality written to %s\n", boot_path))
