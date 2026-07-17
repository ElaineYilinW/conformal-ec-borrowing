source("engine_longitudinal.R")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[1] else "main"
nsim <- if (length(args) >= 2) as.integer(args[2]) else if (mode == "main") 200L else 100L
Bboot <- if (length(args) >= 3) as.integer(args[3]) else if (mode == "main") 500L else 100L
detected_cores <- parallel::detectCores()
if (!is.finite(detected_cores)) detected_cores <- 4L
ncore <- if (length(args) >= 4) as.integer(args[4]) else max(1L, detected_cores - 2L)
out_dir <- if (length(args) >= 5) args[5] else "."
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cfg <- long_default_config()
if (mode == "misspec") {
  cfg$source_strength <- "strong_x1"
  cfg$outcome_x1 <- 1.20
}
if (mode == "clean") cfg$pi_bad <- 0

spec_names <- if (mode == "misspec") c("CC", "MC", "CM", "MM") else "CC"
score_configs <- if (mode == "main") long_score_configs(TRUE) else long_score_configs(FALSE)
inference_scores <- if (mode == "main") c("avg", "max", "quad") else "avg"

jobs <- expand.grid(spec = spec_names, rep = seq_len(nsim), stringsAsFactors = FALSE)
cat(sprintf(
  "run_longitudinal_sim mode=%s nsim=%d per spec Bboot=%d ncore=%d jobs=%d\n",
  mode, nsim, Bboot, ncore, nrow(jobs)
))
flush.console()

worker <- function(j) {
  sp <- jobs$spec[j]
  rr <- jobs$rep[j]
  try(long_one_rep(
    rr, cfg, long_working_spec(sp), Bboot, score_configs, inference_scores,
    base_seed = 20260717L + match(sp, c("CC", "MC", "CM", "MM")) * 1000000L
  ), silent = TRUE)
}

t0 <- Sys.time()
results <- parallel::mclapply(seq_len(nrow(jobs)), worker, mc.cores = ncore)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
raw_path <- file.path(out_dir, paste0("longitudinal_", mode, "_raw.rds"))
saveRDS(list(mode = mode, config = cfg, jobs = jobs, results = results, elapsed = elapsed), raw_path)

all_summaries <- list()
for (sp in spec_names) {
  idx <- which(jobs$spec == sp)
  z <- results[idx]
  good <- vapply(z, is.list, logical(1))
  z <- z[good]
  cat(sprintf("spec=%s valid=%d/%d\n", sp, length(z), length(idx)))
  all_summaries[[sp]] <- long_summarize(z, cfg$tau, sp, cfg)
}

estimation <- do.call(rbind, lapply(all_summaries, `[[`, "estimation"))
screening <- do.call(rbind, lapply(all_summaries, `[[`, "screening"))
borrowing <- do.call(rbind, lapply(all_summaries, `[[`, "borrowing"))
nuisance <- do.call(rbind, lapply(all_summaries, `[[`, "nuisance"))

write.csv(estimation, file.path(out_dir, paste0("longitudinal_", mode, "_estimation.csv")), row.names = FALSE)
write.csv(screening, file.path(out_dir, paste0("longitudinal_", mode, "_screening.csv")), row.names = FALSE)
write.csv(borrowing, file.path(out_dir, paste0("longitudinal_", mode, "_borrowing.csv")), row.names = FALSE)
write.csv(nuisance, file.path(out_dir, paste0("longitudinal_", mode, "_nuisance.csv")), row.names = FALSE)

cat(sprintf("elapsed %.1f minutes; raw results: %s\n", elapsed / 60, raw_path))
print(estimation[, c("spec", "method", "bias", "rmse", "coverage", "type1", "power")], row.names = FALSE)
print(screening[, c("spec", "score", "detection", "wrong_detection")], row.names = FALSE)
