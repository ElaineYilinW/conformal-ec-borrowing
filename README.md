# Safe selective borrowing of external controls via covariate-shift-adjusted longitudinal residual-rank screening

Simulation code for the paper *"Safe selective borrowing of external controls via
covariate-shift-adjusted longitudinal residual-rank screening"* (`borrowing_summary.tex`).

We rank each external control (EC) by the **rank of its longitudinal residual
nonconformity score** among the randomized controls — **weighted** to absorb baseline
covariate shift — and borrow only the conforming ECs through a transported AIPW. The
default method is **sym-ada**: a *symmetric* two-sided adaptive screen with *CV+* rank
calibration and an *estimand-matched* score.

## Requirements

- **R ≥ 3.6** with only the **base** distribution (`stats`, `parallel`). **No CRAN
  packages are needed** — all model fits are closed-form / base `lm`/`glm.fit`.
- A multi-core machine. The scripts parallelize over Monte-Carlo replicates with
  `parallel::mclapply` (forking; Linux/macOS). On Windows use a single core.

## Quick start

```bash
git clone <your-repo-url> conformal-ec-borrowing
cd conformal-ec-borrowing
Rscript make_figs.R                 # figures (seconds)
Rscript reference_cell.R 300 500 8  # Table 1   (nsim=300, Bboot=500, 8 cores)
Rscript robustness.R     300 500 8  # robustness table
```

Every script is self-contained: it `source()`s `engine2.R` (which `source()`s
`engine.R`) from the current directory, so **run them from the repo root**.

## File map

| File | Produces | Notes |
|------|----------|-------|
| `engine.R`          | —                    | Base DGP1 (hetero binary-X), 12 nonconformity scores, transported AIPW. |
| `engine2.R`         | —                    | Builds on `engine.R`: DGP2/DGP3 generators, symmetric/non-symmetric tails, split **and** CV+ calibration, screening driver. **Sourced by everything.** |
| **`reference_cell.R`** | **Table 1** (`tab:ref` / `tab:refcell`) | 10 methods at the reference cell; writes `reference_cell_results.csv`. |
| **`robustness.R`**  | **Robustness table** (`tab:robust`) | 7 one-axis stress rows off the sym-ada/CV+ baseline; writes `robustness_results.csv`. Consolidates the old `newsetup6/7/8`. |
| `make_figs.R`       | `fig_outcome.pdf`, `fig_score.pdf` | The mechanism figure (outcome vs. score). |
| `typeI_check.R`     | prints a 3-row table | Shows the type-I column is a **Bboot floor**, not a borrowing effect (B=100→0.068, B=999→0.050, normal→0.045). |
| `run_all.R`         | `run_all_results.csv` | Full design-space grid (~5900 rows); large, optional. |
| `audit.R`, `reweight_check.R`, `logistic_ipw_check.R`, `calib_check.R` | print checks | Back specific claims: selection bias is estimator-invariant; g-formula = IPW = AIPW; weighted conformal $p$ is super-uniform. |
| `borrowing_summary.tex` | the paper | `\includegraphics` expects `fig_*.pdf` in the same dir. |

## Reproducing every number in the paper

All commands take positional args `nsim Bboot ncore [Delta]`. **Set `ncore` to the
cores you allocate.** `Delta` defaults to 8.

```bash
# --- Table 1 (both the concise doc and the blueprint) ---
Rscript reference_cell.R 300 500 <ncore>     # -> reference_cell_results.csv
#   rows: rct, full, oracle, quad, max, c2raw, c2src, c2sym(=sym-ada),
#         c2fix(=nonsym-fix), c2ada(=nonsym-ada)

# --- Robustness table (concise doc, tab:robust) ---
Rscript robustness.R     300 500 <ncore>     # -> robustness_results.csv
#   7 settings: baseline / severe misspec / unweighted / split /
#               C-slope c2 / C-slope Mahalanobis / C-slope max(trap)

# --- Mechanism figure ---
Rscript make_figs.R                          # -> fig_outcome.pdf, fig_score.pdf

# --- (optional) type-I decomposition discussed in the text ---
Rscript typeI_check.R                        # prints the B=100/999/normal table
```

### On the university Linux cluster

`ncore` is the third argument. Point it at your allocation and run detached:

```bash
# interactive multi-core node: use all but two cores
NC=$(( $(nproc) - 2 ))
nohup Rscript reference_cell.R 300 500 $NC > refcell.log 2>&1 &
nohup Rscript robustness.R     300 500 $NC > robust.log  2>&1 &
tail -f refcell.log        # watch progress; each prints a completion line + the table
```

SLURM batch (adjust partition/account):

```bash
#!/bin/bash
#SBATCH -J ecborrow
#SBATCH -c 32                 # cores -> becomes ncore below
#SBATCH -t 02:00:00
#SBATCH --mem=8G
module load R                 # or your cluster's R module
cd $SLURM_SUBMIT_DIR
Rscript reference_cell.R 300 500 $SLURM_CPUS_PER_TASK
Rscript robustness.R     300 500 $SLURM_CPUS_PER_TASK
```

### Runtime & scaling

Cost is dominated by the **full-pipeline bootstrap** (each replicate re-runs the
entire screen `Bboot` times), so wall-time ≈ `nsim × Bboot / ncore`. As a rough guide,
`reference_cell.R` and `robustness.R` each take **~10–30 min on 32 cores** at
`(nsim=300, Bboot=500)`, and scale close to `1/ncore`. For **publication-grade
Monte-Carlo error** (type-I SE ≈ 0.007 instead of 0.013), bump `nsim` to `1000`
(~3× longer) — the cluster makes this cheap. `run_all.R` is much heavier (hours).

### Why `Bboot = 500`

Type-I is a percentile-bootstrap CI, which carries a shared floor set by the number of
bootstrap draws: with `Bboot = 100` a *correctly calibrated* estimator already reads
~0.07 (not 0.05) purely from coarse tail quantiles — see `typeI_check.R`. `Bboot = 500`
pulls that floor close to the nominal 0.05, so the type-I column reads on an absolute
scale. It does **not** change the *relative* comparison (every method shares the floor).

## Reference packages

The engine is a faithful reimplementation cross-checked against **intFRT** (Zhu et al.,
Conformal Selective Borrowing) and **rdborrow** (Zhou et al., JRSS-A 2025), extended
with covariate-shift weighting and longitudinal estimand-matched scores.
