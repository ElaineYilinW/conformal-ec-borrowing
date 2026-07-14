# Safe borrowing of external controls via covariate-shift-adjusted longitudinal residual-rank screening

Methodology, simulation code, and presentation materials for a hybrid-controlled-trial
design that **borrows only the external controls (ECs) that look compatible with the
randomized controls**.

Each EC is ranked by the **rank of its longitudinal residual nonconformity score** among
the randomized controls — **weighted** to absorb baseline covariate shift — and only the
conforming ECs are borrowed through a transported AIPW. The default method is **sym-ada**:
a *symmetric* two-sided adaptive screen with *CV+* rank calibration and an
*estimand-matched* score.

## Repository layout

| Path | What it is |
|------|------------|
| **`paper/`** | The write-up. `longitudinal_residual_rank_blueprint.{tex,pdf}` (full ~41 pp blueprint), `borrowing_summary.{tex,pdf}` (concise 6 pp), and `longitudinal_residual_rank_blueprint_zh.{tex,pdf}` (Chinese). All figures live in `paper/figures/`. |
| **`sim/`** | Self-contained R simulation. `engine.R` / `engine2.R` (DGPs, scores, transported AIPW), `reference_cell.R` → **Table 1**, `robustness.R` → robustness table, plus diagnostic checks. |
| **`CC_ppt/`** | Google-Slides-ready presentation deck (`CC_presentation.pptx`), its Python build script, and the compiled formula/figure images. |
| **`codex/`** | Alternative *illustrated* blueprints (concise + detailed) with a shared TikZ visual language and a methodology/visual audit note. |

`reference/` (third-party packages, reading PDFs, scratch) and `sim/exploratory/`
(exploratory runs) are kept locally but git-ignored.

## The documents

- **`paper/longitudinal_residual_rank_blueprint.pdf`** — the full methodology blueprint:
  identification assumptions, the residual-rank screen, covariate-shift weighting, CV+
  cross-fitting, the adaptive symmetric screen, the transported-AIPW estimator, and the
  simulation study.
- **`paper/borrowing_summary.pdf`** — a 6-page concise version with the headline results.
- **`CC_ppt/CC_presentation.pptx`** — a 32-slide talk covering motivation, assumptions,
  the core conformal idea, the full method, and the simulation results.
- **`codex/output/pdf/`** — the illustrated concise (7 pp) and detailed (38 pp) renderings.

## Requirements

- **R ≥ 3.6** with only the **base** distribution (`stats`, `parallel`). **No CRAN
  packages are needed** — all model fits are closed-form / base `lm`/`glm.fit`.
- A multi-core machine. The scripts parallelize over Monte-Carlo replicates with
  `parallel::mclapply` (forking; Linux/macOS). On Windows use a single core.
- For the documents: a TeX distribution (`pdflatex`). For the deck build: Python with
  `python-pptx`, and Ghostscript for PDF→PNG.

## Quick start (simulation)

Every script `source()`s `engine2.R` (which `source()`s `engine.R`) **from the current
directory**, so run them from inside `sim/`:

```bash
cd sim
Rscript make_figs.R                 # mechanism figures -> ../paper/figures/ (seconds)
Rscript reference_cell.R 300 500 8  # Table 1   (nsim=300, Bboot=500, 8 cores)
Rscript robustness.R     300 500 8  # robustness table
```

### Simulation file map

| File | Produces | Notes |
|------|----------|-------|
| `engine.R`          | —                    | Base DGP1 (hetero binary-X), 12 nonconformity scores, transported AIPW. |
| `engine2.R`         | —                    | DGP2/DGP3 generators, symmetric/non-symmetric tails, split **and** CV+ calibration, screening driver. **Sourced by everything.** |
| **`reference_cell.R`** | **Table 1** | 10 methods at the reference cell; writes `reference_cell_results.csv`. |
| **`robustness.R`**  | **Robustness table** | 7 one-axis stress rows off the sym-ada/CV+ baseline; writes `robustness_results.csv`. |
| `make_figs.R`       | `../paper/figures/fig_outcome.pdf`, `fig_pval.pdf` | The mechanism figure (outcome vs. score). |
| `typeI_check.R`     | prints a 3-row table | Shows the type-I column is a **Bboot floor**, not a borrowing effect (B=100→0.068, B=999→0.050, normal→0.045). |
| `run_all.R`         | `run_all_results.csv` | Full design-space grid (~5900 rows); large, optional. |
| `audit.R`, `reweight_check.R`, `logistic_ipw_check.R`, `calib_check.R` | print checks | Back specific claims: selection bias is estimator-invariant; g-formula = IPW = AIPW; weighted conformal *p* is super-uniform. |

## Reproducing the tables

All commands take positional args `nsim Bboot ncore [Delta]`. **Set `ncore` to the cores
you allocate.** `Delta` defaults to 8.

```bash
cd sim
Rscript reference_cell.R 300 500 <ncore>     # Table 1  -> reference_cell_results.csv
Rscript robustness.R     300 500 <ncore>     # robustness table -> robustness_results.csv
Rscript typeI_check.R                        # the B=100/999/normal type-I decomposition
```

### On a Linux cluster

```bash
cd sim
NC=$(( $(nproc) - 2 ))
nohup Rscript reference_cell.R 300 500 $NC > refcell.log 2>&1 &
nohup Rscript robustness.R     300 500 $NC > robust.log  2>&1 &
tail -f refcell.log
```

SLURM batch (adjust partition/account):

```bash
#!/bin/bash
#SBATCH -J ecborrow
#SBATCH -c 32
#SBATCH -t 02:00:00
#SBATCH --mem=8G
module load R
cd $SLURM_SUBMIT_DIR/sim
Rscript reference_cell.R 300 500 $SLURM_CPUS_PER_TASK
Rscript robustness.R     300 500 $SLURM_CPUS_PER_TASK
```

Cost is dominated by the **full-pipeline bootstrap** (each replicate re-runs the entire
screen `Bboot` times), so wall-time ≈ `nsim × Bboot / ncore`. `reference_cell.R` and
`robustness.R` each take **~10–30 min on 32 cores** at `(nsim=300, Bboot=500)`. For
publication-grade Monte-Carlo error, bump `nsim` to `1000`. `run_all.R` is much heavier.

### Why `Bboot = 500`

Type-I is a percentile-bootstrap CI, which carries a shared floor set by the number of
bootstrap draws: with `Bboot = 100` a *correctly calibrated* estimator already reads ~0.07
(not 0.05) purely from coarse tail quantiles — see `typeI_check.R`. `Bboot = 500` pulls
that floor close to the nominal 0.05, so the type-I column reads on an absolute scale. It
does **not** change the *relative* comparison (every method shares the floor).

## Building the documents

```bash
cd paper && pdflatex longitudinal_residual_rank_blueprint.tex   # figures resolve from figures/
cd CC_ppt && python3 build_ppt.py                               # rebuild the deck
```

## Reference packages

The engine is a faithful reimplementation cross-checked against **intFRT** (Zhu et al.,
Conformal Selective Borrowing) and **rdborrow** (Zhou et al., JRSS-A 2025), extended with
covariate-shift weighting and longitudinal estimand-matched scores. Local copies live under
`reference/` and are not redistributed here.
