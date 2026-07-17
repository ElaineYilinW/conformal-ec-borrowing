# Twelve-visit longitudinal conformal-selection simulation

## Bottom line

The frozen primary method (`avg`: absolute standardized average longitudinal
residual, weighted five-fold cross-fitted conformal rank) met every prespecified
gate in 200 Monte Carlo replications with 500 full-pipeline subject bootstraps per
replication. It detected 93.4% of incompatible external controls, wrongly excluded
8.34% of compatible external controls, had bias 0.0157 for a true visit-12 ATE of
0.70, reduced RMSE by 11.9% relative to the RCT-only longitudinal estimator,
maintained type-I error at 0.055, and increased power from 0.555 to 0.715.

## Frozen main design

- Samples: 120 RCT treated, 60 RCT controls, and 100 external controls.
- Endpoint: RCT-population ATE at visit 12; true alternative effect `tau = 0.70`.
- Repeated outcomes: 12 post-baseline visits. The conditional mean has one set of
  covariate coefficients shared across visits and a numeric time slope of 0.05.
- Baseline adjustment: realized baseline outcome `Y0`, two continuous covariates,
  and four binary covariates, following the construction in
  `another_draft_main.md`.
- Residual covariance: compound symmetry with marginal variance 4 and correlation
  0.30. This is deliberately smaller than the variance 15 in the draft paper.
- Source shift: an explicit logistic source model creates mild covariate shift.
  Across 100 diagnostic datasets the mean maximum absolute SMD was 0.224 before
  transport and 0.052 after transport; source probabilities had approximate 1% and
  99% quantiles 0.54 and 0.74.
- Contamination: 30% of external controls are incompatible. Their mean trajectory is
  shifted downward by a covariate-dependent level plus `0.1 * time`. The mean shift
  is 6.43 for the visit-average and 6.98 at visit 12, giving standardized signal-to-
  noise ratios 5.37 and 3.49, respectively.
- Working outcome model: `nlme::gls` with a shared covariate vector, numeric time,
  REML, and `corCompSymm(~ 1 | id)`. All 12 outcomes enter estimation of the shared
  coefficients used in visit-12 standardization/AIPW.
- Screening: weighted subject-level five-fold cross-fitted residual ranks with
  threshold `gamma = 0.10`. The primary score is symmetric, so it does not delete
  only one outcome tail.
- Estimation: RCT treated outcomes always estimate the treated mean. Retained
  external controls supplement only the control mean through a transported AIPW
  estimate and an inverse-variance borrowing coefficient.
- Inference: subject bootstrap independently within RCT-treated, RCT-control, and
  external-control groups; every bootstrap reruns screening, source/outcome models,
  transport, and borrowing. Null and alternative results use an exact additive
  location shift, which is valid here because treatment adds the same constant at
  every visit and screening uses controls only.

## Main operating characteristics

| Method | Bias | Empirical SD | RMSE | Coverage | Type I | Power |
|---|---:|---:|---:|---:|---:|---:|
| RCT final only | -0.0217 | 0.3311 | 0.3310 | 0.960 | 0.040 | 0.525 |
| RCT longitudinal | -0.0239 | 0.3234 | 0.3235 | 0.945 | 0.055 | 0.555 |
| Full external borrowing | 0.5774 | 0.3107 | 0.6553 | 0.540 | 0.460 | 0.975 |
| Oracle clean-EC screen | -0.0089 | 0.2731 | 0.2725 | 0.915 | 0.085 | 0.760 |
| **Conformal avg (primary)** | **0.0157** | **0.2852** | **0.2849** | **0.945** | **0.055** | **0.715** |
| Conformal max | 0.0372 | 0.2763 | 0.2781 | 0.935 | 0.065 | 0.750 |
| Conformal diagonal quadratic | 0.0192 | 0.2818 | 0.2817 | 0.945 | 0.055 | 0.720 |

For the primary method the Monte Carlo SE is 0.0202 for bias, 0.0154 for a
type-I probability near 0.05, and 0.0319 for the observed power. Its average
bootstrap SE was 0.2809 versus empirical SD 0.2852 (ratio 0.985). Every method in
every replication had a bootstrap-valid fraction of 1.00.

The RCT-longitudinal versus RCT-final comparison isolates the value of the shared
longitudinal outcome structure: RMSE decreases by 2.3% and power rises from 0.525 to
0.555 before any external data are borrowed. Conformal selection plus borrowing
then reduces RMSE by a further 11.9% and raises power by 0.160 relative to RCT-long.

## Screening and score variants

| Score / covariance treatment | Detection | Wrong detection | Bad fraction retained |
|---|---:|---:|---:|
| **Average / internal CS** | **0.9340** | **0.0834** | **0.0303** |
| Average / raw scale | 0.9340 | 0.0834 | 0.0303 |
| Average / unweighted ranks | 0.9376 | 0.0919 | 0.0289 |
| Final visit only | 0.8515 | 0.0872 | 0.0648 |
| First-to-last contrast `c1` | 0.1169 | 0.0872 | 0.2943 |
| Max / internal diagonal | 0.9069 | 0.0890 | 0.0423 |
| Quadratic / internal diagonal | 0.9327 | 0.0833 | 0.0310 |
| Quadratic / internal AR(1) | 0.9209 | 0.0870 | 0.0362 |
| Quadratic / internal CS | 0.8607 | 0.0895 | 0.0615 |
| Quadratic / shrinkage covariance | 0.8914 | 0.0854 | 0.0483 |
| Average / source CS | 0.4650 | 0.0000 | 0.1885 |
| Max / source CS | 0.3246 | 0.0001 | 0.2265 |
| Quadratic / source CS | 0.2091 | 0.0490 | 0.2643 |

The primary screen retained on average 64.05 clean and 2.00 bad external controls,
with mean borrowing coefficient 0.507 and external effective sample size 55.3.
Clean and bad mean conformal p-values were 0.509 and 0.0449. Estimating the score
covariance from the contaminated external source masks the signal and is therefore
not recommended.

The first-to-last contrast has low detection for a structural reason: it cancels the
large shared level discrepancy and retains only the 1.1 drift contrast. Its
standardized signal-to-noise ratio is 0.465, compared with 5.37 for the average
trajectory score.

## Model misspecification sensitivity

The sensitivity DGP strengthens the X1 outcome coefficient from 0.30 to 1.20 and its
source-logit coefficient from 0.15 to 0.30. `M` means that X1 is omitted from the
corresponding working model. There are 50 replications per specification and 200
full-pipeline bootstraps per replication.

| Outcome/source specification | Avg bias | Avg RMSE | Coverage | Type I | Power | Detection | Wrong detection |
|---|---:|---:|---:|---:|---:|---:|---:|
| CC | -0.0022 | 0.2635 | 0.960 | 0.040 | 0.700 | 0.9351 | 0.0962 |
| MC: omit X1 from outcome | 0.0235 | 0.2741 | 0.940 | 0.060 | 0.780 | 0.9060 | 0.0876 |
| CM: omit X1 from source | 0.0787 | 0.2789 | 0.940 | 0.060 | 0.840 | 0.9167 | 0.0724 |
| MM: omit X1 from both | 0.0844 | 0.2633 | 0.980 | 0.020 | 0.760 | 0.9176 | 0.0859 |

This is a genuine stress test, not evidence of universal robustness: omitting the
source-confounding X1 creates noticeable positive bias even though the small
50-replication bootstrap table happens to retain nominal coverage. All sensitivity
replications and all bootstrap estimates were valid.

## No-contamination negative control

With `pi_bad = 0`, 200 point-estimation replications gave wrong-detection rates
0.0847, 0.0853, and 0.0832 for average, max, and quadratic scores. The primary
estimate had bias -0.0045 and RMSE 0.2635. Full/oracle borrowing had bias -0.0028
and RMSE 0.2546. Thus symmetric outcome-based selection did not create a material
mean shift when every external subject was compatible.

## Implementation and external checks

- Unit checks cover dimensions, score sign symmetry, p-value ranges, subject-level
  resampling sizes, exact null/alternative location shift, density-ratio direction,
  and large-sample recovery of the shared coefficients and CS covariance.
- Main estimates recovered time slope 0.04978 (truth 0.05), residual variance 4.004
  (truth 4), and CS correlation 0.2972 (truth 0.30). Outcome fallback and screen
  fallback rates were zero; the source logistic model converged in every replication.
- Official R documentation confirms that [`gls`](https://stat.ethz.ch/R-manual/R-patched/library/nlme/html/gls.html)
  fits correlated-error linear models by REML/ML and that
  [`corCompSymm`](https://stat.ethz.ch/R-manual/R-patched/library/nlme/html/corCompSymm.html)
  estimates a uniform within-group correlation.
- The weighted-rank self-mass and density-ratio direction follow
  [Conformal Prediction Under Covariate Shift](https://papers.neurips.cc/paper_files/paper/2019/hash/8fb21ee7a2207526da55a679f0332de2-Abstract.html).
- The two-layer transport/borrowing motivation and operating-characteristic
  benchmarks are consistent with the
  [JRSS-A longitudinal external-control paper](https://academic.oup.com/jrsssa/article/188/3/791/7742118).

The simulation uses estimated density ratios and a cross-fitted rank aggregation,
so the empirical calibration results should not be overstated as an exact
finite-sample theorem for this complete adaptive pipeline.

## Reproducible files

- `engine_longitudinal.R`: DGP, GLS, scores, conformal screen, transported AIPW,
  and subject bootstrap.
- `run_longitudinal_sim.R`: parallel main, misspecification, and clean runners.
- `test_longitudinal_engine.R`: deterministic implementation checks.
- `audit_longitudinal_results.R`: prespecified gate and bootstrap-quality audit.
- `results_longitudinal_main/`: 200-replication, 500-bootstrap main raw and CSVs.
- `results_longitudinal_misspec/`: four 50-replication, 200-bootstrap sensitivity
  specifications.
- `results_longitudinal_clean/`: 200-replication no-contamination negative control.

The ideal 1,000 bootstraps would have exceeded the four-hour simulation budget by
the pilot throughput estimate. The main run therefore used 500; it took 147.1
minutes. The misspecification run took 48.4 minutes and the clean check 0.3 minutes.
