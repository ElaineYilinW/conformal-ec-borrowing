# Methodology and visual audit

## The statistical object

The target is a randomized-trial population treatment-effect contrast
`tau_c = E_R[c'Y(1) - c'Y(0)]`. Only the untreated mean can be supplemented
with external controls; the treated mean remains randomized-trial-only.

The proposal has two distinct stages:

1. **Screening.** Randomized controls define an empirical reference distribution
   for a longitudinal nonconformity score. External subjects with extreme ranks
   are not borrowed.
2. **Estimation.** After screening, all nuisance components are refit. Retained
   external controls are transported back to the randomized-trial covariate mix
   and combined with the randomized controls through an inverse-variance weight.

These stages use density ratios in opposite directions:

- screening: `omega_0(x) = f_E,compatible(x) / f_R(x)`;
- estimation: `rho_gamma(x) = f_R(x) / f_E,retained(x)`.

That reversal must be visible; otherwise a non-specialist can easily conclude
that the same weight is reused.

## What the screen does and does not establish

Four targets form a one-way ladder:

`conditional distribution compatibility -> vector mean compatibility ->
contrast-specific mean compatibility`.

A chosen scalar score adds a different, many-to-one diagnostic target:
distribution compatibility implies compatibility of any fixed score, but passing
one score-based screen does not imply full outcome-law exchangeability. A score
can also react to variance/covariance differences that are harmless for a mean
contrast. The document therefore needs a visual "claim boundary."

## Why the rank is useful

The fitted longitudinal model defines an ordering, not a parametric reference
law. If randomized and compatible external outcomes have the same conditional
law and the same fixed score map is used, shared model error appears in both
sources. Ranking against out-of-fold randomized-control scores empirically
calibrates an otherwise arbitrary residual scale. Model quality still controls
separation and hence detection power.

## Why weighting is needed

Randomized and compatible external controls may have different baseline mixes.
When score distributions depend on baseline covariates, the two marginal score
mixtures differ even under conditional compatibility. Weighting randomized
calibration scores by `f_E,compatible / f_R` makes their empirical mixture mimic
the compatible external covariate mix. This is calibration under covariate shift,
not adjustment for outcome incompatibility.

The practical null ratio is directly estimable from all external covariates only
under covariate-independent contamination. If contamination depends on `X`, the
compatible-external ratio is latent; this remains an open problem.

## Why the default is symmetric

When the screen and the estimand use the same signed contrast, deleting the high
tail of clean external scores also lowers the retained clean control mean. Since
the treatment effect is treated mean minus control mean, this creates positive
treatment-effect bias. Screening `|V|` with one threshold trims both clean tails
symmetrically and preserves a symmetric clean mean, at the cost of reduced power
against one-directional contamination. Thus `sym-ada` is a bias-control design,
not merely a generic two-sided outlier test.

## What CV+ contributes

A one-time split is the cleanest validity story but spends only part of an already
small randomized-control sample on calibration. K-fold cross-fitting produces an
out-of-fold score for every randomized control and fold-specific external scores.
It improves data use and simulated power. The weighted CV statistic does not yet
inherit an exact finite-sample theorem from weighted split conformal; the document
must label that gap honestly.

## What adaptive gamma contributes

For each candidate threshold, the method balances an estimated squared bias and
bootstrap variance. The squared difference from the RCT-only estimator is
de-biased by subtracting the paired bootstrap variance of that difference and
flooring at zero. The threshold is label-free, but interval validity after
selection is empirical because the point-selected threshold is reused in the
bootstrap rather than nested/reselected.

## Visual gaps in the old versions

- no trial/external-data story before notation;
- no longitudinal trajectory-to-score illustration;
- no separation of overlap screening from outcome-compatibility screening;
- no covariate-mixture illustration for weighting;
- no visual distinction between split conformal and CV+;
- no illustration of symmetric trimming versus orientation bias;
- no end-to-end diagram showing screen, refit, reverse transport, and estimation;
- no compatibility/claim-boundary diagram;
- no visual map from drift shape to score choice;
- main simulation result presented only as a dense table;
- old density plot uses overlapping fills and small legends, making group identity
  and the keep/drop rule difficult to read;
- concise PDF begins with an equation before building intuition and ends with an
  orphaned references page.

## Planned figure set

Shared by both versions:

1. hybrid-trial problem and three borrowing choices;
2. longitudinal trajectory -> residual -> score -> rank;
3. covariate-shift reweighting;
4. symmetric versus one-sided selection;
5. end-to-end screen-and-refit workflow;
6. four-criterion reference-cell result graphic.

Additional in the detailed version:

7. overlap gate versus outcome-compatibility gate;
8. compatibility and claim-boundary ladder;
9. split calibration versus CV+ data reuse;
10. common-reference versus source-specific standardization;
11. four drift patterns and score/estimand alignment;
12. simulation design map and open-problem boundary.
