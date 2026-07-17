# 1 Introduction {#introduction .unnumbered}

Randomized clinical trials (RCTs) are the gold standard for evaluating
treatment effects, but their sample sizes are often limited by ethical,
logistical, or financial constraints. In recent years, there has been
growing interest in supplementing RCT data with external control (EC)
data drawn from observational databases, natural history studies, or
historical trials to improve the precision of treatment effect
estimation \[Pocock, 1976, Valancius et al., 2024, Li et al., 2023, Zhou
et al., 2025, Hupf et al., 2021, Zhu et al., 2026\].

A fundamental challenge in using external controls is the potential for
systematic differences between trial and external populations. Even
after adjusting for observed baseline covariates, external controls may
differ from trial controls in ways that bias treatment effect
estimates-for instance, through differences in hospital environment,
cohort eligibility criteria, or unmeasured confounders that the observed
covariate set cannot capture. Most existing strategies for pooling trial
and external data-from naïve pooling to semiparametric efficient pooling
rely on a mean-exchangeability assumption between trial and external
controls. When this assumption fails, the resulting estimators can
suffer from severely inflated type I error and misleading inference
\[Pocock, 1976, Li et al., 2023\].

Our approach differs from these prior efforts in two deliberate ways.
First, where Valancius et al. \[2024\] and Li et al. \[2023\] operate
under a mean-exchangeability assumption between trial and external
controls, so that any residual discrepancy is absorbed into the nuisance
model without being named, we instead introduce an explicit bias
function $b(Z)$, defined as the conditional mean difference between
trial and external controls given $Z$, and learned directly from
longitudinal data. Second, where influence-function-based approaches in
this line of work require a constrained optimization over the efficient
influence function, our estimator is a closed-form plug-in that applies
directly to the bias-corrected outcomes, and the full semiparametric
efficient bound is recovered at rate $1 / T$ as longitudinal replication
grows. A parallel line of work by Zhou et al. \[2025\] also leverages
longitudinal outcome data in the external-control setting, but targets a
different estimand and develops weightingbased estimators for a
different design rather than a bias-modeling plug-in; our framework is
therefore complementary, occupying a different point on the tradeoff
between structural assumptions and estimator complexity. This
combination of explicit bias modeling with plug-in computation keeps the
estimator interpretable and directly implementable with existing
hybrid-trial software such as the rdborrow R package \[Shi et al.,
2025\].

Within the broader landscape of external-control methodology surveyed in
Zhu et al. \[2026\], our approach belongs to the family of
bias-model-based methods that explicitly model the conditional
discrepancy between trial and external controls \[Stuart and Rubin,
2008,

Kallus et al., 2018, Yang et al., 2025\], in contrast with
selective-borrowing strategies that retain only exchangeable external
units \[Zhu et al., 2025\], prognostic-adjustment approaches that use
external data only for nuisance fitting \[Schuler et al., 2022\],
test-then-pool methods that aggregate when compatibility tests do not
reject \[Yang et al., 2023\], estimator averaging that interpolates
between RCT-only and pooled estimators \[Cheng and Cai, 2021\], and
Bayesian dynamic borrowing through informative, commensurate, or robust
priors \[Ibrahim and Chen, 2000, Hobbs et al., 2011, Schmidli et al.,
2014\].

In this paper, we propose a semiparametric framework for incorporating
external controls that explicitly models the systematic discrepancy,
what we call the bias function, between the conditional means of trial
and external controls. Rather than assuming that external controls are
exchangeable with trial controls after covariate adjustment, we allow
the discrepancy to be nonzero but learnable. If this bias function were
completely unrestricted, external data would provide no efficiency gain
(formalized as Theorem 4 in the Supplementary Appendix); conversely, the
more parsimoniously the bias can be modeled, the greater the potential
benefit.

A key insight of our approach is that longitudinal data provides a
natural mechanism for learning the bias function. In many clinical
settings, subjects contribute repeated outcome measurements over the
course of follow-up. When the bias structure is stable across
timepoints, these repeated measurements can be pooled to estimate the
bias far more precisely than cross-sectional data alone would allow.
Intuitively, each additional timepoint provides another \"view\" of the
same underlying discrepancy between trial and external controls,
increasing the effective sample size for bias estimation. This leads to
our main theoretical finding: a plug-in estimator that uses an estimated
bias function can achieve near-oracle efficiency, that is, it performs
nearly as well as if the bias were known, when the mumber of repeated
assessments is large relative to the complexity of the bias model,

The more parsimoniously the bias function can be modeled, the faster it
can be learned from longitudinal data, and the smaller the residual gap
between the plug-in estimator and the oracle benchmark. This trade-off
between structural restriction and robustness gives practitioners a
principled way to navigate the bias-variance trade-off in external
control analyses.

# Our contributions are as follows: {#our-contributions-are-as-follows .unnumbered}

1.  We develop a plug-in doubly robust estimator for the average
    treatment effect in the trial population that incorporates
    bias-corrected external controls via a variance-optimal weighting
    scheme, automatically adapting to covariate overlap and relative
    noise across sources.

2.  We establish theoretical results showing that the variance gap
    between the plug-in estimator and the oracle benchmark decreases as
    more repeated assessments become\
    available, and we characterize this rate for both parametric and
    nonparametric bias models.

3.  We propose a two-step bias estimation procedure that leverages
    longitudinal data with sample splitting in time, ensuring that bias
    learning and treatment effect estimation use independent outcome
    information.

4.  We demonstrate through simulation studies that the proposed
    estimators maintain near-nominal type I error and achieve meaningful
    power gains over RCT-only analyses, and illustrate the method in an
    application to spinal muscular atrophy (SMA) using data from the
    Phase 3 SUNFISH trial augmented by natural history controls
    \[Mercuri et al., 2022\].\
    The remainder of the paper is organized as follows. Section 2 sets
    up the notation, assumptions, and plug-in estimator. Section 3
    develops the theoretical motivation for the plug-in choice and the
    near-oracle behavior it inherits as longitudinal information
    accumulates. Section 4 reports simulation results, Section 5
    illustrates the method on the SMA example, and Section 6 discusses
    implications and limitations.

# 2 Setup, notation, and plug-in estimator {#setup-notation-and-plug-in-estimator .unnumbered}

This section formalizes a realistic external-control setup: a randomized
trial supplemented by a control-only external sample, with subjects
measured longitudinally at $t=1, \ldots, T$. Throughout the paper, the
operational target is the trial-population average treatment effect at
the final (analysis) time point $T$; earlier time points enter only
through the bias-learning machinery developed below.

## 2.1 Notation and target estimand {#notation-and-target-estimand .unnumbered}

For each independent subject $i$, we observe
$O_{i}=\left(X_{i}, S_{i}, A_{i}, Y_{i}^{*}\right)$, where
$X_{i}=\left(X_{i 1}, \ldots, X_{i p}\right)$ are baseline covariates,
$S_{i} \in\{0,1\}$ is a data-source indicator ( $S=1$ : randomized
trial; $S=0$ : external/observational sample), $A_{i} \in\{0,1\}$ is
treatment assignment, and
$Y_{i}^{*}= \left(Y_{i 0}^{*}, Y_{i 1}^{*}, \ldots, Y_{i T}^{*}\right)$
denotes the raw longitudinal outcome with baseline measurement at $t=0$.
The external sample contains controls only, so
$\operatorname{Pr}(A=1 \mid S=0)=0$.

Define the working outcome as the change from baseline:

$$Y_{i t}=Y_{i t}^{*}-Y_{i 0}^{*}, \quad t=1, \ldots, T$$

and write $Y_{i}=\left(Y_{i 1}, \ldots, Y_{i T}\right)$. The
analysis-timepoint outcome is $Y_{i T}$. We define the adjust-\
ment set

$$Z_{i}=\left(Y_{i 0}^{*}, X_{i 1}, \ldots, X_{i p}\right)$$

which includes the baseline outcome measurement and all baseline
covariates.\
Define the trial participation propensity and trial fraction

$$p(z)=\operatorname{Pr}(S=1 \mid Z=z), \quad \rho=\operatorname{Pr}(S=1),$$

and the trial propensity score

$$e(z)=\operatorname{Pr}(A=1 \mid S=1, Z=z)$$

Target estimand. The target is the average treatment effect (ATE) in the
trial population at the analysis timepoint,

$$\tau_{\mathrm{rct}}=\mathbb{E}\left[Y_{T}^{(1)}-Y_{T}^{(0)} \mid S=1\right],$$

where $Y_{T}^{(a)}$ denotes the potential outcome at timepoint $T$ under
treatment $a \in\{0,1\}$.

Identification via outcome regressions. Throughout, the subscript $i$
indexes subjects (e.g., $Y_{i T}, X_{i}, S_{i}, A_{i}$ ); inside
expectations and variances we drop the $i$ and treat $Y_{T}$ as a
generic random variable from the source population, so that conditional
moments such as $\mathbb{E}\left(Y_{T} \mid S, A, Z\right)$ and
$\operatorname{Var}\left(Y_{T} \mid S, A, Z\right)$ refer to the
conditional distribution of $Y_{T}$ for a generic subject given
$(S, A, Z)$. Define the source- and arm-specific conditional means at
the analysis timepoint:

$$\begin{aligned}
m_{11}(z) & =\mathbb{E}\left(Y_{T} \mid S=1, A=1, Z=z\right), \\
m_{10}(z) & =\mathbb{E}\left(Y_{T} \mid S=1, A=0, Z=z\right), \\
m_{00}(z) & =\mathbb{E}\left(Y_{T} \mid S=0, A=0, Z=z\right),
\end{aligned}$$

with conditional variances

$$\begin{aligned}
& V_{11}(z)=\operatorname{Var}\left(Y_{T} \mid S=1, A=1, Z=z\right) \\
& V_{10}(z)=\operatorname{Var}\left(Y_{T} \mid S=1, A=0, Z=z\right) \\
& V_{00}(z)=\operatorname{Var}\left(Y_{T} \mid S=0, A=0, Z=z\right)
\end{aligned}$$

Under the causal-identification conditions formalized in Section 2.2
below (consistency, trial randomization, and trial-internal positivity),
$\mathbb{E}\left(Y_{T}^{(a)} \mid S=1, Z=z\right)=m_{1 a}(z)$ for
$a \in\{0,1\}$,\
and the trial ATE reduces to

$$\tau_{\mathrm{rct}}=\mathbb{E}\left[m_{11}(Z)-m_{10}(Z) \mid S=1\right]=\frac{1}{\rho} \mathbb{E}\left[p(Z)\left\{m_{11}(Z)-m_{10}(Z)\right\}\right]$$

The first form averages over the trial covariate distribution
$P_{Z \mid S=1}$; the second is the equivalent full-data representation
and makes the $p(Z) / \rho$ weight that appears in the trial-only
influence function transparent. The external-control regression $m_{00}$
does not enter identification of $\tau_{\text {rct }}$; it enters the
estimator only through the bias function in Assumption 2, which is the
mechanism by which external controls improve precision rather than
identification.

## 2.2 Assumptions for identification {#assumptions-for-identification .unnumbered}

Two structural assumptions support identification of
$\tau_{\mathrm{rct}}$ and the efficiency gains from external controls.\
Assumption 1 (Causal identification and source overlap). (a)
Consistency: $Y_{T}=Y_{T}^{(a)}$ whenever $A=a$. (b) Trial
randomization: $Y_{T}^{(a)} \perp\!\!\!\perp A \mid(S=1, Z)$ for
$a \in\{0,1\}$. (c) Trialinternal positivity: $0<e(Z)<1$ almost surely.
(d) Source overlap: $0<p(Z)<1$ almost surely.

Parts (a)-(c) hold by RCT design and are standard. Part (d) ensures that
the trial and external samples share a common covariate region within
which the bias function can be learned and information can be
transferred; in practice we assess (d) empirically using estimated
$\hat{p}(Z)$ and restrict to common support if needed.\
Assumption 2 (Bias restriction). Conditional on the adjustment set $Z$,
the mean difference between trial controls and external controls is
restricted to a known functional form indexed by an unknown parameter
$\theta \in \Theta$ :

$$\begin{equation*}
m_{10}(z)-m_{00}(z)=b(z ; \theta), \quad \theta \in \Theta \tag{1}
\end{equation*}$$

Assumption 2 encodes the substantive belief that the systematic
discrepancy between trial and external controls, after conditioning on
observed covariates, can be captured by a structured function
$b(\cdot ; \theta)$. Without such a restriction, external data cannot
reduce the variance of $\tau_{\text {rct }}$ estimators relative to the
trial-only bound (Theorem 4 in the Supplementary Appendix). In this
paper we focus on two tiers: a linear parametric specification

$$\begin{equation*}
b(z ; \theta)=\psi(z)^{\top} \theta, \quad \psi(z) \in \mathbb{R}^{q} \tag{2}
\end{equation*}$$

with a fixed, low-dimensional basis $\psi(z)$ and fixed $q$; and a
nonparametric specification in which $b(\cdot)$ belongs to a function
class $\mathcal{F}$ (e.g., splines or basis expansions) with effective
dimension $q_{n}$ that may grow with $n$. In longitudinal settings,
$b(\cdot ; \theta)$ may additionally incorporate time as an argument
while maintaining a shared parameter $\theta$ across time (see Section
2.5).

A separate working assumption-that the longitudinal outcomes
$Y_{1}, \ldots, Y_{T}$ are mutually conditionally independent given (
$S, A, Z$ )-is used only in the appendix's finite- $T$ efficiency
analysis (Theorem 5). The main-text estimator and its near-oracle
properties developed in Section 3 do not require it, so we defer its
formal statement to the Supplementary Appendix.

## 2.3 Oracle bias correction and the \"no-bias\" equivalence {#oracle-bias-correction-and-the-no-bias-equivalence .unnumbered}

A straightforward approach is to separate the structural role of the
bias restriction from the estimation of $\theta$. If $\theta$ were
known, we define the bias-corrected outcome

$$\begin{equation*}
Y_{T}^{\dagger}=Y_{T}+(1-S) b(Z ; \theta) \tag{3}
\end{equation*}$$

For trial units ( $S=1$ ) this leaves outcomes unchanged. For external
controls ( $S=0, A=0$ ),

$$\mathbb{E}\left(Y_{T}^{\dagger} \mid S=0, A=0, Z\right)=m_{00}(Z)+b(Z ; \theta)=m_{10}(Z)$$

Thus, under (1), knowing the bias function is equivalent to removing the
between-source discrepancy among controls via the deterministic shift
(3). After this oracle correction, the external controls are
mean-exchangeable with trial controls given $Z$ for the purpose of
learning the trial control mean. This equivalence clarifies why oracle
procedures can be read as \"no-bias\" external control methods applied
to $Y_{T}^{\dagger}$.

## 2.4 Influence-function and variance-optimal weighting of controls {#influence-function-and-variance-optimal-weighting-of-controls .unnumbered}

We next describe the influence-function that motivates our plug-in
estimator. Let $\eta$ denote the collection of nuisance functions
$\left(p, e, \rho, m_{11}, m_{10}, m_{00}, V_{10}, V_{00}\right)$
defined in Section 2.

Nonparametric trial efficient influence function (EIF). The
nonparametric efficient influence function for $\tau_{\text {rct }}$
(using trial information only) is the familiar AIPW form, where
$w=(s, a, z, y)$ denotes a realized observation (lowercase letters are
realizations of the\
corresponding random variables $S, A, Z, Y_{T}$ ):

$$\begin{equation*}
\phi_{0, \mathrm{rct}}(w)=\frac{s}{\rho}\left[\frac{a\left\{y-m_{11}(z)\right\}}{e(z)}-\frac{(1-a)\left\{y-m_{10}(z)\right\}}{1-e(z)}+m_{11}(z)-m_{10}(z)-\tau_{\mathrm{rct}}\right] \tag{4}
\end{equation*}$$

This expression makes clear that the only cell where control outcomes
enter (4) is $(S, A)=(1,0)$ through the residual $Y_{T}-m_{10}(Z)$.
External controls can reduce variance by contributing information about
the control mean, but only to the extent that the discrepancy
$m_{10}-m_{00}$ is handled.

Restriction-model correction term. Under (2), the restriction model
implies that the orthogonal complement of the mean-function tangent
space consists of functions of the form

$$w \mapsto h(r ; \eta) \nu(z)\{y-m(r)\}$$

where $m(r)=m(s, a, z)$ is the cellwise mean on the observed support and

$$\begin{equation*}
h(r ; \eta)=(1-a)\left\{\frac{s}{p(z)\{1-e(z)\}}-\frac{1-s}{1-p(z)}\right\} \tag{5}
\end{equation*}$$

The factor $h(r ; \eta)$ is a signed two-cell weight that contrasts
trial controls and external controls after adjusting for the trial
participation propensity $p(z)$ and for randomization through $e(z)$.
The function $\nu(\cdot)$ controls how strongly residual information is
transferred between sources.

Oracle influence function and the optimal $\nu_{\text {orc }}(z)$. If
$\theta$ were known (equivalently, if one works in the oracle-corrected
no-bias control model), the variance-minimizing choice of $\nu(z)$ has a
closed form:

$$\begin{equation*}
\nu_{\mathrm{orc}}(z)=-\frac{V_{10}(z) /\{\rho(1-e(z))\}}{V_{10}(z) /\{p(z)(1-e(z))\}+V_{00}(z) /\{1-p(z)\}} \tag{6}
\end{equation*}$$

The oracle influence function is then

$$\begin{equation*}
\phi_{\mathrm{orc}}(w)=\phi_{0, \mathrm{ret}}(w)-h(r ; \eta) y_{\mathrm{orc}}(z)\{y-m(r)\} . \tag{7}
\end{equation*}$$

Expression (6) has a transparent interpretation: the denominator is the
effective conditional variance of the signed residual contrast across
the two control cells, and the numerator is the contribution of trial
controls to the covariance between $\phi_{0, \text { rct }}$ and that
contrast, As a result, $\nu_{\text {orc }}(z)$ implements a generalized
least-squares combination of trial and external control\
information, automatically downweighting (i) regions with poor overlap
(extreme $p(z)$ ) and (ii) sources with large conditional variance.

## 2.5 Longitudinal bias modeling {#longitudinal-bias-modeling .unnumbered}

Suppose each subject contributes repeated outcomes
$\left(Y_{1}, \ldots, Y_{T}\right)$ at timepoints $t=1, \ldots, T$. The
estimand targets the analysis timepoint $Y_{T}$; earlier timepoints
serve exclusively to learn the bias function. We allow time-specific
mean functions
$m_{s a}(z, t)=\mathbb{E}\left(Y_{t} \mid S=s, A=a, Z=z\right)$ and
impose a shared-parameter bias restriction

$$m_{10}(z, t)-m_{00}(z, t)=b_{t}(z ; \theta), \quad t=1, \ldots, T$$

where $\theta$ is stable across time. The same parameter $\theta$ links
the trial-control and external-control means across time, so repeated
timepoints give repeated information about $\theta$.

We package time inside a time-indexed basis $\psi_{t}(z)$ so that the
linear specification reads uniformly

$$b_{t}(z ; \theta)=\psi_{t}(z)^{\top} \theta$$

with the same $\theta \in \mathbb{R}^{q}$ across $t$. For example,
$\psi_{t}(z)=(\psi(z), t)^{\top}$ recovers the additive
cross-sectional-plus-time-drift form with no separate parameter, and
$\psi_{t}(z)=(\psi(z), t, t \cdot \psi(z))^{\top}$ allows
covariate-by-time interactions. For the nonparametric tier,
$b_{t} \in \mathcal{F}$ with an effective dimension $q_{n}$ that may
grow with $n$.

Two-step bias estimation procedure. Let
$\mathcal{I}_{b} \subseteq\{1, \ldots, T-1\}$ denote the set of
auxiliary timepoints used for bias learning, with
$L=\left|\mathcal{I}_{b}\right|$; the analysis timepoint $T$ is reserved
for treatment-effect estimation. Estimating
$b_{T}(Z)=m_{10}(Z, T)-m_{00}(Z, T)$ proceeds in two stages.

Step 1 (nuisance fits). For each $t \in \mathcal{I}_{b}$, fit flexible
regressions $\hat{m}_{10}(z, t)$ on trial controls
$\left\{i: S_{i}=1, A_{i}=0\right\}$ and $\hat{m}_{00}(z, t)$ on
external controls $\left\{i: S_{i}=0\right\}$. We recommend generalized
additive models (GAMs) with penalized splines for continuous covariates
so that the bias-response retains any nonlinear structure; if both
first-stage fits are linear, the bias-response is forced to be linear
regardless of the second-stage model.

Step 2 (working bias model). Form the bias-response
$\widetilde{b}(z, t)=\hat{m}_{10}(z, t)-\hat{m}_{00}(z, t)$ and fit the
working model by pooling auxiliary timepoints:

$$\begin{equation*}
\hat{\theta}=\arg \min _{\theta \in \Theta} \sum_{t \in \mathcal{I}_{b}} \sum_{i=1}^{n}\left\{\widetilde{b}\left(Z_{i}, t\right)-b_{t}\left(Z_{i} ; \theta\right)\right\}^{2}+\lambda J(\theta), \tag{8}
\end{equation*}$$

where $J(\theta)$ is a tier-appropriate penalty: $J \equiv 0$ for the
linear specification (ordinary least squares with
$b_{t}=\psi_{t}^{\top} \theta$ ) and a roughness penalty for the
nonparametric tier. The analysistimepoint bias estimate is then

$$\hat{b}_{T}(z)=b_{T}(z ; \hat{\theta})$$

and the bias-corrected outcome is

$$\begin{equation*}
\hat{Y}_{T}^{\dagger}=Y_{T}+(1-S) \hat{b}_{T}(Z) \tag{9}
\end{equation*}$$

The next subsection plugs $\hat{Y}_{T}^{\dagger}$ into the augmented
trial-population estimator developed in Section 2.4.

T-learner vs. meta-learner estimation of $b_{t}$. The two-step procedure
above estimates the source contrast $b_{t}(z)=m_{10, t}(z)-m_{00, t}(z)$
by fitting the two source-specific control means separately and
subtracting them; this is, in effect, a T-learner \[Künzel et al.,
2019\] applied to the conditional discrepancy between trial and external
controls. Meta-learner constructions in particular, an X-learner
\[Künzel et al., 2019\] adapted to the externalcontrol setting can
sharpen the finite-sample estimate of $b_{t}$ by using imputed
source-contrast pseudo-outcomes and by adapting to imbalance between
trial and external control samples. We adopt the simpler T-learner form
here for transparency, but the framework is compatible with meta-learner
first-stage fits at no change to the downstream plug-in.

Sample splitting in time. Restricting the bias-learning timepoints to
$\mathcal{I}_{b} \subseteq\{1, \ldots, T-1\}$ ensures that
$\hat{\theta}$ does not directly reuse the analysis-timepoint outcome
$Y_{T}$. Under the regularity working assumption deferred to the
Supplementary Appendix, this further removes the crosstime covariance
between $\hat{\theta}$ and the oracle influence-function contribution at
timepoint $T$, giving the clean variance decomposition in Section 3;
under serial correlation, a residual cross-covariance of smaller order
may remain but does not affect the leading-order conclusion.

## 2.6 Oracle-style plug-in estimator under unknown bias {#oracle-style-plug-in-estimator-under-unknown-bias .unnumbered}

With $\hat{\theta}$ from (8) and the bias-corrected outcome
$\hat{Y}_{T}^{\dagger}$ from (9), we plug into the oracle-form augmented
estimator developed in Section 2.4. As a benchmark, when longitudinal
data are unavailable, a single-timepoint special case obtains by
regressing $\left\{Y_{i T}-\hat{m}_{00}\left(Z_{i}\right)\right\}$ on
$\psi\left(Z_{i}\right)$ among trial controls,

$$\begin{equation*}
\hat{\theta}_{1 \mathrm{tp}}=\arg \min _{\theta \in \mathbf{R}^{4}} \sum_{i, S_{i}=1, A_{i}=0}\left\{Y_{i T}-\hat{m}_{00}\left(Z_{i}\right)-\psi\left(Z_{i}\right)^{\top} \theta\right\}^{2} \tag{10}
\end{equation*}$$

which reuses the analysis-timepoint outcome and is therefore not the
recommended estimator when longitudinal data are available.

Outcome regressions for treated and controls. We estimate $m_{11}(z)$
using treated trial subjects $(S=1, A=1)$. For the control mean, we work
with the bias-corrected outcome and define

$$m_{0}(z) \equiv m_{10}(z)=\mathbb{E}\left(Y_{T}^{\dagger} \mid A=0, Z=z\right)$$

which is common across trial controls and oracle-corrected external
controls. In practice we fit $\hat{m}_{0}(z)$ by regressing
$\widehat{Y}_{T}^{\dagger}$ on $Z$ using the pooled controls
$\left\{i: A_{i}=0\right\}$, thereby allowing both sources to inform the
control mean after bias correction.

Control-residual weight and variance ratio. To implement the
variance-optimal combination implicit in (6), it is convenient to work
with the conditional variance ratio

$$r(z)=\frac{V_{10}(z)}{V_{00}(z)}$$

estimated by $\hat{r}(z)=\hat{V}_{10}(z) / \hat{V}_{00}(z)$ using trial
controls and external controls, respectively. (A pragmatic alternative
sets $r(z) \equiv 1$ under approximate homoskedasticity across control
sources.) Define the resulting control-residual weight

$$\begin{equation*}
\widehat{W}^{d r}(S, A, Z)=\frac{S(1-A) \hat{p}(Z)+(1-S) \hat{p}(Z) \hat{r}(Z)}{\hat{p}(Z)\{1-\hat{e}(Z)\}+\{1-\hat{p}(Z)\} \hat{r}(Z)} \tag{11}
\end{equation*}$$

which is identically zero for treated trial units and combines trial and
external controls in a way that adapts to overlap through $\hat{p}(Z)$
and to relative noise through $\hat{r}(Z)$. Algebraically, (11)
corresponds to the oracle correction term in (7) written in a
single-weight form on the control residual; it is simply a re-expression
of the same influence-function calculation.

Final plug-in estimator. Let $\hat{\rho}$ denote the empirical fraction
of trial units (or the known design fraction when fixed by sampling).
Our oracle-style plug-in estimator is

$$\begin{equation*}
\hat{\tau}_{\text {plug }}=\frac{1}{n} \sum_{i=1}^{n} \frac{1}{\hat{\rho}}\left[S_{i}\left\{\hat{m}_{11}\left(Z_{i}\right)-\hat{m}_{0}\left(Z_{i}\right)\right\}+\frac{S_{i} A_{i}}{\hat{e}\left(Z_{i}\right)}\left\{Y_{i r}-\hat{m}_{11}\left(Z_{i}\right)\right\}-\widehat{W}^{d r}\left(S_{i}, A_{i}, Z_{i}\right)\left\{\hat{Y}_{i T}^{\dagger}-\hat{m}_{0}\left(Z_{i}\right)\right\}\right] \tag{12}
\end{equation*}$$

The estimator (12) is explicit; it is the sample analog of the oracle
influence-function representation with $\theta$ replaced by the
longitudinal $\hat{\theta}$ in (8) and with nuisance regressions
substituted by their estimates. This is the operational estimator used
in the simulation and

SMA application; Section 3 shows it attains the known-bias oracle
benchmark as longitudinal information accumulates.

## 2.7 Nuisance estimation {#nuisance-estimation .unnumbered}

Equation (12) requires estimates of the nuisance functions
$\eta=\left(p, e, \rho, m_{11}, m_{10}, m_{00}, V_{10}, V_{00}\right)$
and, through $\widehat{Y}_{T}^{\dagger}$, of the bias parameter
$\hat{\theta}$. Three practical strategies are available depending on
the application.\
(a) Cross-fitted nuisance estimation. The default is sample-splitting in
subjects: partition the dataset into $K$ folds and, for each fold $k$,
fit all nuisance functions on the complement of fold $k$ using flexible
methods (e.g., generalized additive models or random forests for the
outcome regressions and propensities) and evaluate the resulting
$\eta^{(-k)}$ on fold $k$. Cross-fitting removes own-observation bias
and, under standard rate conditions on the individual nuisance fits,
allows the plug-in to inherit the $\sqrt{n}$ asymptotics of Section 3
without requiring the nuisance estimators to lie in a Donsker class.\
(b) Without cross-fitting under Donsker-class conditions. When the
nuisance estimators are restricted to a sufficiently small class (e.g.,
a fixed parametric or low-complexity GAM specification),
empirical-process arguments yield $o_{p}\left(n^{-1 / 2}\right)$
remainder terms without sample splitting. Cross-fitting can then be
omitted, at the cost of constraining the nuisance flexibility.\
(c) Simple parametric nuisance fits. In rare-disease applications such
as the SMA example, the trial size may be too small to support stable
cross-fitted machine-learning nuisance fits, or fast bootstrap inference
may be desired. In these regimes we fit all nuisances by simple
parametric models (e.g., logistic regression for $p$ and $e$, ordinary
least squares for $m_{11}$, $m_{10}, m_{00}$ ), recomputed within each
bootstrap resample. This loses the flexibility of (a) but preserves the
bias-correction mechanism and keeps each bootstrap iteration
computationally light.

## 2.8 Inference and diagnostics {#inference-and-diagnostics .unnumbered}

We report percentile bootstrap confidence intervals obtained by
rerunning the full estimation pipeline (including nuisance estimation
and $\hat{\theta}$ ) within each bootstrap resample; cross-fitting
details are deferred to the Supplementary Appendix. As diagnostics, we
examine overlap via\
$\hat{p}(Z)$, and compare $\hat{\tau}_{\text {plug }}$ to trial-only
estimators to assess sensitivity to the external-control component.

# 3 Theoretical investigation {#theoretical-investigation .unnumbered}

## 3.1 Theoretical Motivation {#theoretical-motivation .unnumbered}

Several observations motivate the plug-in estimator developed in Section
2.

1.  Unrestricted bias yields no efficiency gain. If the bias function
    $b(z)=m_{10}(z)- m_{00}(z)$ is left completely unrestricted, the
    semiparametric efficient influence function for
    $\tau_{\text {rct }}$ reduces to the trial-only AIPW influence
    function and the efficiency bound equals the trial-only bound
    $\sigma_{\text {ret }}^{2}$ (Theorem 4 in the Supplementary
    Appendix). No regular asymptotically linear estimator that
    incorporates external controls without a structural link between
    $m_{10}$ and $m_{00}$ can attain asymptotic variance smaller than
    $\sigma_{\text {ret }}^{2} / n$. This forces some structural
    restriction on the bias function before external data can deliver
    any precision gain, and as point 4 below makes precise-further
    motivates using more than a single time point so that the
    restriction can be learned cheaply.

2.  Even with the restriction, fixed- $T$ efficiency is not easy to
    implement. Under the bias restriction (1), the semiparametric
    efficient single-timepoint estimator can be characterized by
    projecting $\phi_{0, \text { rct }}$ onto the tangent space induced
    by the restriction, yielding a closed-form weight
    $\nu^{*}(z)=\nu_{\text {orc }}(z)+\psi(z)^{\top} \lambda^{*} / \Sigma(z)$
    with an explicit Lagrange multiplier enforcing
    $\mathbb{E}\left[\nu^{*}(Z) \psi(Z)\right]=0$ (Theorem 3 in the
    Supplementary Appendix). At any fixed $T$, the longitudinal
    extension (Theorem 5 in the Supplementary Appendix) couples a single
    multiplier $\lambda^{*}$ across all $T$ time-specific efficient
    weights $\nu_{t}^{*}$ and therefore requires a joint $T$-fold
    solver. The benchmarks satisfy

$$\sigma_{\mathrm{ret}}^{2} \geq \sigma_{\mathrm{eff}}^{2} \geq \sigma_{\mathrm{orc}}^{2}$$

where $\sigma_{\text {eff }}^{2}$ is the feasible efficient bound when
$\theta$ must be learned at the analysis timepoint and
$\sigma_{\text {orc }}^{2}$ is the infeasible known-bias benchmark.\
3. The plug-in is the practical alternative. Rather than solve the fully
efficient finite- $T$ problem, we estimate the bias function from
auxiliary timepoints $t<T$ and substitute the estimate into the
oracle-form estimator at the analysis timepoint $T$. The construction is
meaningful only because the bias structure is shared across time:
earlier outcomes provide\
repeated information about the same discrepancy between trial and
external controls, while the final timepoint is reserved for
treatment-effect estimation.\
4. The variance gap to oracle vanishes with longitudinal replication.
The asymptotic regime is therefore not merely $n \rightarrow \infty$ at
fixed $T$; it is a longitudinal regime in which the bias-learning budget
grows. If $L \leq T-1$ auxiliary timepoints are used to estimate a
fixed-dimensional $\theta$, the additional first-order variance from
bias estimation is of order $1 / L$ (Theorem 1). For a nonparametric
bias class with effective dimension $q_{n}$, the corresponding excess is
of order $q_{n} / T$ (Theorem 2). Consequently, as $T$ grows relative to
the complexity of the bias model, the plug-in estimator approaches the
known-bias oracle benchmark $\sigma_{\mathrm{orc}}^{2}$. Increasingly
flexible bias models can then be used without sacrificing first-order
efficiency, provided their effective complexity grows slowly enough
relative to $T$.

The single-timepoint efficient influence function, the decomposition of
the gain $\sigma_{\text {ret }}^{2}-\sigma_{\text {eff }}^{2}$, and the
conditional-independence-based finite- $T$ formula are derived in the
Supplementary Appendix as Theorem 3, Corollary 2, and Theorem 5,
respectively. Section 3.2 below summarizes the WLS decomposition only
enough to interpret the plug-in's behavior, and Sections 3.4-3.5
formalize the vanishing-gap statements in point 4.

## 3.2 Efficiency gain as a weighted least-squares projection {#efficiency-gain-as-a-weighted-least-squares-projection .unnumbered}

Under the linear bias restriction (2), the gap
$\sigma_{\text {rct }}^{2}-\sigma_{\text {eff }}^{2}$ admits a compact
weighted least-squares (WLS) representation: it equals the residual of a
WLS projection of the oracle gain profile
$g_{\text {ore }}(z)=V_{10}(z) /[\rho\{1-e(z)\}]$ onto the bias basis
$\psi(z)$, with observation weights $1 / \Sigma(z)$, where
$\Sigma(z)=V_{10}(z) /[p(z)\{1-e(z)\}]+V_{00}(z) /[1-p(z)]$. The
Supplementary Appendix's Corollary 2 makes this explicit, building on
the efficient-IF characterization of Theorem 3. Two consequences suffice
for what follows: the gain is nonnegative, and strict gains arise only
when $g_{\text {ore }}$ varies with $z$ in directions not spanned by
$\psi$.

The plug-in estimator (12) targets the oracle-form estimating equation
with $\theta$ replaced by $\hat{\theta}$. At a fixed single timepoint,
this estimator need not attain the unknown- $\theta$ efficient bound
$\sigma_{\text {eff }}^{2}$. Its justification is instead longitudinal:
when the bias parameter or function is learned from $L$ auxiliary
timepoints, the first-stage estimation penalty is $O(1 / L)$ in
fixed-dimensional models and $O\left(q_{n} / T\right)$ in nonparametric
models. Hence the plug-in approaches the oracle benchmark
$\sigma_{\text {ore }}^{2}$ as longitudinal information accumulates, as
formalized in Sections 3.4-3.5.

## 3.3 Asymptotic linearity and the role of estimating $\theta$ {#asymptotic-linearity-and-the-role-of-estimating-theta .unnumbered}

Let $\hat{\tau}_{\text {plug }, T}$ denote the plug-in estimator at
analysis timepoint $T$ constructed from (12) with $Y=Y_{T}$, and let
$\hat{\tau}_{\text {orc }, T}$ denote the same estimator computed with
the true $\theta$ (i.e., using oracle-corrected outcomes). Under overlap
$(0<p(Z)<1$ a.s. and $0<e(Z)<1$ a.s. in the trial) and standard
cross-fitting regularity (nuisance estimation errors contribute
$o_{p}\left(n^{-1 / 2}\right)$ ), the oracle estimator satisfies

$$\sqrt{n}\left(\hat{\tau}_{\mathrm{orc}, T}-\tau_{\mathrm{rct}}\right)=\frac{1}{\sqrt{n}} \sum_{i=1}^{n} \phi_{\mathrm{orc}}\left(O_{i}\right)+o_{p}(1)$$

so
$\sqrt{n}\left(\hat{\tau}_{\text {orc }, T}-\tau_{\text {rct }}\right) \xrightarrow{d} \mathcal{N}\left(0, \sigma_{\text {orc }}^{2}\right)$.\
For the plug-in estimator, a standard functional delta-method argument
yields the central expansion below.

Proposition 1 (Plug-in expansion). Under the regularity conditions of
the appendix (Conditions (A1)-(A5)), the longitudinal plug-in estimator
satisfies

$$\begin{equation*}
\sqrt{n}\left(\hat{\tau}_{\mathrm{plug}, T}-\tau_{\mathrm{rct}}\right)=\frac{1}{\sqrt{n}} \sum_{i=1}^{n} \phi_{\mathrm{orc}}\left(O_{i}\right)+\sqrt{n} \Gamma^{\top}(\hat{\theta}-\theta)+o_{p}(1), \tag{13}
\end{equation*}$$

where $\Gamma \in \mathbb{R}^{q}$ is the first-order sensitivity of the
oracle-style influence representation to perturbations of $\theta$.

Proposition 1 is the central bridge between the oracle and the plug-in:
the plug-in equals the oracle at first order plus the cost of estimating
$\theta$. The first term contributes the oracle variance
$\sigma_{\text {orc }}^{2} / n$; the second term is the sole
contribution of bias-parameter estimation at first order and is the
\"variance gap\" that Theorems 1 and 2 below quantify. The proof is in
Appendix A. 2 (the same expansion underpins both near-oracle theorems).

## 3.4 Parametric bias models: variance gap decays with number of time points {#parametric-bias-models-variance-gap-decays-with-number-of-time-points .unnumbered}

We now consider the longitudinal setting where $\theta$ is stable across
time and is estimated by pooling information across a set of auxiliary
timepoints that is disjoint from the analysis timepoint $T$. Write $L$
for the number of auxiliary timepoints used to fit $\widehat{\theta}$;
under sample splitting in time, $L$ can be as large as $T-1$. The key
requirement is that the effective information for $\theta$ scales
linearly in $L$, so that $\hat{\theta}$ concentrates at rate
$(n L)^{-1 / 2}$.\
Theorem 1 (Near-oracle behavior with longitudinal replication:
parametric bias).

Assume the bias restriction (2) is correctly specified with
$\theta \in \mathbb{R}^{q}$ fixed. Suppose $\hat{\theta}$ is obtained by
pooling information across $L$ auriliary timepoints disjoint from the
analysis timepoint $T$ (so $L \leq T-1$ ), in a way that yields

$$\sqrt{n L}(\hat{\theta}-\theta) \xrightarrow{d} \mathcal{N}\left(0, \Sigma_{\theta}\right)$$

for some positive semidefinite matrix $\Sigma_{\theta}$. Assume
cross-fitted nuisance estimation is sufficiently accurate that the
nuisance remainder is $o_{p}\left(n^{-1 / 2}\right)$. Then

$$\sqrt{n}\left(\hat{\tau}_{\mathrm{plug}, T}-\tau_{\mathrm{rct}}\right) \xrightarrow{d} \mathcal{N}\left(0, \sigma_{\mathrm{orc}}^{2}+\frac{1}{L} \Gamma^{\top} \Sigma_{\theta} \Gamma\right),$$

and equivalently,

$$\operatorname{Var}\left(\hat{\tau}_{\mathrm{plug}, T}\right)=\frac{\sigma_{\mathrm{orc}}^{2}}{n}+\frac{C_{q}}{n L}+o\left(\frac{1}{n L}\right), \quad C_{q}=\Gamma^{\top} \Sigma_{\theta} \Gamma$$

In particular, if $L \rightarrow \infty$ as $n \rightarrow \infty$
(which is permitted whenever $T \rightarrow \infty$ ), then
$\sqrt{n}\left(\hat{\tau}_{\text {plug }, T^{-}}-\right. \tau_{\text {ret }} \xrightarrow{d} \mathcal{N}\left(0, \sigma_{\text {orc }}^{2}\right)$.

Heuristic scaling in $q$. When $q$ increases (e.g., richer bases $\psi$
), the constant $C_{q}$ typically grows with $q$, yielding the heuristic
relative efficiency scaling

$$\frac{\operatorname{Var}\left(\hat{\tau}_{\mathrm{plug}, T}\right)}{\sigma_{\mathrm{orc}}^{2} / n}=1+O\left(\frac{q}{T}\right)$$

provided the remaining nuisance components remain well behaved.

## 3.5 Nonparametric bias models: effective dimension and the $q_{n} / T$ tradeoff {#nonparametric-bias-models-effective-dimension-and-the-q_n-t-tradeoff .unnumbered}

We next allow the bias function to be estimated in a nonparametric class
whose effective complexity grows with n.

Theorem 2 (Near-oracle behavior with longitudinal replication:
nonparametric bias). Assume $b(\cdot)$ belongs to a function class
$\mathcal{F}$ that admits an approximation with effective dimension
$q_{n}$ (e.g., splines or basis expansions). Suppose the pooled
longitudinal estimator satisfies the mean-square rate

$$\mathrm{E}\left[\|\hat{b}-b\|_{L_{2}\left(P_{z}\right)}^{2}\right]=O\left(\frac{q_{n}}{n T}\right)$$

or equivalently
$\|\hat{b}-b\|_{L_{2}\left(P_{z}\right)}=O_{p}\left(\sqrt{q_{n} /(n T)}\right)$
with uniform integrability. Assume further that cross-fitted nuisance
estimation errors other than $\hat{b}$ contribute
$o_{p}\left(n^{-1 / 2}\right)$. Then the plug-in estimator satisfies

$$\operatorname{Var}\left(\hat{\tau}_{\mathrm{plug}, T}\right)=\frac{\sigma_{\mathrm{orc}}^{2}}{n}+O\left(\frac{q_{n}}{n T}\right)$$

and consequently, if $T / q_{n} \rightarrow \infty$, then

$$\sqrt{n}\left(\hat{\tau}_{\text {plug }, T}-\tau_{\text {rct }}\right) \xrightarrow{d} \mathcal{N}\left(0, \sigma_{\text {orc }}^{2}\right) .$$

Interpretation. Theorems 1 and 2 formalize the central mechanism:
longitudinal replication increases the effective sample size for
learning the bias, making the oracle-style plug-in estimator behave as
if the bias were known at the $\sqrt{n}$ scale. For small $T$, the
oracle-style choice can be less efficient than the semiparametric
efficient estimator under unknown $\theta$; however, the variance gap
induced by plug-in bias learning decreases with $T$ and vanishes when
$T$ dominates the bias-model complexity.

# 4 Simulation Study {#simulation-study .unnumbered}

We conducted a Monte Carlo simulation study to evaluate the proposed
bias-corrected external control (EC) estimator, examining empirical
bias, variability (Monte Carlo standard deviation), $95 \%$ confidence
interval coverage, and power/type-I error over a grid of treatment
effect sizes.

Relation to the SMA design. While the simulation design is anchored to
the real SMA data analyzed in Section 5, we deliberately expand it over
a range of timepoints ( $T \in\{2,12\}$ ), sample-size configurations,
and bias regimes (linear and nonlinear) to probe the theoretical
predictions of Section 3 beyond what any single applied dataset can
exercise.

## 4.1 Design, data-generating mechanism, and estimators {#design-data-generating-mechanism-and-estimators .unnumbered}

In each replicate we generated a hybrid dataset consisting of an
internal randomized clinical trial (RCT) and an external control-only
dataset. Let $S \in\{0,1\}$ indicate data source ( $S=1$ RCT, $S=0$
external), $A \in\{0,1\}$ treatment assignment, and $X$ baseline
covariates observed in both sources. The RCT used a $2: 1$ randomization
with $n_{\mathrm{trt}}$ treated and $n_{\mathrm{ctrl}}$ internal
controls, and the external sample contributed $n_{\mathrm{ec}}$
additional controls. We considered two sample-size configurations,
$\left(n_{\mathrm{trt}}: n_{\mathrm{ctr}}: n_{\mathrm{ec}}\right) \in\{(100: 50: 50),(150: 75: 75)\}$,
corresponding to total\
sample sizes $n \in\{200,300\}$. In the external dataset
$\operatorname{Pr}(A=0 \mid S=0)=1$; in the RCT treatment was assigned
with $\operatorname{Pr}(A=1 \mid S=1)=2 / 3$.

Covariates. Baseline covariates included $p=2$ continuous variables (
$X_{1}, X_{2}$ ) and four binary indicators
$\left(B_{1}, \ldots, B_{4}\right)$. The $p+4=6$ latent variables were
generated jointly as multivariate normal with common mean 0.3 and
exponentially decaying covariance

$$\Sigma_{i j}=0.6 \exp \{-0.3|i-j|\},$$

and the binary covariates were obtained by thresholding:
$B_{j}=\mathbb{1}\left\{U_{j}>c_{j}\right\}$ with
$\left(c_{1}, c_{2}, c_{3}, c_{4}\right)= (0,0.3,-0.2,0.1)$. Study
membership was assigned independently of covariates by fixing sample
sizes, so $p(Z)=\operatorname{Pr}(S=1 \mid Z)=n_{\text {rct }} / n$ (no
covariate shift by construction).

Baseline and longitudinal outcomes. For each subject, the baseline raw
outcome was drawn as\
$Y_{i 0}^{\cdot} \sim \mathcal{N}\left(\mu\left(X_{i}\right), \sigma^{2}\right), \quad \mu(X)=0.5+0.3 X_{1}+0.2 X_{2}+0.2 B_{1}+0.15 B_{2}+0.1 B_{3}+0.15 B_{4}$,\
with $\sigma^{2}=15$. Conditional on the realized baseline and
covariates, the $T$ post-baseline outcomes were generated jointly:

$$\left(Y_{i 1}^{*}, \ldots, Y_{i T}^{*}\right) \mid Y_{i 0}^{*}, X_{i}, S_{i}, A_{i} \sim \mathcal{N}_{T}\left(\mu_{i}, \Sigma_{\mathrm{cs}}\right)$$

where $\mu_{i}=\left(\mu_{i 1}, \ldots, \mu_{i T}\right)^{\top}$ with

$$\mu_{i t}=f\left(Z_{i}\right)+\tau A_{i}-\left(1-S_{i}\right) b_{t}\left(Z_{i}\right)$$

and $\Sigma$ is the $T \times T$ matrix

$$\Sigma=\sigma^{2}\left(\begin{array}{cccc}
1 & \rho & \cdots & \rho \\
\rho & 1 & \cdots & \rho \\
\vdots & \vdots & \ddots & \vdots \\
\rho & \rho & \cdots & 1
\end{array}\right)$$

with $\sigma^{2}=15$ and $\rho=0.3$. The prognostic function is

$$f(Z)=0.5 Y_{0}^{*}+0.3 X_{1}-0.3 X_{2}+0.4 B_{1}+0.3 B_{2}+0.25 B_{3}+0.35 B_{4},$$

which depends on the realized baseline $Y_{0}^{*}$.

Change from baseline and adjustment set. As defined in Section 2, the
analysis endpoint is the change from baseline
$Y_{t}=Y_{t}^{*}-Y_{0}^{*}$ and the adjustment set is
$Z= \left(Y_{0}^{*}, X_{1}, X_{2}, B_{1}, \ldots, B_{4}\right)$. Since
$Y_{0}^{*}$ is a fixed quantity in $Z$, the change score inherits the
same covariance $\Sigma$ and has conditional mean
$\mathbb{E}\left(Y_{t} \mid Z, S, A\right)=\vec{f}(Z)+\tau A-(1-S) b_{t}(Z)$,
where $\tilde{f}(Z)=f(Z)-Y_{0}^{*}$. The baseline measurement remains
prognostic for the change score via a standard regression-to-the-mean
effect.

Target estimand. The target estimand
$\tau_{\text {rct }}=\mathbb{E}\left\{Y_{T}^{(1)}-Y_{T}^{(0)} \mid S=1\right\}$
as defined in Section 2; and the treatment effect $\tau$ was varied over
a grid from 0 to 3 .

EC bias functions. To represent systematic differences between internal
and external controls, we introduced a time-indexed EC bias function

$$b_{t}(Z) \equiv \mathbb{E}\left(Y_{t} \mid S=1, A=0, Z\right)-\mathbb{E}\left(Y_{t} \mid S=0, A=0, Z\right)$$

combining a covariate-dependent component with a linear time drift. Both
settings represent realistic bias regimes: the linear setting is
recoverable by a linear working specification, while the nonlinear
setting introduces smooth quadratic and additive nonlinearities that no
linear working model can capture but that a generalized additive working
bias model handles directly. The leading coefficients are calibrated so
that the bias-induced shift is of comparable magnitude to the prognostic
signal $f(Z)$, ensuring that ignoring it has a measurable impact on
naive pooling estimators.

In the linear-bias setting (Tables 1 and 2), the bias was additive in
all components of the adjustment set:

$$b_{t}(Z)=0.5+0.3 Y_{0}^{*}+0.25 X_{1}+0.25 X_{2}+0.3 B_{1}+0.2 B_{2}+0.15 B_{3}+0.25 B_{4}+0.1 t$$

In the nonlinear-bias setting (Tables 3 and 4), the bias involved
quadratic dependence on the baseline outcome and a continuous covariate:

$$b_{t}(Z)=0.3\left(Y_{0}^{*}\right)^{2}+0.25 X_{1}^{2}+0.25 X_{2}+0.3 B_{1}+0.2 B_{2}+0.15 B_{3}+0.25 B_{4}+0.1 t$$

Because $b_{t}(Z)$ here decomposes into smooth univariate functions of
each component of $Z$ plus a linear time drift, this setting is exactly
the regime targeted by a GAM-based working bias model.

We compared the following estimators. The proposed estimators are
AIPW-external-linear-bias and AIPW-external-nonlinear-bias, which
estimate $b(Z)$ under a working bias model (lin-\
ear or nonlinear, respectively) using the two-step procedure described
in Section 2.5. As benchmarks, AIPW-rct only is a doubly robust
estimator using only RCT participants, and AIPW-oracle is an augmented
estimator that uses the true $b(Z)$ as if known. Unsurprisingly, AIPW-no
bias correction augments the RCT with external controls under an
untestable mean-exchangeability assumption between internal and external
controls given $Z$ (i.e., no explicit $b(Z)$ ) and is included to
illustrate the consequences of ignoring source heterogeneity.

Nuisance estimation. The trial participation propensity
$\operatorname{Pr}(S=1 \mid Z)$ was estimated by logistic regression on
the full adjustment set. The cross-sectional outcome regressions
$\hat{m}_{s a}(z)$ at the analysis timepoint were fit by OLS on the
respective treatment-source subsets. The two-step bias estimation
procedure described in Section 2.5 was implemented with GAMbased
first-stage nuisance fits (penalized thin-plate regression splines for
each continuous covariate and timepoint, linear terms for binary
covariates, estimated by REML). For the linear working bias model, OLS
was used in the second step; for the nonlinear specification, natural
spline bases were used for continuous covariates with linear terms for
binary covariates and timepoint.

Confidence intervals were constructed using a percentile bootstrap with
$B=200$ resamples; within each resample the full estimation pipeline
(including nuisance fitting) was rerun as in the code; cross-fitting
details are deferred to the Supplementary Appendix. Power was computed
as the proportion of replicates in which the two-sided $95 \%$ bootstrap
CI excluded 0 ; under $\tau=0$, this corresponds to empirical type-I
error.

## 4.2 Simulation settings and research questions {#simulation-settings-and-research-questions .unnumbered}

We organized the Monte Carlo study around two bias-generating mechanisms
(linear versus nonlinear), two sample-size configurations
$\left(n_{\mathrm{trt}}: n_{\mathrm{ctr}}: n_{\mathrm{ec}}\right) \in\{(100: 50: 50)$,
(150:75:75)}, and the amount of longitudinal information available to
learn that bias (the number of repeated assessments $T \in\{2,12\}$ ).
Results are reported for all four combinations of bias type and
sample-size configuration.

Common setup. Unless stated otherwise, we used the data-generating
mechanism described above with $p=2$ continuous and four binary
covariates, randomization probability $e(X)=2 / 3$ in the RCT, and no
covariate shift by construction (study membership independent of $Z$ ).
The EC bias followed the time-indexed form $b_{t}(Z)$ described above,
combining a covariatedependent component and a systematic time
component, thereby allowing discrepancies between internal and external
controls to accumulate over follow-up. We evaluated empirical\
bias, Monte Carlo standard deviation, $95 \%$ confidence interval
coverage, and power/type-I error over a grid of treatment effects
$\tau \in[0,3]$.

Linear-bias setting. We first generated EC bias using the linear
specification to assess whether the proposed bias-corrected estimators
(with either a linear or nonlinear working bias model) recover valid
inference when the bias structure is simple and consistent with the
parametric restrictions imposed in Li et al. \[2023\] and Valancius
[et.al](http://et.al). \[2024\]. The primary outputs are power curves
across $\tau$ and empirical type-I error at $\tau=0$. This setting
addresses: (i) whether bias correction eliminates the extreme type-I
error inflation of uncorrected pooling, (ii) how close the proposed
procedures come to the oracle benchmark, and (iii) the magnitude of
efficiency/power gains relative to RCT-only estimators.

Nonlinear-bias setting. We next generated EC bias using the nonlinear
specification to target robustness and functional-form sensitivity. In
particular, this contrasts the behavior of the linear working bias model
versus the GAM-based working bias model when the truth is nonlinear. The
corresponding research questions are: (i) how much bias and
miscalibration results from imposing a linear working bias model under
nonlinear truth, (ii) whether the nonlinear working bias model restores
nominal coverage/type-I error, and (iii) whether valid bias correction
continues to translate into meaningful power gains.

Varying the number of repeated assessments $T$. To study the value of
longitudinal information for bias learning, we compared $T=2$ and $T=12$
repeated assessments while holding other aspects of the design fixed. In
both cases the target estimand is evaluated at the final assessment
$\left(Y_{T}\right)$. This comparison evaluates whether additional
repeated measurements improve estimation of the bias trajectory and
thereby improve the accuracy and precision of $\hat{\tau}$.

## 4.3 Results {#results .unnumbered}

Tables 1-4 report the operating characteristics of each estimator
separately for the linear-bias and nonlinear-bias data-generating
settings, with two tables per setting: one under the null ( $\tau=0$ )
reporting the type-I error and one under the alternative ( $\tau=1.7$ )
reporting the power. Within every table the rows are stacked into two
blocks corresponding to the sample-size configurations
$\left(n_{\mathrm{trt}}: n_{\mathrm{ctr}}: n_{\mathrm{ec}}\right)=(100: 50: 50)$
and $(150: 75: 75)$, and the columns are split between $T=2$ and $T=12$
repeated assessments. Power curves across the full grid of treatment
effects are deferred to the Supplementary Appendix.

::: center
+------------------------------+---------------------------------------------------------------+-----------------------------------------------------+
| Estimator                    | $T=2$                                                         | $T=12$                                              |
+:=============================+:==========================+:================+:================+:================+:================+:================+
|                              | Mean                      | SD              | Type I error    | Mean            | SD              | Type I error    |
|                              | $\hat{\boldsymbol{\tau}}$ |                 |                 | $\hat{\tau}$    |                 |                 |
+------------------------------+---------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| ( $n_{\text {trt }}: n_{\text {ctrl }}: n_{\text {ec }}$ ) = (100:50:50)                                                                           |
+------------------------------+---------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-no bias correction      | 0.839                     | 0.577           | 0.306           | 1.339           | 0.569           | 0.594           |
+------------------------------+---------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-oracle                  | -0.002                    | 0.573           | 0.057           | 0.011           | 0.559           | 0.048           |
+------------------------------+---------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-external-linear-bias    | 0.043                     | 0.703           | 0.051           | 0.011           | 0.630           | 0.034           |
+------------------------------+---------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-external-nonlinear-bias | 0.043                     | 0.705           | 0.047           | 0.009           | 0.634           | 0.033           |
+------------------------------+---------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-rct only                | -0.004                    | 0.704           | 0.056           | 0.008           | 0.689           | 0.047           |
+------------------------------+---------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| $\left(n_{\text {trt }}: n_{\text {ctri }}: n_{\text {ec }}\right)=(150: 75: 75)$                                                                  |
+------------------------------+---------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-no bias correction      | 0.830                     | 0.459           | 0.414           | 1.339           | 0.467           | 0.795           |
+------------------------------+---------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-oracle                  | -0.011                    | 0.455           | 0.047           | 0.006           | 0.464           | 0.056           |
+------------------------------+---------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-external-linear-bias    | 0.030                     | 0.550           | 0.048           | 0.006           | 0.522           | 0.051           |
+------------------------------+---------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-external-nonlinear-bias | 0.031                     | 0.552           | 0.047           | 0.006           | 0.523           | 0.053           |
+------------------------------+---------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-rct only                | -0.016                    | 0.565           | 0.055           | 0.006           | 0.573           | 0.057           |
+------------------------------+---------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+

: Table 1: Linear-bias setting, under the null ( $\tau=0$ ): mean of
$\hat{\tau}$, standard deviation, and type-I error of each estimator at
$T=2$ and $T=12$, stacked by sample-size configuration (
$n_{\mathrm{trt}}: n_{\mathrm{ctrl}}: n_{\mathrm{ec}}$ ).
:::

Linear-bias setting. Under linear EC bias the uncorrected hybrid
estimator AIPW-no bias correction exhibits substantial type-I error
inflation that worsens with both $T$ and $n_{\text {ec }}$ (Table 1).
With $n_{\text {ec }}=50$ its type-I error rises from 0.306 at $T=2$ to
0.594 at $T=12$; increasing the external-control sample to
$n_{\mathrm{ec}}=75$ pushes these values to 0.414 and 0.795 ,
respectively. The two proposed bias-corrected estimators (linear and
nonlinear working bias models) keep the type-I error within 0.033-0.053
across all four configurations, comparable to the RCT-only and oracle
benchmarks. Under the alternative $\tau=1.7$ (Table 2), both
bias-corrected estimators improve power over the RCT-only baseline; the
gap is modest at $T=2$ and grows at $T=12$ (e.g., power of 0.718 vs .0
.664 for RCT-only with $n_{\mathrm{ec}}=50$, and 0.895 vs .0 .843 with
$n_{\mathrm{ec}}=75$, for the linear working model), narrowing the
distance to the oracle as $T$ increases.

Nonlinear-bias setting. Under nonlinear EC bias the type-I error
inflation of the uncorrected estimator is even more severe (Table 3),
reaching 0.673 at $T=2$ and 0.721 at $T=12$ when $n_{\mathrm{ec}}=50$,
and 0.871 at $T=2$ and 0.904 at $T=12$ when $n_{\mathrm{ec}}=75$. As
before, both bias-corrected estimators maintain type-I error close to
the nominal level in every configuration. Under the alternative
$\tau=1.7$ (Table 4), the nonlinear working bias model is consistently
at least as powerful as the linear working model and approaches the
oracle as $T$ grows: at $T=12$ with $n_{\mathrm{ec}}=50$ it achieves
power 0.714 versus 0.704 for the linear working model and 0.670 for
RCT-only, and at $T=12$ with $n_{\mathrm{oc}}=75$ it reaches 0.894
versus 0.873 for

::: center
+------------------------------+-----------------------------------------------------+-----------------------------------------------------+
| Estimator                    | $T=2$                                               | $T=12$                                              |
+:=============================+:================+:================+:================+:================+:================+:================+
|                              | Mean            | SD              | Power           | Mean            | SD              | Power           |
|                              | $\hat{\tau}$    |                 |                 | $\hat{\tau}$    |                 |                 |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| $\left(n_{\text {trt }}: n_{\text {ctrl }}: n_{\text {ec }}\right)=(100: 50: 50)$                                                        |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-no bias correction      | 2.556           | 0.573           | 0.993           | 3.033           | 0.573           | 1.000           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-oracle                  | 1.716           | 0.569           | 0.849           | 1.705           | 0.563           | 0.857           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-external-linear-bias    | 1.763           | 0.701           | 0.685           | 1.703           | 0.638           | 0.718           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-external-nonlinear-bias | 1.764           | 0.702           | 0.674           | 1.703           | 0.642           | 0.707           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-rct only                | 1.717           | 0.707           | 0.672           | 1.693           | 0.695           | 0.664           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| $\left(n_{\text {trt }}: n_{\text {ctrl }}: n_{\text {ec }}\right)=(150: 75: 75)$                                                        |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-no bias correction      | 2.537           | 0.458           | 0.999           | 3.027           | 0.469           | 1.000           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-oracle                  | 1.697           | 0.455           | 0.959           | 1.695           | 0.466           | 0.962           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-external-linear-bias    | 1.745           | 0.555           | 0.875           | 1.692           | 0.524           | 0.895           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-external-nonlinear-bias | 1.747           | 0.557           | 0.871           | 1.692           | 0.525           | 0.889           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-rct only                | 1.693           | 0.567           | 0.854           | 1.691           | 0.578           | 0.843           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+

: Table 2: Linear-bias setting, under the alternative ( $T=1.7$ ): mean
of $\hat{\tau}$, standard deviation, and power of each estimator at
$T=2$ and $T=12$, stacked by sample-size configuration (
$n_{\mathrm{trt}}: n_{\mathrm{ctrl}}: n_{\mathrm{ec}}$ ).
:::

::: center
+------------------------------+-----------------------------------------------------+-----------------------------------------------------+
| Estimator                    | $T=2$                                               | $T=12$                                              |
+:=============================+:================+:================+:================+:================+:================+:================+
|                              | Mean            | SD              | Type I error    | Mean            | SD              | Type I error    |
|                              | $\hat{\tau}$    |                 |                 | $\hat{\tau}$    |                 |                 |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| ( $n_{\text {trt }}: n_{\text {ctri }}: n_{\text {ec }}$ ) $=(100: 50: 50)$                                                              |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-no bias correction      | 1.639           | 0.651           | 0.673           | 1.732           | 0.642           | 0.721           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-oracle                  | 0.017           | 0.569           | 0.050           | 0.013           | 0.561           | 0.047           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-external-linear-bias    | 0.023           | 0.731           | 0.047           | 0.012           | 0.658           | 0.033           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-external-nonlinear-bias | 0.026           | 0.702           | 0.043           | 0.015           | 0.639           | 0.035           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-rct only                | 0.016           | 0.699           | 0.052           | 0.007           | 0.692           | 0.053           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| $\left(n_{\text {trt }}: n_{\text {ctri }}: n_{\text {oc }}\right)=(150: 75: 75)$                                                        |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-no bias correction      | 1.639           | 0.517           | 0.871           | 1.750           | 0.535           | 0.904           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-oracle                  | -0.010          | 0.451           | 0.045           | 0.003           | 0.468           | 0.058           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-external-linear-bias    | 0.001           | 0.566           | 0.042           | 0.000           | 0.554           | 0.057           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-external-nonlinear-bias | 0.003           | 0.548           | 0.043           | 0.004           | 0.527           | 0.052           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-rct only                | -0.010          | 0.560           | 0.053           | -0.002          | 0.578           | 0.058           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+

: Table 3: Nonlinear-bias setting, under the null ( $\tau=0$ ): mean of
$\hat{\tau}$, standard deviation, and type-I error of each estimator at
$T=2$ and $T=12$, stacked by sample-size configuration (
$n_{\mathrm{trt}}: n_{\mathrm{ctr}}: n_{\mathrm{ec}}$ ).
:::

the linear working model and 0.843 for RCT-only, consistent with the
prediction that richer bias modeling pays off when the underlying bias
is genuinely nonlinear.

::: center
+------------------------------+-----------------------------------------------------+-----------------------------------------------------+
| Estimator                    | $T=2$                                               | $T=12$                                              |
+:=============================+:================+:================+:================+:================+:================+:================+
|                              | Mean            | SD              | Power           | Mean            | SD              | Power           |
|                              | $\hat{\tau}$    |                 |                 | $\hat{\tau}$    |                 |                 |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| $\left(n_{\text {trt }}: n_{\text {ctrl }}: n_{\text {ec }}\right)=(100: 50: 50)$                                                        |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-no bias correction      | 3.338           | 0.653           | 1.000           | 3.433           | 0.641           | 1.000           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-oracle                  | 1.714           | 0.569           | 0.849           | 1.711           | 0.558           | 0.862           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-external-linear-bias    | 1.722           | 0.731           | 0.648           | 1.713           | 0.653           | 0.704           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-external-nonlinear-bias | 1.724           | 0.702           | 0.651           | 1.714           | 0.632           | 0.714           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-rct only                | 1.712           | 0.702           | 0.671           | 1.706           | 0.689           | 0.670           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| $\left(n_{\text {trt }}: n_{\text {ctrl }}: n_{\text {ec }}\right)=(150: 75: 75)$                                                        |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-no bias correction      | 3.345           | 0.517           | 1.000           | 3.448           | 0.534           | 1.000           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-oracle                  | 1.694           | 0.455           | 0.961           | 1.698           | 0.468           | 0.961           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-external-linear-bias    | 1.705           | 0.570           | 0.845           | 1.700           | 0.554           | 0.873           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-external-nonlinear-bias | 1.707           | 0.554           | 0.858           | 1.703           | 0.527           | 0.894           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| AIPW-rct only                | 1.690           | 0.564           | 0.853           | 1.695           | 0.579           | 0.843           |
+------------------------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+

: Table 4: Nonlinear-bias setting, under the alternative $(\tau=1.7)$ :
mean of $\hat{\tau}$, standard deviation, and power of each estimator at
$T=2$ and $T=12$, stacked by sample-size configuration (
$n_{\text {trt }}: n_{\text {ctri }}: n_{\text {ec }}$ ).
:::

Effect of increasing $T$. Across all four tables, moving from $T=2$ to
$T=12$ consistently improves the bias-corrected estimators: standard
deviations decrease and power increases, while type-I error remains
controlled. The improvement is most pronounced under nonlinear bias,
where the proposed estimators at $T=2$ offer limited gains over RCT-only
but at $T=12$ meaningfully approach the oracle, consistent with the
theoretical $q / T$ scaling of the variance gap.

Effect of sample size. The two sample-size configurations
$\left(n_{\mathrm{trt}}: n_{\mathrm{ctr}}: n_{\mathrm{ec}}\right) \in\{(100: 50: 50),(150: 75: 7!$
allow us to separate the role of internal trial size from the role of
the external sample. The internal sample size primarily controls the
oracle floor: comparing the $(100: 50: 50)$ row to the (150:75:75) row,
the RCT-only and oracle standard deviations decrease by roughly the
expected $\sqrt{n_{\text {rct }}}$ factor, and the power of every
estimator increases correspondingly. The external sample size, in
contrast, controls the precision with which the bias function can be
learned and amplifies the cost of misspecification: increasing
$n_{\mathrm{ec}}$ from 50 to 75 aggravates the type-I error inflation of
the uncorrected hybrid estimator (from 0.594 to 0.795 in the linear
setting at $T=12$, and from 0.721 to 0.904 in the nonlinear setting at
$T=12$ ), while the bias-corrected estimators remain well calibrated and
gain modest additional power. Taken together, the two sample-size knobs
act roughly independently: internal $n_{\text {rct }}$ tightens the
oracle benchmark, external $n_{\mathrm{ec}}$ improves bias learning and
helps the plug-in close the gap to that benchmark-exactly the
decomposition predicted by Section 3.

# 5 SMA Example {#sma-example .unnumbered}

We illustrate the proposed bias-corrected EC estimator in an application
to spinal muscular atrophy (SMA) using an RCT augmented by a
control-only external dataset. The analysis horizon was $T=2$, and the
target estimand was the average treatment effect in the RCT population
at that final timepoint, where the outcome is the change from baseline
in MFM
$\left(Y_{T}=Y_{T}^{*}-Y_{0}^{*}=\mathrm{MFM}_{T}-\mathrm{MFM}_{0}\right)$
measured at $t=T=2, S=1$ indicates RCT participation, and $A=1$ denotes
the experimental treatment.

## 5.1 Data, covariates, and analysis pipeline {#data-covariates-and-analysis-pipeline .unnumbered}

The internal data source was the RCT, and the external controls were
drawn from a natural history study (NatHis) under eligibility criteria
harmonized to the RCT. Baseline covariates
$X_{i}=\left(X_{i 1}, \ldots, X_{i 5}\right)$ available in both sources
included:

- $X_{i 1}=$ age,

- $X_{i 2}=$ sex,

- $X_{i 3}=$ SMA type (Type II vs Type III),

- $X_{i 4}=$ scoliosis development $(1 / 0)$,

- $X_{i 5}=$ SMN2 gene copies (2 vs 3+4).

The adjustment set was defined as
$Z_{i}=\left(Y_{i 0}^{*}, X_{i 1}, \ldots, X_{i 5}\right)=\left(\mathrm{MFM}_{i 0}, X_{i 1}, \ldots, X_{i 5}\right)$,
and these variables were used for covariate adjustment and for modeling
study selection/overlap. To support the analysis, we evaluated
overlap/positivity using summaries of estimated trial participation
propensity $\hat{p}(Z)=\widehat{\operatorname{Pr}}(S=1 \mid Z)$.

The trial participation propensity $\operatorname{Pr}(S=1 \mid Z)$ was
estimated by logistic regression on the full adjustment set, and the
cross-sectional outcome regressions $\widehat{m}_{s a}(z)$ at the
analysis timepoint were fit by OLS on the respective treatment-source
subsets, matching the specification used in the simulation study. The
two-step bias estimation procedure described in Section 2.5 was
implemented with GAM-based first-stage nuisance fits (penalized
thin-plate regression splines for each continuous covariate and
timepoint, linear terms for binary covariates, estimated by REML). The
EC bias function at the analysis timepoint was modeled using two working
specifications, corresponding to the linear and nonlinear bias-corrected
estimators. The linear bias model used ordinary least squares with
covariates and timepoint as main effects, while the nonlinear bias model
used a GAM with smooth terms for continuous covariates and a
random-effect smooth for timepoint, fit by REML.

Inference used a percentile bootstrap with $B=200$ resamples; within
each resample the full estimation pipeline (including nuisance fitting)
was rerun as in the simulation code.

## 5.2 Results {#results-1 .unnumbered}

Table 5 reports the treatment effect estimates and uncertainty summaries
for the candidate estimators using external controls from a natural
history study (NatHis). All four estimators fall in a narrow range and
their $95 \%$ bootstrap intervals overlap substantially, indicating that
in this SMA dataset the systematic bias between the RCT controls and the
NatHis controls is modest after adjustment. The two bias-corrected
estimators (AIPW-external-linear-bias and AIPW-external-nonlinear-bias)
yield point estimates close to the RCT-only analysis with comparable
uncertainty, while the uncorrected hybrid (AIPW-no bias correction)
yields a slightly smaller point estimate and a somewhat narrower
interval. Because the bias signal here is small, the main role of the
bias-corrected framework in this dataset is to confirm that the
conclusions drawn from the RCT alone are not driven by a remaining
discrepancy with the external controls; the broader benefits of bias
correction are more clearly seen in the simulation study (Section 4).

::: center
  Estimator                      Mean    SD      CI Lower   CI Upper
  ------------------------------ ------- ------- ---------- ----------
  AIPW-external-linear-bias      1.970   0.675   0.418      3.166
  AIPW-external-nonlinear-bias   1.967   0.648   0.650      3.106
  AIPW-rct only                  2.089   0.676   0.908      3.546
  AIPW-no bias correction        1.732   0.647   0.502      3.058

  : Table 5: SMA example: treatment effect estimates (Mean), bootstrap
  standard deviation (SD), and $95 \%$ percentile bootstrap confidence
  intervals (CI Lower, CI Upper). EC from NatHis.
:::

# 6 Discussion {#discussion .unnumbered}

We have proposed a semiparametric framework for incorporating external
controls into randomized trial analyses via explicit bias function
modeling. The proposed methodology is analytically efficient and offers
a strong alternative to what is currently available in the literature:
the plug-in estimator has a closed-form, doubly robust structure that
combines trial and external control information through a
variance-optimal weighting scheme. A key\
practical advantage is the framework's two-stage architecture: the bias
learning step, which estimates the discrepancy between trial and
external controls from longitudinal data, is entirely separate from the
treatment effect estimation step, which applies the bias-corrected
outcomes in a standard augmented inverse-probability-weighted estimator.
The explicit biasfunction modeling distinguishes this framework from the
recent causal-inference approaches of Valancius et al. \[2024\] and Li
et al. \[2023\], which assume mean exchangeability between trial and
external controls, and it operationalizes the general
influence-function-based approach into a simple plug-in form that
practitioners can run without solving a separate efficient-score
optimization problem. Each stage can be independently checked and
validated: the bias model can be assessed through goodness-of-fit
diagnostics on the longitudinal residuals, while the treatment effect
estimator inherits the familiar doubly robust guarantees. This
modularity makes the framework straightforward to implement with
existing software and easy to adapt to different clinical settings, as
either stage can be modified (e.g., swapping in a more flexible bias
model or a different outcome regression) without altering the other. A
key advantage of the bias-learning step is that it does not require a
pre-specified functional form for the bias; instead, the analyst simply
increases the number of longitudinal timepoints to handle greater bias
complexity. Moreover, the plug-in estimators that achieve the
near-oracle property are themselves straightforward to implement,
requiring no bespoke influence-function solver.

A central theoretical contribution is the demonstration that different
levels of bias model complexity can be accommodated by increasing the
number of longitudinal timepoints. Theorems 1 and 2 formalize this
mechanism: for a parametric bias model of dimension $q$, the variance
gap between the plug-in and oracle estimators scales as $O(q / T)$,
while for a nonparametric bias class with effective dimension
$\boldsymbol{q}_{\boldsymbol{n}}$, the gap scales as
$O\left(\boldsymbol{q}_{\boldsymbol{n}} / T\right)$. In both cases, the
plug-in estimator converges to oracle efficiency as $T$ grows relative
to the bias model complexity. The practical implication is that
investigators with richer longitudinal data can safely employ more
flexible bias specifications, such as generalized additive models with
smooth nonlinear covariate effects, or richer nonparametric estimators,
without sacrificing inferential precision. Our simulation results
corroborate this tradeoff: under a correctly specified linear bias, the
linear working model closely approaches the oracle benchmark; under
nonlinear bias, the GAM-based working model maintains valid inference
while the linear specification, though still controlling type I error,
yields somewhat lower power. The $q / T$ scaling thus provides concrete
guidance for practitioners: one need not agonize over bias model
selection when sufficient longitudinal measurements are available, as
the framework tolerates flexible modeling choices that would otherwise
introduce excessive estimation variability. Empirically, all proposed
estimators maintained near-nominal type I error across our simulation
settings, and the GAM-based working model provided robustness to bias
misspecification at modest\
cost to power; we therefore recommend the GAM-based bias-corrected
estimator as a default when longitudinal replication permits, with the
linear specification reserved for settings where domain knowledge
supports a simple bias structure. While the $q / T$ guidance focuses on
longitudinal replication, the cross-sectional sample size $n$ continues
to govern the variance of the treatment effect estimator itself, so the
largest practical gains arise when both dimensions are sufficient.

Beyond the asymptotic results, the simulation makes the practical
picture concrete. Across every configuration the bias-corrected
estimators kept type-I error within roughly $0.03-0.05$, while the
uncorrected hybrid inflated it to 0.79 in the linear setting and 0.90 in
the nonlinear setting at $T=12$ with $n_{\mathrm{ec}}=75$. For
practitioners we suggests using the linear working bias model when
domain knowledge supports simpler bias, when the trial control arm is
small, or when fast bootstrap inference is needed, and the GAM-based
working model otherwise: both maintain near-nominal type-I error, so the
choice between them is mostly a power trade-off. The two sample sizes in
the design play complementary roles. The internal RCT size
$n_{\text {ret }}$ tightens the variance of every estimator, while the
external size $n_{\text {ec }}$ and $T$ controls how precisely the bias
function can be learned and amplifies the type-I cost of any remaining
bias. This is most useful in rare-disease applications, where the trial
is small and external borrowing offers the largest precision gain, as
long as we have enough longitudinal replication supports the estimation
bias model.

Several directions remain for future work. First, adaptive procedures
that select the bias model complexity in a data-driven manner, for
instance by cross-validation over the longitudinal timepoints, would
reduce the burden on the analyst. Second, extending the framework to
accommodate time-varying treatments represents a natural generalization.
In longitudinal settings where treatment may change over time,
g-computation provides the natural estimation strategy for identifying
causal effects under sequential treatment regimes; integrating the
bias-correction approach developed here with MMRM or g-computation
methods would broaden the applicability of the framework to a wider
class of clinical trial designs. Finally, the framework could be
extended to settings with multiple external data sources, each with its
own bias structure, by introducing source-specific bias functions linked
through shared or hierarchical parameters.

# Acknowledgements {#acknowledgements .unnumbered}

We would like to express our gratitude to Winnie Yeung and Tammy Mclver
who provided valuable insights and expertise that made this research
possible.

# Author Contributions {#author-contributions .unnumbered}

All authors designed the study. J.S. implemented the methodology and
performed the simulations and analyses. J.Z. and H.P. provided
supervision. J.S. drafted the manuscript, and J.Z. and H.P. contributed
to the writing and critically reviewed the manuscript. All authors
contributed to interpretation of the findings and critical revision of
the manuscript.

# Funding {#funding .unnumbered}

This work is supported by the Food and Drug Administration (FDA) of the
U.S. Department of Health and Human Services as part of a contract, BAA
75F40125, totaling $\$ 1,085,530$. The contents are those of the authors
and do not necessarily represent the official view of, nor an
endorsement by, the FDA, HHS, or the U.S. government.

# Conflict of Interest {#conflict-of-interest .unnumbered}

Drs Zhu and Pang are employees of Genentech/Roche. They own Roche
stocks. However, the published work is methodology focused.

# References {#references .unnumbered}

David Cheng and Tianxi Cai. Adaptive combination of randomized and
observational data. arXiv preprint arXiv:2111.15012, 2021. URL
<https://arxiv.org/abs/2111.15012>.

Brian P. Hobbs, Bradley P. Carlin, Sumithra J. Mandrekar, and Daniel J.
Sargent. Hierarchical commensurate and power prior models for adaptive
incorporation of historical information in clinical trials. Biometrics,
67(3):1047-1056, 2011. doi:10.1111/j.1541-0420.2011. 01564. x.

Bradley Hupf, Veronica Bunn, Jianchang Lin, and Cheng Dong. Bayesian
semiparametric meta-analytic-predictive prior for historical control
borrowing in clinical trials. Statistics in Medicine, 40(14):3385-3399,
June 2021. doi:10. 1002/sim. 8970.

Joseph G. Ibrahim and Ming-Hui Chen. Power prior distributions for
regression models. Statistical Science, 15(1):46-60, 2000. doi:10.
1214/ss/1009212673.

Nathan Kallus, Aahlad Manas Puli, and Uri Shalit. Removing hidden
confounding by experimental grounding. In Advances in Neural Information
Processing Systems, volume 31, 2018. URL
<https://arxiv.org/abs/1810.11646>.

Sören R. Künzel, Jasjeet S. Sekhon, Peter J. Bickel, and Bin Yu.
Metalearners for estimating heterogeneous treatment effects using
machine learning. Proceedings of the National Academy of Sciences,
116(10):4156-4165, 2019. doi:10.1073/pnas. 1804597116.\
Xinyu Li, Wang Miao, Fang Lu, and Xiao-Hua Zhou. Improving efficiency of
inference in clinical trials with external control data. Biometrics,
79(1):394-403, March 2023. doi:10.1111/biom. 13583.

Eugenio Mercuri, Nicolas Deconinck, Elena S. Mazzone, Andrés Nascimento,
Maryam Oskoui, Kayoko Saito, Carole Vuillerot, Giovanni Baranello, Odile
Boespflug-Tanguy, Nathalie Goemans, Janbernd Kirschner, Anna
Kostera-Pruszczyk, Laurent Servais, Marianne Gerber, Ksenija Gorni, Omar
Khwaja, Heidemarie Kletzl, Renata S. Scalco, Hannah Staunton, Wai Yin
Yeung, Carmen Martin, Paulo Fontoura, and John W. Day. Safety and
efficacy of once-daily risdiplam in type 2 and non-ambulant type 3
spinal muscular atrophy (SUNFISH part 2): A phase 3, double-blind,
randomised, placebo-controlled trial. The Lancet Neurology, 21(1):42-52,
2022. doi:10.1016/S1474-4422(21) 00367-7.\
Stuart J. Pocock. The combination of randomized and historical controls
in clinical trials. Journal of Chronic Diseases, 29(3):175-188, March
1976. doi:10.1016/0021-9681(76) 90044-8.

Heinz Schmidli, Sandro Gsteiger, Satrajit Roychoudhury, Anthony O'Hagan,
David Spiegelhalter, and Beat Neuenschwander. Robust
meta-analytic-predictive priors in clinical trials with historical
control information. Biometrics, 70(4):1023-1032, 2014.
doi:10.1111/biom. 12242.

Alejandro Schuler, David Walsh, Diana Hall, Jon Walsh, and Charles
Fisher. Increasing the efficiency of randomized trial estimates via
linear adjustment for a prognostic score. The International Journal of
Biostatistics, 18(2):329-356, 2022. doi:10.1515/ijb-2021-0072.\
Lei Shi, Herbert Pang, Chen Chen, and Jiawen Zhu. rdborrow: An R package
for causal inference incorporating external controls in randomized
controlled trials with longitudinal outcomes. Journal of
Biopharmaceutical Statistics, 35(6):1043-1066, 2025. doi:10.1080/
10543406.2025 .2489283.

Elizabeth A. Stuart and Donald B. Rubin. Matching with multiple control
groups with adjustment for group differences. Journal of Educational and
Behavioral Statistics, 33(3): 279-306, 2008.
doi:10.3102/1076998607306078.

Michael Valancius, Herbert Pang, Jiawen Zhu, Stephen R. Cole, Michele
Jonsson Funk, and Michael R. Kosorok. A causal inference framework for
leveraging external controls in hybrid trials. Biometrics,
80(4):ujae095, November 2024. doi:10.1093/biomtc/ujae095.

Shu Yang, Chenyin Gao, Donglin Zeng, and Xiaofei Wang. Elastic
integrative analysis of randomised trial and real-world data for
treatment heterogeneity estimation. Journal of the Royal Statistical
Society Series B: Statistical Methodology, 85(3):575-596, 2023.
doi:10.1093/jrsssb/qkad017.

Shu Yang, Siyi Liu, Donglin Zeng, and Xiaofei Wang. Data fusion methods
for the heterogeneity of treatment effect and confounding function.
Bernoulli, 31(4):2987-3012, 2025. doi:10.3150/24-BEJ1835.

Xiner Zhou, Jiawen Zhu, Christiana Drake, and Herbert Pang. Causal
estimators for incorporating external controls in randomized trials with
longitudinal outcomes. Journal of the Royal Statistical Society Series
A: Statistics in Society, 188(3):791-818, July 2025.
doi:10.1093/jrsssa/qnae075.

Ke Zhu, Shu Yang, and Xiaofei Wang. Enhancing statistical validity and
power in hybrid controlled trials: A randomization inference approach
with conformal selective borrowing. Proceedings of the 42nd
International Conference on Machine Learning (ICML), 2025. URL
<https://arxiv.org/abs/2410>. 11713.

Ke Zhu, Rima Izem, Peng Yang, Ying Yuan, Herbert Pang, Mark van der
Laan, Lei Nie, Birol Emir, Pallavi Mishra-Kalyani, Hana Lee, and Shu
Yang. Externally controlled trials: A review of design and borrowing
through a causal lens. arXiv preprint arXiv:2605.03282, 2026. URL
<https://arxiv.org/abs/2605>. 03282.
