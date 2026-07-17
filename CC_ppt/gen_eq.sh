#!/bin/bash
# Render every formula / table used in the deck as a transparent high-res PNG.
# Each is a standalone LaTeX snippet -> pdflatex -> gs (pngalpha, 300 dpi).
cd "$(dirname "$0")/eq" || exit 1
PRE='\documentclass[border=8pt]{standalone}
\usepackage{amsmath,amssymb,bm,booktabs,array}
\begin{document}'
POST='\end{document}'
mk () {
  printf '%s\n%s\n%s\n' "$PRE" "$2" "$POST" > "$1.tex"
  pdflatex -interaction=nonstopmode -halt-on-error "$1.tex" >/dev/null 2>&1
  gs -q -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r300 -dGraphicsAlphaBits=4 -dTextAlphaBits=4 -sOutputFile="$1.png" "$1.pdf" 2>/dev/null
  [ -f "$1.png" ] && echo "ok $1" || echo "FAIL $1"
}

mk eq_estimand '$\displaystyle \tau_c=c^{\top}\bm\tau=\mathbb{E}_{\mathrm{R}}\!\left[c^{\top}\{\bm Y(1)-\bm Y(0)\}\right]$'
mk eq_ps '$\displaystyle e(X)=\Pr(S{=}1\mid X),\qquad \frac{f_{\mathrm{R}}(x)}{f_{\mathrm{E}}(x)}=\frac{1-\Pr(S{=}1)}{\Pr(S{=}1)}\cdot\frac{e(x)}{1-e(x)}$'
mk eq_compat '$\displaystyle \begin{aligned} \textbf{ME-}\bm{c}:&\quad \mathbb{E}\{c^{\top}\bm Y(0)\mid X,S{=}1\}=\mathbb{E}\{c^{\top}\bm Y(0)\mid X,S{=}0\}\\[2pt] \textbf{VME}:&\quad \mathbb{E}\{\bm Y(0)\mid X,S{=}1\}=\mathbb{E}\{\bm Y(0)\mid X,S{=}0\}\\[2pt] \textbf{DE}:&\quad \bm Y(0)\mid X,S{=}1\ \overset{d}{=}\ \bm Y(0)\mid X,S{=}0 \end{aligned}$'
mk eq_nest '$\displaystyle \text{DE}\ \Longrightarrow\ \text{VME}\ \Longrightarrow\ \text{ME-}c\qquad(\text{the estimand needs only ME-}c)$'
mk eq_overlap '$\displaystyle \operatorname{supp}\{X\mid S{=}1\}\subseteq\operatorname{supp}\{X\mid S{=}0\};\qquad f_{\mathrm{R}}(x)\neq f_{\mathrm{E}}(x)\ \text{is allowed}$'
mk eq_pbasic '$\displaystyle p_j=\frac{1+\sum_{i=1}^{n_C}\mathbb{I}\!\left(V_i^{\mathrm{R}}\ge V_j^{\mathrm{E}}\right)}{n_C+1}$'
mk eq_scoredef '$\displaystyle V=s(\bm Y,X;\widehat\eta)\qquad(\text{a \emph{nonconformity score}: how unusual a point looks given the fitted model }\widehat\eta)$'
mk eq_keepband '$\displaystyle \text{keep EC } j\ \Longleftrightarrow\ \gamma<p_j<1-\gamma\qquad(\text{drop when the rank falls in \emph{either} tail --- both too small and too large})$'
mk eq_exch '$\displaystyle \Pr\!\left(V^{\mathrm{R}}\le v\mid X{=}x\right)=\Pr\!\left(V^{\mathrm{E}}\le v\mid X{=}x\right)\qquad(\text{true for any fixed score, even a misspecified one})$'
mk eq_scores '$\displaystyle \bm r=\bm Y-\widehat m(X);\qquad V^{\text{contrast}}=\frac{c^{\top}\bm r}{\widehat s(X)},\qquad V^{\max}=\max_t\left|\frac{r_t}{\widehat\sigma_t(X)}\right|,\qquad V^{\text{Mah}}=\bm r^{\top}\widehat\Sigma^{-1}\bm r$'
mk eq_afunc '$\displaystyle \widehat\mu_0=a^{\top}\overline{\bm Y}_0=\sum_{t=1}^{T}a_t\,\overline{Y}_{0t} \qquad\Longrightarrow\qquad \text{take}\ \ c=a$'
mk eq_mec_a '$\displaystyle \mathbb{E}\{a^{\top}\bm Y(0)\mid X,S{=}1\}=\mathbb{E}\{a^{\top}\bm Y(0)\mid X,S{=}0\}\qquad(\text{ME-}c\ \text{must be indexed by the \emph{estimator}{}'"'"'s } a)$'
mk eq_omega '$\displaystyle \omega_0(x)=\frac{f_{\mathrm{E},0}(x)}{f_{\mathrm{R}}(x)}\qquad(\text{compatible-EC to RCT-control density ratio})$'
mk eq_wpval '$\displaystyle \widehat p_j^{\,w}=\frac{\omega_0(X_j)+\sum_i \omega_0(X_i)\,\mathbb{I}\{V_i^{\mathrm{R}}\ge V_j^{\mathrm{E}}\}}{\omega_0(X_j)+\sum_i \omega_0(X_i)}$'
mk eq_transport '$\displaystyle \mathbb{E}_{\mathrm{R}}\!\left\{\omega_0(X)\,\mathbb{I}(V\le v)\right\}=\Pr_{\mathrm{E},0}(V\le v)\qquad(\text{weighted RCT scores represent the compatible-EC score law})$'
mk eq_cv '$\displaystyle V_i=s\{\bm Y_i,X_i;\widehat\eta^{(-k)}\},\quad i\in\mathcal C_k\qquad(\text{out-of-fold score; }K\text{-fold cross-fit, CV}^{+})$'
mk eq_symscreen '$\displaystyle \text{retain}\ \{\,j:\ \widehat p^{\,w}_{|V|,j}>\gamma\,\},\qquad |V_j|=\left|c^{\top}\bm r_j\right|\quad(\text{symmetric: trim both tails equally})$'
mk eq_adaptive '$\displaystyle \widehat\gamma=\operatorname*{arg\,min}_{\gamma}\ \Big[\max\!\big\{(\widehat\tau(\gamma)-\widehat\tau_{\mathrm{RCT}})^2-\widehat{\mathrm{Var}}_b(\gamma),\,0\big\}+\widehat{\mathrm{Var}}_b(\gamma)\Big]$'
mk eq_estimator '$\displaystyle \widehat\tau=\widehat\mu_1-\big\{(1-\lambda_\gamma)\,\widehat\mu_0^{\mathrm{R}}+\lambda_\gamma\,\widehat\mu_{0,\gamma}^{\mathrm{E}}\big\},\qquad \lambda_\gamma=\frac{\widehat V_{\mathrm{R}}}{\widehat V_{\mathrm{R}}+\widehat V_{\mathrm{E},\gamma}}$'
mk eq_dgp_shift '$\displaystyle X\mid S{=}1\sim\mathrm{Bernoulli}(0.30),\qquad X\mid S{=}0\sim\mathrm{Bernoulli}(0.70)\quad(\text{covariate shift, overlap preserved})$'
mk eq_dgp_y '$\displaystyle Y_t(0)=2+0.5t+X(1+0.4t)+b+\sigma(X)\,\varepsilon_t,\qquad b\sim N(0,0.5^2),\ \ \bm\varepsilon\sim N_4(\bm 0,R),\ \ R_{uv}=0.5^{|u-v|}$'
mk eq_dgp_sig '$\displaystyle \sigma(0)=1,\quad \sigma(1)=1.8;\qquad \bm\tau=(0,0.25,0.60,1.00)^{\top},\quad \bm Y(1)=\bm Y(0)+\bm\tau$'
mk eq_dgp_contam '$\displaystyle B\sim\mathrm{Bernoulli}(\pi_B):\quad B{=}0\ \text{compatible (above)},\qquad B{=}1\ \text{incompatible: add drift }\ \Delta\cdot(\text{pattern})_t$'
mk eq_estimands '$\displaystyle c_1=(-1,0,0,1)^{\top}\ (\text{last}-\text{baseline}),\qquad c_2=\tfrac{1}{4}(1,1,1,1)^{\top}\ (\text{visit average})$'

# ---- reduced result tables (good methods + anchors only) ----
mk tab_ref '\renewcommand{\arraystretch}{1.35}\large\begin{tabular}{lcccccc}
\toprule
Method & Detection & False-excl. & Bias & RMSE & Type-I & Power\\
\midrule
RCT-only & --- & --- & $-0.00$ & $0.20$ & $0.06$ & $0.64$\\
Full EC-AIPW & --- & --- & $-0.42$ & $0.49$ & $0.42$ & $0.06$\\
Oracle (unattainable) & $1.00$ & $0.00$ & $-0.00$ & $0.17$ & $0.06$ & $0.78$\\
\midrule
\textbf{sym-ada (ours)} & $1.00$ & $0.10$ & $\bm{0.00}$ & $\bm{0.18}$ & $\bm{0.06}$ & $\bm{0.73}$\\
\bottomrule
\end{tabular}'
mk tab_robust '\renewcommand{\arraystretch}{1.35}\large\begin{tabular}{lcccccc}
\toprule
Setting (one axis changed) & Bias & RMSE & Type-I & Power & Detection & False-excl.\\
\midrule
baseline (sym-ada, CV$^+$) & $-0.00$ & $0.18$ & $0.04$ & $0.71$ & $0.99$ & $0.10$\\
severe model misspecification & $+0.00$ & $0.19$ & $0.04$ & $0.70$ & $1.00$ & $0.11$\\
unweighted (correct model) & $-0.00$ & $0.18$ & $0.05$ & $0.71$ & $1.00$ & $0.09$\\
one-time split (vs.\ CV$^+$) & $-0.00$ & $0.19$ & $0.05$ & $0.65$ & $0.99$ & $0.11$\\
\bottomrule
\end{tabular}'

rm -f *.aux *.log
echo "== done; PNGs in eq/ =="
