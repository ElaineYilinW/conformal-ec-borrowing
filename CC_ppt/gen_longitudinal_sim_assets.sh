#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/eq"

PRE='\documentclass[border=8pt]{standalone}
\usepackage{amsmath,amssymb,bm,booktabs,array,xcolor}
\definecolor{deckblue}{HTML}{2F75B5}
\definecolor{deckred}{HTML}{C74440}
\begin{document}'
POST='\end{document}'

mk () {
  printf '%s\n%s\n%s\n' "$PRE" "$2" "$POST" > "$1.tex"
  pdflatex -interaction=nonstopmode -halt-on-error "$1.tex" >/dev/null
  gs -q -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r300 -dGraphicsAlphaBits=4 -dTextAlphaBits=4 \
    -sOutputFile="$1.png" "$1.pdf"
}

mk eq_long_sample_source '$\displaystyle (n_1,n_0,n_E)=(120,60,100),\quad T=12,\quad A\mid S{=}1\sim\mathrm{Bernoulli}(2/3),\ A=0\mid S{=}0;\qquad \operatorname{logit}\Pr(S{=}1\mid Z)=\alpha-0.15z(Y_0)+0.15z(X_1)-0.10z(X_2)+0.20\widetilde B_1-0.15\widetilde B_2+0.10\widetilde B_3.$'

mk eq_long_baseline '$\displaystyle U\sim N_6(0.3\bm 1,\Sigma_X),\quad (\Sigma_X)_{jk}=0.6e^{-0.3|j-k|},\quad B_j=\mathbb I(U_{j+2}>c_j),\ c=(0,0.3,-0.2,0.1);\qquad Y_0\mid U\sim N\!\left(0.5+0.3X_1+0.2X_2+0.2B_1+0.15B_2+0.10B_3+0.15B_4,\ 4\right).$'

mk eq_long_outcome '$\displaystyle \begin{aligned}Y_{it}&=\beta_0-0.50Y_{0i}+0.30X_{1i}-0.30X_{2i}+0.40B_{1i}+0.30B_{2i}+0.25B_{3i}+0.35B_{4i}+0.05t+0.70A_i-C_i\,b_t(Z_i)+\varepsilon_{it},\quad t=1,\ldots,12,\\ \bm\varepsilon_i&\sim N_{12}(\bm0,\Sigma),\qquad \Sigma=4\{(1-0.30)I_{12}+0.30\bm1\bm1^{\!\top}\}.\end{aligned}$'

mk eq_long_contam '$\displaystyle C_i\mid S_i{=}0\sim\mathrm{Bernoulli}(0.30),\quad C_i=0\mid S_i{=}1;\qquad b_t(Z)=2.5\{1.30+0.30Y_0+0.25X_1+0.25X_2+0.30B_1+0.20B_2+0.15B_3+0.25B_4\}+0.10t.$'

mk eq_long_score_p '$\displaystyle \begin{aligned}\bm r_i^{(-k)}&=\bm Y_i-\widehat{\bm m}^{(-k)}(Z_i),& V_i^{\mathrm{avg}}&=\frac{\color{deckred}{\left|\bm1^{\!\top}\bm r_i/12\right|}}{\sqrt{c^{\top}\widehat\Sigma_Rc}},\quad c=\bm1/12,\\[4pt] \widehat p_j^{,w}&=\frac{\omega(Z_j)+\sum_i\omega(Z_i)\mathbb I\{V_i\ge V_j\}}{\omega(Z_j)+\sum_i\omega(Z_i)},& \text{retain }j&\Longleftrightarrow \widehat p_j^{,w}>\gamma=0.10.\qquad\color{deckblue}{\textbf{Two-sided: }V\text{ uses }|\cdot|.}\end{aligned}$'

mk tab_long_main '\renewcommand{\arraystretch}{1.30}\large\begin{tabular}{lrrrrrr}
\toprule
Method & Detection & Wrong excl. & Bias & RMSE & Type-I & Power\\
\midrule
RCT final only & --- & --- & $-0.022$ & $0.331$ & $0.040$ & $0.525$\\
RCT longitudinal & --- & --- & $-0.024$ & $0.323$ & $0.055$ & $0.555$\\
Full EC borrowing & --- & --- & $+0.577$ & $0.655$ & $0.460$ & $0.975$\\
Oracle clean EC & $1.000$ & $0.000$ & $-0.009$ & $0.273$ & $0.085$ & $0.760$\\
\midrule
\textbf{Conformal average} & $\bm{0.934}$ & $\bm{0.083}$ & $\bm{+0.016}$ & $\bm{0.285}$ & $\bm{0.055}$ & $\bm{0.715}$\\
Conformal maximum & $0.907$ & $0.089$ & $+0.037$ & $0.278$ & $0.065$ & $0.750$\\
Conformal quadratic & $0.933$ & $0.083$ & $+0.019$ & $0.282$ & $0.055$ & $0.720$\\
\bottomrule
\end{tabular}'

rm -f eq_long_sample_source.aux eq_long_sample_source.log \
  eq_long_baseline.aux eq_long_baseline.log eq_long_outcome.aux eq_long_outcome.log \
  eq_long_contam.aux eq_long_contam.log eq_long_score_p.aux eq_long_score_p.log \
  tab_long_main.aux tab_long_main.log

echo "Longitudinal equation/table PNGs written to $(pwd)"
