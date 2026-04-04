# Theoretical Framework for Hidden Population Estimation

## Introduction

Estimating the size of unauthorized (hidden) populations is a
fundamental problem in official statistics and migration research. By
definition, these populations avoid administrative registration, so no
single data source provides a complete count. However, when multiple
*partial* administrative sources are available – for example,
border-guard apprehensions, police identifications, and social-insurance
records – statistical models can exploit cross-source variation to infer
the size of the population that escapes all sources.

The `uncounted` package implements the estimation framework of
Beresewicz and Pawlukiewicz (2020), which builds on the dual-system
ideas of Zhang (2008). The core insight is that the *observed* count of
unauthorized migrants in a given group (e.g., a country-of-origin
$\times$ year cell) is a power-law function of the *reference* (known)
population, modulated by an auxiliary registration rate. Fitting this
power-law model and inverting it yields an estimate of the total hidden
population.

This vignette describes the theory behind the model, the available
estimation methods, and the inference procedures implemented in the
package.

## The Model

### Mean structure

Let $i = 1,\ldots,n$ index observational units (e.g., country-of-origin
$\times$ year $\times$ sex cells). Denote

- $m_{i}$: the observed count of unauthorized migrants (e.g.,
  apprehensions),
- $N_{i}$: the reference population size (e.g., insured foreigners from
  the same origin),
- $n_{i}$: an auxiliary count from a second administrative source (e.g.,
  police identifications).

The conditional mean of $m_{i}$ is modelled as

$$E\left\lbrack m_{i} \mid N_{i},n_{i} \right\rbrack\; = \; N_{i}^{\alpha}\,\left( \frac{n_{i} + \gamma}{N_{i}} \right)^{\!\beta},$$

where $\alpha$, $\beta$, and $\gamma \geq 0$ are parameters to be
estimated. On the log scale, this becomes a linear model:

$$\log E\left\lbrack m_{i} \right\rbrack\; = \;\alpha\,\log N_{i}\; + \;\beta\,\log\!\left( \frac{n_{i} + \gamma}{N_{i}} \right).$$

Note that the argument of the second logarithm can be rewritten as
$\gamma + n_{i}/N_{i}$ when the ratio form is used inside the code.

### Parameter interpretation

**$\alpha$ – elasticity with respect to reference population.** The
parameter $\alpha$ governs how the expected observed count scales with
the reference population. When $\alpha = 1$, the relationship is
proportional; when $\alpha \in (0,1)$, the scaling is sublinear, meaning
that larger reference populations yield proportionally fewer detections.
Values of $\alpha$ near zero indicate weak dependence on population
size. The constraint $\alpha \in (0,1)$ is theoretically motivated but
not always imposed during estimation (see Section “Covariates” below).

**$\beta$ – elasticity with respect to registration rate.** The
parameter $\beta$ captures the association between the auxiliary
registration rate $\left( \gamma + n_{i}/N_{i} \right)$ and the observed
unauthorized count. A positive $\beta$ indicates that groups with higher
auxiliary detection rates also tend to have more observed unauthorized
migrants. This is expected when both sources partially detect the same
latent population.

**$\gamma$ – offset for zero auxiliary counts.** Many groups may have
$n_{i} = 0$ (no auxiliary detections), which would make
$\log\left( n_{i}/N_{i} \right)$ undefined. The offset $\gamma \geq 0$
ensures that the rate term remains positive. In practice, $\gamma$ is
typically small (close to zero). It can be fixed at a known value,
estimated from the data, or omitted entirely when all $n_{i} > 0$.

### Population size estimator

Given $\widehat{\alpha}$, the total unauthorized population is estimated
as

$$\widehat{\xi}\; = \;\sum\limits_{i = 1}^{n}N_{i}^{\widehat{\alpha}}.$$

When $\widehat{\alpha} < 1$ and $N_{i} > 1$ for all $i$, we have
$\widehat{\xi} < \sum_{i}N_{i}$, reflecting the fact that only a
fraction of the reference population is unauthorized. The population
size $\widehat{\xi}$ is a monotone increasing function of $\alpha$ (for
$N_{i} \geq 1$), a property exploited for confidence interval
construction.

``` r
library(uncounted)
data(irregular_migration)

## Basic Poisson model with estimated gamma
fit <- estimate_hidden_pop(
  data      = irregular_migration,
  observed  = ~ m,
  auxiliary = ~ n,
  reference_pop = ~ N,
  method    = "poisson"
)
summary(fit)
#> Unauthorized population estimation
#> Method: POISSON | vcov: HC3 
#> N obs: 1382 
#> Gamma: 0.001831 (estimated) 
#> Log-likelihood: -11380.65 
#> AIC: 22767.3  BIC: 22782.99 
#> Deviance: 20158.5 
#> 
#> Coefficients:
#>       Estimate Std. Error z value  Pr(>|z|)    
#> alpha 0.736646   0.040258  18.298 < 2.2e-16 ***
#> beta  0.575046   0.091861   6.260  3.85e-10 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>       Observed Estimate Estimate (BC) CI lower CI upper
#> (all)   27,105  290,597       290,498  127,726  660,703

## Population size estimate with 95% CI
popsize(fit)
#>   group observed estimate estimate_bc    lower    upper share_pct
#> 1 (all)    27105 290597.5    290498.3 127726.5 660703.2       100
```

## Covariates

### Covariate-varying parameters

The scalar parameters $\alpha$ and $\beta$ can be made group-specific
through design matrices. Let $\mathbf{X}_{\alpha,i}$ and
$\mathbf{X}_{\beta,i}$ be covariate row vectors for observation $i$. The
unconstrained specification is

$$\alpha_{i} = \mathbf{X}_{\alpha,i}^{\top}\mathbf{a},\qquad\beta_{i} = \mathbf{X}_{\beta,i}^{\top}\mathbf{b},$$

where $\mathbf{a}$ and $\mathbf{b}$ are coefficient vectors.

### Constrained estimation

In the unconstrained case, nothing prevents ${\widehat{\alpha}}_{i}$
from falling outside $(0,1)$ or ${\widehat{\beta}}_{i}$ from being
negative. When theoretical constraints are desired, link functions can
be applied:

$$\alpha_{i} = {logit^{- 1}}\left( \mathbf{X}_{\alpha,i}^{\top}\mathbf{a} \right) = \frac{1}{1 + \exp\left( - \mathbf{X}_{\alpha,i}^{\top}\mathbf{a} \right)},\qquad\beta_{i} = \exp\left( \mathbf{X}_{\beta,i}^{\top}\mathbf{b} \right),$$

ensuring $\alpha_{i} \in (0,1)$ and $\beta_{i} > 0$. Constrained
estimation is available for the Poisson and Negative Binomial methods
(set `constrained = TRUE`). Coefficients are then reported on the link
scale; response-scale values are shown by
[`summary()`](https://rdrr.io/r/base/summary.html).

``` r
## Alpha and beta varying by sex
fit_cov <- estimate_hidden_pop(
  data      = irregular_migration,
  observed  = ~ m,
  auxiliary = ~ n,
  reference_pop = ~ N,
  method    = "poisson",
  cov_alpha = ~ sex,
  cov_beta  = ~ sex
)
summary(fit_cov)
#> Unauthorized population estimation
#> Method: POISSON | vcov: HC3 
#> N obs: 1382 
#> Gamma: 0.004998 (estimated) 
#> Log-likelihood: -11113.89 
#> AIC: 22237.78  BIC: 22263.94 
#> Deviance: 19624.98 
#> 
#> Coefficients:
#>                   Estimate Std. Error z value  Pr(>|z|)    
#> alpha:(Intercept) 0.670239   0.077201  8.6817 < 2.2e-16 ***
#> alpha:sexMale     0.075696   0.085822  0.8820 0.3777710    
#> beta:(Intercept)  0.528611   0.140139  3.7721 0.0001619 ***
#> beta:sexMale      0.082153   0.139500  0.5889 0.5559223    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>            Observed Estimate Estimate (BC) CI lower CI upper
#> sex=Female    4,512   55,913        55,837   12,446  250,514
#> sex=Male     22,593  198,638       198,553   78,555  501,855

## Constrained Poisson: alpha in (0,1), beta > 0
fit_constr <- estimate_hidden_pop(
  data      = irregular_migration,
  observed  = ~ m,
  auxiliary = ~ n,
  reference_pop = ~ N,
  method    = "poisson",
  cov_alpha = ~ sex,
  constrained = TRUE
)
summary(fit_constr)
#> Unauthorized population estimation
#> Method: POISSON | vcov: HC3 
#> N obs: 1382 
#> Gamma: 0.005432 (estimated) 
#> Log-likelihood: -11133.17 
#> AIC: 22274.33  BIC: 22295.26 
#> Deviance: 19663.53 
#> 
#> Coefficients (link scale: logit for alpha, log for beta):
#>                   Estimate Std. Error z value  Pr(>|z|)    
#> alpha:(Intercept)  0.82810    0.24580  3.3691 0.0007542 ***
#> alpha:sexMale      0.21610    0.18499  1.1682 0.2427387    
#> beta              -0.52089    0.17947 -2.9025 0.0037025 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Response-scale parameters (alpha in (0,1), beta > 0):
#>   Alpha (response scale):
#>             alpha SE(alpha)
#> sex=Female 0.6960    0.0520
#> sex=Male   0.7397    0.0407
#>   Beta (response scale):
#>    beta SE(beta)
#> 1 0.594   0.1066
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>            Observed Estimate Estimate (BC) CI lower CI upper
#> sex=Female    4,512   72,385        72,355   25,688  203,803
#> sex=Male     22,593  186,088       186,035   81,314  425,626
```

## Estimation Methods

The `uncounted` package implements four estimation methods. All share
the same mean structure; they differ in the objective function and
distributional assumptions.

### OLS on the log scale

Ordinary least squares applied to the log-linearized model:

$$\log m_{i} = \alpha\log N_{i} + \beta\log\!\left( \frac{n_{i} + \gamma}{N_{i}} \right) + \varepsilon_{i}.$$

When some $m_{i} = 0$, the response is replaced by
$\log\left( m_{i} + 1 \right)$. OLS is fast and transparent but ignores
the count nature of $m_{i}$ and may be biased under heteroscedasticity
(the log-transformation changes the conditional mean).

When $\gamma$ is estimated, it is profiled out via grid search: for each
candidate $\gamma$, the OLS coefficients are obtained in closed form,
and the value minimizing the residual sum of squares is selected.

### NLS on the original scale

Nonlinear least squares minimizes

$$\sum\limits_{i = 1}^{n}\left( m_{i} - N_{i}^{\alpha}\left( \frac{n_{i} + \gamma}{N_{i}} \right)^{\!\beta} \right)^{2}$$

using the Levenberg–Marquardt algorithm. This avoids the
log-transformation bias of OLS but still treats $m_{i}$ as continuous.

### Poisson pseudo-maximum likelihood (PPML)

The Poisson log-likelihood is

$$\ell({\mathbf{θ}}) = \sum\limits_{i = 1}^{n}\lbrack m_{i}\log\mu_{i} - \mu_{i}\rbrack,$$

where
$\mu_{i} = N_{i}^{\alpha_{i}}\,\left( \gamma + n_{i}/N_{i} \right)^{\beta_{i}}$.
As shown by Santos Silva and Tenreyro (2006), the Poisson PML estimator
is consistent for the parameters of the conditional mean
$E\left\lbrack m_{i} \mid N_{i},n_{i} \right\rbrack$ even if the true
distribution is not Poisson, provided the mean is correctly specified.
This makes PPML robust to heteroscedasticity and is the recommended
default in `uncounted`.

### Negative Binomial ML

The Negative Binomial extends the Poisson by adding a dispersion
parameter $\theta > 0$. The log-likelihood is

\$\$ \ell(\boldsymbol{\theta}, \theta) = \sum\_{i=1}^{n} \Bigl\[
\log\Gamma(m_i + \theta) - \log\Gamma(\theta) - \log(m_i!) +
\theta\log\theta - \theta\log(\theta + \mu_i) + m_i\log\mu_i -
m_i\log(\theta + \mu_i) \Bigr\]. \$\$

This accounts for overdispersion (variance exceeding the Poisson mean)
and can provide more efficient estimates when the data are truly
overdispersed. The parameter $\theta$ is jointly estimated alongside
$\alpha$, $\beta$, and $\gamma$.

### Gamma estimation

The offset parameter $\gamma$ is handled differently depending on the
method:

- **Poisson / NB:** $\gamma$ is jointly optimized within the likelihood
  using L-BFGS-B with box constraints (`gamma_bounds`, default
  $\left\lbrack 10^{- 10},\, 0.5 \right\rbrack$).
- **OLS:** $\gamma$ is profiled out via a grid search over the bounds,
  followed by local optimization
  ([`optimize()`](https://rdrr.io/r/stats/optimize.html)).
- **All methods:** $\gamma$ can be fixed at a known value
  (`gamma = 0.01`) or omitted entirely (`gamma = NULL`), the latter
  requiring $n_{i} > 0$ for all observations.

``` r
## Compare estimation methods
fit_ols <- estimate_hidden_pop(
  data = irregular_migration,
  observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
  method = "ols"
)

fit_pois <- estimate_hidden_pop(
  data = irregular_migration,
  observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
  method = "poisson"
)

fit_nb <- estimate_hidden_pop(
  data = irregular_migration,
  observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
  method = "nb"
)

## Extract population size estimates
popsize(fit_ols)
#>   group observed estimate estimate_bc    lower    upper share_pct
#> 1 (all)    27105 23991.01    23864.32 18430.34 30900.45       100
popsize(fit_pois)
#>   group observed estimate estimate_bc    lower    upper share_pct
#> 1 (all)    27105 290597.5    290498.3 127726.5 660703.2       100
popsize(fit_nb)
#>   group observed estimate estimate_bc    lower   upper share_pct
#> 1 (all)    27105  1259565     1231405 731347.6 2073376       100
```

## Inference

### Robust standard errors

All estimation methods report heteroscedasticity-consistent (HC) robust
standard errors computed via the sandwich estimator. The general form is

$$\widehat{Var}\left( \widehat{\mathbf{θ}} \right) = \left( \mathbf{J}^{\top}\mathbf{J} \right)^{- 1}\,\mathbf{J}^{\top}\,{diag}\left( {\widetilde{e}}_{1}^{2},\ldots,{\widetilde{e}}_{n}^{2} \right)\,\mathbf{J}\,\left( \mathbf{J}^{\top}\mathbf{J} \right)^{- 1},$$

where $\mathbf{J}$ is the Jacobian (or score matrix) and
${\widetilde{e}}_{i}$ are adjusted residuals. The adjustment depends on
the HC variant:

| Variant | Residual adjustment                                                                     |
|---------|-----------------------------------------------------------------------------------------|
| HC0     | ${\widetilde{e}}_{i} = e_{i}$ (no correction)                                           |
| HC1     | ${\widetilde{e}}_{i} = e_{i}\sqrt{n/(n - p)}$                                           |
| HC2     | ${\widetilde{e}}_{i} = e_{i}/\sqrt{1 - h_{ii}}$                                         |
| HC3     | ${\widetilde{e}}_{i} = e_{i}/\left( 1 - h_{ii} \right)$ (default)                       |
| HC4     | ${\widetilde{e}}_{i} = e_{i}/\left( 1 - h_{ii} \right)^{\delta_{i}/2}$                  |
| HC5     | ${\widetilde{e}}_{i} = e_{i}/\left( 1 - h_{ii} \right)^{\min{(\delta_{i},h_{\max})}/2}$ |

Here $h_{ii}$ are the hat-matrix diagonal elements (leverages) and
$\delta_{i}$ is a data-dependent exponent. HC3 (the default) is
recommended for moderate sample sizes as it provides a good balance
between finite-sample bias correction and stability.

### Delta-method confidence intervals for population size

The population size $\widehat{\xi} = \sum_{i}N_{i}^{\widehat{\alpha}}$
is a nonlinear function of $\widehat{\alpha}$. Confidence intervals are
constructed via a *monotone transformation* of the Wald interval on the
link scale.

**Step 1.** Build a Wald interval for the linear predictor
$\widehat{\eta}$ (which equals $\widehat{\alpha}$ when unconstrained, or
the logit of $\widehat{\alpha}$ when constrained):

$$\widehat{\eta}\; \pm \; z_{\alpha/2} \cdot {se}\left( \widehat{\eta} \right),$$

where
${se}\left( \widehat{\eta} \right) = \sqrt{\mathbf{x}^{\top}\mathbf{V}\mathbf{x}}$
uses the HC-robust covariance $\mathbf{V}$.

**Step 2.** Map the interval endpoints through the monotone function
$g(\alpha) = \sum_{i}N_{i}^{\alpha}$:

$${\widehat{\xi}}_{L} = \sum\limits_{i}N_{i}^{\alpha_{L}},\qquad{\widehat{\xi}}_{U} = \sum\limits_{i}N_{i}^{\alpha_{U}}.$$

Since $g$ is monotone increasing for $N_{i} \geq 1$, the resulting
interval
$\left\lbrack {\widehat{\xi}}_{L},{\widehat{\xi}}_{U} \right\rbrack$ has
the correct coverage probability (asymptotically).

**Total across groups.** When the model includes covariates in $\alpha$,
the total $\widehat{\xi} = \sum_{g}{\widehat{\xi}}_{g}$ has its own
delta-method CI computed via the gradient $\nabla_{\alpha}\xi$ and a
log-normal approximation to ensure positivity.

### Bias correction

Because $\xi(\alpha) = \sum_{i}N_{i}^{\alpha}$ is convex in $\alpha$ for
$N_{i} > 1$, Jensen’s inequality implies
$E\left\lbrack \widehat{\xi} \right\rbrack \geq \xi$. A second-order
Taylor expansion of $h(\alpha) = \sum_{i}N_{i}^{\alpha}$ around the true
$\alpha$ gives the approximate bias:

$${Bias}\left( \widehat{\xi} \right)\; \approx \;\frac{1}{2}\,\sum\limits_{i = 1}^{n}N_{i}^{\alpha}\,\left( \log N_{i} \right)^{2}\; \cdot \;\mathbf{x}^{\top}\mathbf{V}_{model}\,\mathbf{x},$$

where $\mathbf{V}_{model}$ is the model-based (homoscedastic)
variance–covariance matrix. The bias-corrected estimator is

$${\widehat{\xi}}^{BC} = \widehat{\xi} - \widehat{Bias}.$$

Model-based variance is used rather than the HC-robust variance, because
the latter can be inflated by high-leverage observations in skewed data,
leading to overcorrection. Bias correction is also applied to the CI
bounds, evaluated at their respective $\alpha$ values.

### Fractional weighted bootstrap

For finite-sample inference and as a robustness check on the analytical
CIs, the package provides the fractional weighted bootstrap (FWB) of Xu
et al. (2020) via
[`bootstrap_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/bootstrap_popsize.md).

Instead of resampling rows, FWB generates random non-negative weights
$w_{1},\ldots,w_{n}$ from a ${Dirichlet}(1,\ldots,1)$ distribution and
refits the model under these weights. This avoids duplicate-row issues
with discrete data and is naturally suited to weighted GLMs.

**Cluster bootstrap.** When observations are not independent (e.g., the
same country observed across years and sexes), a cluster variable can be
specified. In this case, all observations within the same cluster
receive the same FWB weight, preserving within-cluster correlation.

**Bootstrap CIs.** Two types are available:

- *Percentile*: the $\alpha/2$ and $(1 - \alpha/2)$ quantiles of the
  bootstrap distribution.
- *Bias-corrected percentile (BC)*: adjusts the quantile probabilities
  for median bias using the standard BC formula.

The bootstrap median is the default point estimate because the
transformation $\xi = \sum N_{i}^{\alpha}$ is convex, so the bootstrap
distribution of $\widehat{\xi}$ is right-skewed and the bootstrap mean
tends to overestimate $\xi$.

``` r
## Fit a model with country-level grouping
fit <- estimate_hidden_pop(
  data      = irregular_migration,
  observed  = ~ m,
  auxiliary = ~ n,
  reference_pop = ~ N,
  method    = "poisson",
  countries = ~ country
)

## Cluster bootstrap with 199 replicates
boot_result <- bootstrap_popsize(
  fit,
  R       = 49,
  cluster = ~ country_code,
  seed    = 42
)
boot_result
#> Bootstrap population size estimation
#> R = 49 | CI type: perc | Point estimate: median | Converged: 49 / 49 
#> Cluster bootstrap
#> 95% CI
#> 
#>   Point estimate: bootstrap median (recommended) | CI: bootstrap percentile
#>        Plugin Plugin (BC) Boot median Boot mean CI lower CI upper
#> (all) 290,597     290,498     330,467   388,022  185,961  741,945

## Access the full results table
boot_result$popsize_full
#>   group   plugin plugin_bc boot_median boot_mean    lower    upper
#> 1 (all) 290597.5  290498.3    330466.8  388022.4 185961.5 741945.3
```

## References

- Beręsewicz, M., & Pawlukiewicz, K. (2020). Estimation of the number of
  irregular foreigners in Poland using non-linear count regression
  models. arXiv preprint arXiv:2008.09407.

- Santos Silva, J. M. C. and Tenreyro, S. (2006). The log of gravity.
  *The Review of Economics and Statistics*, 88(4), 641–658.

- Xu, L., Gotwalt, C., Hong, Y., King, C. B., and Meeker, W. Q. (2020).
  Applications of the fractional-random-weight bootstrap. *The American
  Statistician*, 74(4), 345–358.

- Zhang, L.-C. (2008). *Developing methods for determining the number of
  unauthorized foreigners in Norway* (Documents 2008/11). Statistics
  Norway.
  <https://www.ssb.no/a/english/publikasjoner/pdf/doc_200811_en/doc_200811_en.pdf>

- Beręsewicz, M., & Pawlukiewicz, K. (2020). Estimation of the number of
  irregular foreigners in Poland using non-linear count regression
  models. arXiv preprint arXiv:2008.09407.
