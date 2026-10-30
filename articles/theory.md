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
Beręsewicz and Pawlukiewicz (2020), which builds on the dual-system
ideas of Zhang (2008). The core insight is that the *observed* count of
unauthorized migrants in a given group (e.g., a country-of-origin
$`\times`$ year cell) is a power-law function of the *reference* (known)
population, modulated by an auxiliary registration rate. Fitting this
power-law model and inverting it yields an estimate of the total hidden
population.

This vignette describes the theory behind the model, the available
estimation methods, and the inference procedures implemented in the
package.

## Theoretical model

### Model statement

Let $`M_{it}`$ denote the unobserved unauthorized population from
country $`i`$ in period $`t`$, $`N_{it}`$ the corresponding authorized
population, $`m_{it}`$ the observed unauthorized count, $`n_{it}`$ an
auxiliary observed count, and $`p_{it}`$ the detection rate. Denote

- $`M_{it}`$: the unobserved unauthorized population size,
- $`m_{it}`$: the observed count of unauthorized migrants (e.g.,
  apprehensions),
- $`N_{it}`$: the reference population size (e.g., insured foreigners
  from the same origin),
- $`n_{it}`$: an auxiliary count from a second administrative source
  (e.g., police identifications),
- $`p_{it}`$: the detection rate, i.e. the probability that an
  unauthorized migrant is observed.

Following the paper, define the theoretical unauthorized population size
as

``` math
\xi_{it} = E(M_{it} \mid N_{it}),
```

and the theoretical detection rate as

``` math
\rho_{it} = E(p_{it} \mid N_{it}, n_{it}).
```

Under the maintained assumptions that detection does not depend on the
unobserved population size once we condition on observed quantities, and
that $`n_{it}`$ provides no additional information about $`M_{it}`$
beyond $`N_{it}`$, the expected observed count factorises as

``` math
\mu_{it} = E(m_{it} \mid N_{it}, n_{it}) = \xi_{it} \rho_{it}.
```

This is the core theoretical decomposition: the observed count is the
product of a latent population component and a detection component.

For the remainder of the vignette, each observed $`(i,t)`$ cell is
treated as one observational unit, so we suppress the explicit $`t`$
index when discussing the package implementation.

## Empirical model

### Baseline parametric specification

The theoretical model does not itself impose functional forms on
$`\xi_i`$ and $`\rho_i`$. In the baseline empirical specification used
in the paper, the unauthorized population follows a power law

``` math
\xi_i = N_i^{\alpha_i},
```

and the detection rate is modelled as

``` math
\rho_i = \left(\gamma_i + \frac{n_i}{N_i}\right)^{\beta_i},
```

where $`\gamma_i \ge 0`$ is a baseline detection offset that keeps the
model well-defined even when $`n_i = 0`$. Combining the two pieces gives

``` math
\mu_i = N_i^{\alpha_i}\left(\gamma_i + \frac{n_i}{N_i}\right)^{\beta_i}.
```

Under the default `link_rho = "power"`, this is exactly the model fit by
the package. On the log scale it becomes

``` math
\log E[m_i] \;=\; \alpha_i \,\log N_i \;+\; \log \rho_i.
```

### Package extensions to the detection model

To preserve the paper’s baseline specification while allowing bounded
detection rates when desired, the package can rewrite the detection
component as

``` math
\eta_i = \beta_i \log\!\left(\gamma_i + \frac{n_i}{N_i}\right),
\qquad \rho_i = h(\eta_i),
```

with one of four links:

- `link_rho = "power"`: $`\rho_i = \exp(\eta_i)`$,
- `link_rho = "cloglog"`: $`\rho_i = 1 - \exp(-\exp(\eta_i))`$,
- `link_rho = "logit"`: $`\rho_i = (1 + \exp(-\eta_i))^{-1}`$,
- `link_rho = "probit"`: $`\rho_i = \Phi(\eta_i)`$.

The `power` link reproduces the paper’s baseline power-law detection
model exactly. The `cloglog`, `logit`, and `probit` links are package
extensions that constrain $`0 < \rho_i < 1`$.

Under every supported detection link,

``` math
\xi_i = \frac{E[m_i]}{\rho_i} = N_i^{\alpha_i}.
```

This is why the public population-size estimator remains
$`\hat\xi = \sum_i N_i^{\hat\alpha_i}`$ even when a bounded detection
link is used. The choice of `link_rho` affects $`\hat\xi`$ only
indirectly through the fitted values of $`\hat\alpha_i`$.

### Parameter interpretation

**$`\alpha`$ – Community Anchor.** The parameter $`\alpha`$ governs how
the expected observed count scales with the reference population. When
$`\alpha = 1`$, the relationship is proportional; when
$`\alpha \in (0, 1)`$, the scaling is sublinear, meaning that larger
reference populations yield proportionally fewer detections. Values of
$`\alpha`$ near zero indicate weak dependence on population size. The
constraint $`\alpha \in (0, 1)`$ is theoretically motivated but not
always imposed during estimation (see Section “Covariates” below).

**$`\beta`$ – Exposure Elasticity.** The parameter $`\beta`$ captures
the association between the auxiliary registration rate
$`(\gamma + n_i / N_i)`$ and the observed unauthorized count. A positive
$`\beta`$ indicates that groups with higher auxiliary detection rates
also tend to have more observed unauthorized migrants. This is expected
when both sources partially detect the same latent population.

**$`\gamma`$ – baseline detection rate.** Many groups may have
$`n_i = 0`$ (no auxiliary detections), which would make
$`\log(n_i / N_i)`$ undefined. The offset $`\gamma \geq 0`$ ensures that
the rate term remains positive. In practice, $`\gamma`$ is typically
small (close to zero). It can be fixed at a known value, estimated from
the data, or omitted entirely when all $`n_i > 0`$.

### Population size estimator

Given $`\hat\alpha`$, the total unauthorized population is estimated as

``` math
\hat\xi \;=\; \sum_{i=1}^{n} N_i^{\hat\alpha}.
```

When $`\hat\alpha < 1`$ and $`N_i > 1`$ for all $`i`$, we have
$`\hat\xi < \sum_i N_i`$, reflecting the fact that only a fraction of
the reference population is unauthorized. The population size
$`\hat\xi`$ is a monotone increasing function of $`\alpha`$ (for
$`N_i \geq 1`$). Confidence intervals for $`\hat\xi`$ are obtained by
the delta method (see the inference section below).

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
#> Method: POISSON | estimator: MLE | link_rho: power | vcov: HC3 
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
#>   group observed estimate estimate_bc    lower  upper share_pct
#> 1 (all)    27105 290597.4    290498.3 127726.4 660703       100
```

## Covariates

### Covariate-varying parameters

The scalar parameters $`\alpha`$ and $`\beta`$ can be made
group-specific through design matrices. Let $`\mathbf{X}_{\alpha,i}`$
and $`\mathbf{X}_{\beta,i}`$ be covariate row vectors for observation
$`i`$. The unconstrained specification is

``` math
\alpha_i = \mathbf{X}_{\alpha,i}^\top \mathbf{a}, \qquad
\beta_i  = \mathbf{X}_{\beta,i}^\top \mathbf{b},
```

where $`\mathbf{a}`$ and $`\mathbf{b}`$ are coefficient vectors.

### Constrained estimation

In the unconstrained case, nothing prevents $`\hat\alpha_i`$ from
falling outside $`(0, 1)`$ or $`\hat\beta_i`$ from being negative. When
theoretical constraints are desired, link functions can be applied:

``` math
\alpha_i = \mathrm{logit^{-1}}(\mathbf{X}_{\alpha,i}^\top \mathbf{a})
  = \frac{1}{1 + \exp(-\mathbf{X}_{\alpha,i}^\top \mathbf{a})},
  \qquad
\beta_i = \exp(\mathbf{X}_{\beta,i}^\top \mathbf{b}),
```

ensuring $`\alpha_i \in (0, 1)`$ and $`\beta_i > 0`$. Constrained
estimation is available for the Poisson and Negative Binomial methods
(set `constrained = TRUE`). Here the fitted values are
`alpha_values = plogis(X_alpha %*% alpha_coefs)` and
`beta_values = exp(X_beta %*% beta_coefs)`, while the stored coefficient
vectors remain on the transformed scales; response-scale values are
shown by [`summary()`](https://rdrr.io/r/base/summary.html).

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
#> Method: POISSON | estimator: MLE | link_rho: power | vcov: HC3 
#> N obs: 1382 
#> Gamma: 0.004998 (estimated) 
#> Log-likelihood: -11113.89 
#> AIC: 22237.78  BIC: 22263.94 
#> Deviance: 19624.98 
#> 
#> Coefficients:
#>                   Estimate Std. Error z value  Pr(>|z|)    
#> alpha:(Intercept) 0.670241   0.077201  8.6818 < 2.2e-16 ***
#> alpha:sexMale     0.075693   0.085822  0.8820 0.3777948    
#> beta:(Intercept)  0.528613   0.140139  3.7721 0.0001619 ***
#> beta:sexMale      0.082146   0.139501  0.5889 0.5559543    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>            Observed Estimate Estimate (BC) CI lower CI upper
#> sex=Female    4,512   55,915        55,839   12,446  250,521
#> sex=Male     22,593  198,635       198,550   78,554  501,847

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
#> Method: POISSON | estimator: MLE | link_rho: power | vcov: HC3 
#> N obs: 1382 
#> Gamma: 0.005433 (estimated) 
#> Log-likelihood: -11133.17 
#> AIC: 22274.33  BIC: 22295.26 
#> Deviance: 19663.53 
#> 
#> Coefficients (link scale: logit for alpha, log for beta):
#>                   Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)  0.82810    0.24579  3.3691 0.000754 ***
#> alpha:sexMale      0.21608    0.18498  1.1681 0.242759    
#> beta              -0.52088    0.17947 -2.9024 0.003704 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Response-scale fitted parameters (alpha_values in (0,1), beta_values > 0):
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
#> sex=Female    4,512   72,385        72,356   25,689  203,801
#> sex=Male     22,593  186,085       186,032   81,312  425,616
```

## Estimation Methods

The `uncounted` package implements four estimation methods. All share
the same mean structure; they differ in the objective function and
distributional assumptions.

### OLS on the log scale

Ordinary least squares applied to the log-linearized model:

``` math
\log m_i = \alpha \log N_i + \beta \log\!\left(\gamma + \frac{n_i}{N_i}\right) + \varepsilon_i.
```

When some $`m_i = 0`$, the response is replaced by $`\log(m_i + 1)`$.
OLS is fast and transparent but ignores the count nature of $`m_i`$ and
may be biased under heteroscedasticity (the log-transformation changes
the conditional mean).

When $`\gamma`$ is estimated, it is profiled out via grid search: for
each candidate $`\gamma`$, the OLS coefficients are obtained in closed
form, and the value minimizing the residual sum of squares is selected.

### NLS on the original scale

Nonlinear least squares minimizes

``` math
\sum_{i=1}^n \left(m_i - N_i^{\alpha}
  \left(\gamma + \frac{n_i}{N_i}\right)^{\!\beta}\right)^2
```

using the L-BFGS-B quasi-Newton algorithm. This avoids the
log-transformation bias of OLS but still treats $`m_i`$ as continuous.

### Poisson pseudo-maximum likelihood (PPML)

The Poisson log-likelihood is

``` math
\ell(\boldsymbol{\theta}) = \sum_{i=1}^{n}
  \bigl[m_i \log \mu_i - \mu_i\bigr],
```

where $`\mu_i = N_i^{\alpha_i}\,(\gamma + n_i/N_i)^{\beta_i}`$. As shown
by Santos Silva and Tenreyro (2006), the Poisson PML estimator is
consistent for the parameters of the conditional mean
$`E[m_i \mid N_i, n_i]`$ even if the true distribution is not Poisson,
provided the mean is correctly specified. This makes PPML robust to
heteroscedasticity and is the recommended default in `uncounted`.

The PPML point estimates coincide exactly with those of a Poisson GLM —
both solve the same score equations
$`\sum_i (m_i - \mu_i)\mathbf{z}_i = \mathbf{0}`$ — but by default
`uncounted` reports heteroskedasticity-consistent (HC3) standard errors
rather than the model-based Fisher-information variance. This robust
variance is what distinguishes the *pseudo*-ML estimator from naive
Poisson ML; it is also more general than quasi-Poisson, which assumes
$`\mathrm{Var}(m_i) = \phi\mu_i`$.

For the Poisson and Negative Binomial count models, the package can also
use `estimator = "gmm"` or `estimator = "el"`. In those cases the same
per-observation score equations are treated as moment conditions and
fitted with `momentfit`. The model for $`\mu_i`$ is unchanged; only the
estimating criterion differs. For these moment-based count fits the
default robust covariance is HC1 rather than HC3, because leverage-based
HC2/HC3 adjustments are not available generically in that setting.

### Negative Binomial ML

The Negative Binomial extends the Poisson by adding a dispersion
parameter $`\theta > 0`$. The log-likelihood is

``` math
\ell(\boldsymbol{\theta}, \theta) = \sum_{i=1}^{n} \Bigl[
  \log\Gamma(m_i + \theta) - \log\Gamma(\theta) - \log(m_i!)
  + \theta\log\theta - \theta\log(\theta + \mu_i)
  + m_i\log\mu_i - m_i\log(\theta + \mu_i)
\Bigr].
```

This accounts for overdispersion (variance exceeding the Poisson mean)
and can provide more efficient estimates when the data are truly
overdispersed. The parameter $`\theta`$ is jointly estimated alongside
$`\alpha`$, $`\beta`$, and $`\gamma`$.

### iOLS (Gamma PML)

The iOLS estimator (Benatia, Bellégo & Pape, 2025) targets the Gamma
Pseudo-Maximum Likelihood (GPML) score equations

``` math
\sum_i (m_i / \mu_i - 1)\,\mathbf{z}_i = \mathbf{0}
```

via a two-phase iterated OLS algorithm. Unlike PPML, which weights
observations by $`\mu_i`$, GPML gives equal weight to all observations.
This makes it more robust to influential high-count observations. See
[`vignette("iols")`](https://ncn-foreigners.github.io/uncounted/articles/iols.md)
for the full algorithm and comparison.

Currently `gamma = "estimate"` is not supported for iOLS; use a fixed
gamma from a Poisson fit.

### Gamma estimation

The offset parameter $`\gamma`$ is handled differently depending on the
method:

- **Poisson / NB:** $`\gamma`$ is jointly optimized within the
  likelihood using L-BFGS-B with box constraints (`gamma_bounds`,
  default $`[10^{-10},\, 0.5]`$).
- **OLS:** $`\gamma`$ is profiled out via a grid search over the bounds,
  followed by local optimization
  ([`optimize()`](https://rdrr.io/r/stats/optimize.html)).
- **All methods:** $`\gamma`$ can be fixed at a known value
  (`gamma = 0.01`) or omitted entirely (`gamma = NULL`), the latter
  requiring $`n_i > 0`$ for all observations.

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

fit_iols <- estimate_hidden_pop(
  data = irregular_migration,
  observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
  method = "iols", gamma = fit_pois$gamma
)

## Extract population size estimates
popsize(fit_ols)
#>   group observed estimate estimate_bc    lower    upper share_pct
#> 1 (all)    27105 25162.53    25035.32 19527.68 32096.36       100
popsize(fit_pois)
#>   group observed estimate estimate_bc    lower  upper share_pct
#> 1 (all)    27105 290597.4    290498.3 127726.4 660703       100
popsize(fit_nb)
#>   group observed estimate estimate_bc    lower   upper share_pct
#> 1 (all)    27105  1259529     1231371 731342.8 2073274       100
popsize(fit_iols)
#>   group observed estimate estimate_bc    lower    upper share_pct
#> 1 (all)    27105   163948    163244.4 76109.98 350134.4       100
```

## Inference

### Robust standard errors

All estimation methods report heteroscedasticity-consistent (HC) robust
standard errors computed via the sandwich estimator. The general form is

``` math
\widehat{\mathrm{Var}}(\hat{\boldsymbol{\theta}}) =
  (\mathbf{J}^\top \mathbf{J})^{-1}\,
  \mathbf{J}^\top \,\mathrm{diag}(\tilde{e}_1^2, \ldots, \tilde{e}_n^2)\,
  \mathbf{J}\,
  (\mathbf{J}^\top \mathbf{J})^{-1},
```

where $`\mathbf{J}`$ is the Jacobian (or score matrix) and
$`\tilde{e}_i`$ are adjusted residuals. The adjustment depends on the HC
variant:

| Variant | Residual adjustment                                               |
|---------|-------------------------------------------------------------------|
| HC0     | $`\tilde{e}_i = e_i`$ (no correction)                             |
| HC1     | $`\tilde{e}_i = e_i \sqrt{n/(n-p)}`$                              |
| HC2     | $`\tilde{e}_i = e_i / \sqrt{1 - h_{ii}}`$                         |
| HC3     | $`\tilde{e}_i = e_i / (1 - h_{ii})`$ (default)                    |
| HC4     | $`\tilde{e}_i = e_i / (1 - h_{ii})^{\delta_i/2}`$                 |
| HC5     | $`\tilde{e}_i = e_i / (1 - h_{ii})^{\min(\delta_i, h_{\max})/2}`$ |

Here $`h_{ii}`$ are the hat-matrix diagonal elements (leverages) and
$`\delta_i`$ is a data-dependent exponent. HC3 is the default for OLS,
NLS, Poisson MLE, NB MLE, and iOLS, and is recommended for moderate
sample sizes as it provides a good balance between finite-sample bias
correction and stability. For the moment-based count estimators
(`estimator = "gmm"` and `estimator = "el"`), the default robust
covariance is HC1.

### Delta-method confidence intervals for population size

The population size $`\hat\xi = \sum_i N_i^{\hat\alpha}`$ is a nonlinear
function of $`\hat\alpha`$. The implemented subgroup intervals are based
on the delta method for $`\hat\xi_g`$ itself, followed by a log-normal
approximation to preserve positivity.

For a subgroup $`g`$, let $`\mathbf{g}_g`$ denote the gradient of
$`\hat\xi_g = \sum_{i \in g} N_i^{\hat\alpha_i}`$ with respect to the
alpha coefficients. In the unconstrained case,

``` math
\mathbf{g}_g =
  \sum_{i \in g} N_i^{\hat\alpha_i}\log(N_i)\mathbf{x}_i.
```

When `constrained = TRUE`, the same expression is multiplied by the
derivative of the inverse-logit map, $`\hat\alpha_i(1-\hat\alpha_i)`$,
because the model is parameterized on the link scale.

Using the HC-robust covariance $`\mathbf{V}`$ for the alpha
coefficients,

``` math
\widehat{\mathrm{Var}}(\hat\xi_g) =
  \mathbf{g}_g^\top \mathbf{V} \mathbf{g}_g.
```

The subgroup standard error
$`\widehat{\mathrm{se}}(\hat\xi_g)=\sqrt{\widehat{\mathrm{Var}}(\hat\xi_g)}`$
is then mapped to a log-normal interval

``` math
\left(
\hat\xi_g \exp\!\left\{-z_{\alpha/2}
  \frac{\widehat{\mathrm{se}}(\hat\xi_g)}{\hat\xi_g}\right\},
\;
\hat\xi_g \exp\!\left\{z_{\alpha/2}
  \frac{\widehat{\mathrm{se}}(\hat\xi_g)}{\hat\xi_g}\right\}
\right).
```

When analytical bias correction is requested, the package rescales both
bounds by the ratio $`\hat\xi_g^{BC}/\hat\xi_g`$ so that the returned
interval is centered on the bias-corrected estimate.

**Total across groups.** The package does not report totals by default.
Totals are optional and should be requested only when the groups are
non-overlapping subpopulations. They are usually not meaningful for
panel or time-period groupings because the same individuals may
contribute to multiple periods. When a total is requested,
$`\hat\xi = \sum_g \hat\xi_g`$ has its own delta-method CI computed via
the gradient $`\nabla_\alpha \xi`$ and a log-normal approximation to
ensure positivity.

### Bias correction

Because $`\xi(\alpha) = \sum_i N_i^\alpha`$ is convex in $`\alpha`$ for
$`N_i > 1`$, Jensen’s inequality implies $`E[\hat\xi] \geq \xi`$.
Writing $`N_i^{\hat\alpha_i} = \exp(\hat\alpha_i \log N_i)`$ and
treating $`\hat{\boldsymbol{\alpha}}`$ as approximately normal gives

``` math
E[N_i^{\hat\alpha_i}] \approx N_i^{\alpha_i}
  \exp\!\left\{\frac{1}{2}(\log N_i)^2
  \mathbf{x}_i^\top \mathbf{V}_{\mathrm{model}} \mathbf{x}_i\right\}.
```

Inverting this yields the multiplicative lognormal correction used by
the package for unconstrained models:

``` math
\hat\xi^{BC} =
  \sum_{i=1}^{n}
  N_i^{\hat\alpha_i}
  \exp\!\left(
    -\frac{1}{2}(\log N_i)^2
    \mathbf{x}_i^\top \mathbf{V}_{\mathrm{model}} \mathbf{x}_i
  \right).
```

This correction is strictly positive and exact under the Gaussian
approximation for the linear predictor. For constrained fits, where
$`\alpha_i = \operatorname{logit}^{-1}(\eta_i)`$, the corresponding
logistic-normal integral has no closed form, so the package falls back
to a second-order Taylor approximation on the link scale:

``` math
\mathrm{Bias}(\hat\xi) \;\approx\;
  \frac{1}{2}\,\sum_{i=1}^{n} N_i^{\alpha}\,(\log N_i)^2 \;\cdot\;
  \mathbf{x}^\top \mathbf{V}_{\mathrm{model}}\,\mathbf{x},
```

where $`\mathbf{V}_{\mathrm{model}}`$ is the model-based (homoscedastic)
variance–covariance matrix.

Model-based variance is used rather than the HC-robust variance, because
the latter can be inflated by high-leverage observations in skewed data,
leading to overcorrection. For Poisson and NB fits, `vcov_model` is the
Fisher-style inverse information for the mean-model parameters. This
remains true when the coefficients are obtained by `estimator = "gmm"`
or `estimator = "el"`: the package evaluates the same Fisher-style model
variance at the non-MLE estimate and reserves the robust HC or FWB
covariance for interval estimation. Bias correction is also applied to
the CI bounds, evaluated at their respective $`\alpha`$ values.

### Fractional weighted bootstrap

For finite-sample inference and as a robustness check on the analytical
CIs, the package provides the fractional weighted bootstrap (FWB) of Xu
et al. (2020) via
[`bootstrap_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/bootstrap_popsize.md),
implemented with the `fwb` package (Greifer, 2026).

Instead of resampling rows, FWB generates random non-negative weights
$`w_1, \ldots, w_n`$ from a $`\mathrm{Dirichlet}(1, \ldots, 1)`$
distribution and refits the model under these weights. This avoids
duplicate-row issues with discrete data and is naturally suited to
weighted GLMs.

**Cluster bootstrap.** When observations are not independent (e.g., the
same country observed across years and sexes), a cluster variable can be
specified. In this case, all observations within the same cluster
receive the same FWB weight, preserving within-cluster correlation.

**Bootstrap CIs.** Two types are available:

- *Percentile*: the $`\alpha/2`$ and $`(1 - \alpha/2)`$ quantiles of the
  bootstrap distribution.
- *Bias-corrected percentile (BC)*: adjusts the quantile probabilities
  for median bias using the standard BC formula.

The bootstrap median is the default point estimate because the
transformation $`\xi = \sum N_i^\alpha`$ is convex, so the bootstrap
distribution of $`\hat\xi`$ is right-skewed and the bootstrap mean tends
to overestimate $`\xi`$.

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
#> (all) 290,597     290,498     250,522   330,188  178,051  845,465

## Access the full results table
boot_result$popsize_full
#>   group   plugin plugin_bc boot_median boot_mean  lower    upper
#> 1 (all) 290597.4  290498.3    250522.1  330188.4 178051 845464.7
```

### Dependence diagnostics

The decomposition
``` math
\mu_i = \xi_i \rho_i
```
is an identifying assumption. Aggregated observed counts do not test
that separability directly, so the package treats sensitivity analysis
as a separate diagnostic layer rather than as part of the ordinary
confidence interval.

Two complementary tools are available:

- [`dependence_bounds()`](https://ncn-foreigners.github.io/uncounted/reference/dependence_bounds.md)
  keeps the baseline estimate fixed and reports an identification
  envelope under
  ``` math
  \mu_i = \xi_i \rho_i \kappa_i, \qquad 1 / \Gamma \le \kappa_i \le \Gamma.
  ```
- [`profile_dependence()`](https://ncn-foreigners.github.io/uncounted/reference/profile_dependence.md)
  imposes a parametric offset
  ``` math
  \mu_i(\delta) = \exp(\delta)\,\xi_i \rho_i
  ```
  and refits the model, so the plug-in estimate $`\hat\xi(\delta)`$ can
  move with the sensitivity parameter.

``` r

bounds <- dependence_bounds(fit, Gamma = c(1, 1.1, 1.25, 1.5))
prof_dep <- profile_dependence(
  fit,
  delta_grid = seq(-0.5, 0.5, length.out = 9),
  plot = FALSE
)
robustness_dependence(prof_dep, threshold = 500000)
```

For threshold questions on the bootstrap distribution, use
[`exceedance_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/exceedance_popsize.md)
on the result of
[`bootstrap_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/bootstrap_popsize.md)
to obtain an empirical bootstrap tail area.

A third diagnostic,
[`frailty_sensitivity()`](https://ncn-foreigners.github.io/uncounted/reference/frailty_sensitivity.md),
targets the identifying assumption more directly. It treats residual
dependence between hidden stock and detection as an omitted shared
frailty, linearizes the fitted Poisson/NB mean model with its IRLS
working response, and indexes that omitted factor by weighted
partial-$`R^2`$ parameters. The resulting first-order bias is then
mapped into the grouped hidden-population estimand $`\Xi`$, so the
output is on the same scale as `popsize(by = ...)` rather than on the
raw coefficient scale.

``` r

frailty_sensitivity(
  fit,
  by = ~ year,
  r2_d = c(0, 0.05, 0.10),
  r2_y = c(0, 0.05, 0.10),
  plot = FALSE
)
```

The full derivation, including the omitted-frailty model, IRLS
linearization, weighted bias formulas, and robustness-value definitions,
is documented in
[`?frailty_sensitivity`](https://ncn-foreigners.github.io/uncounted/reference/frailty_sensitivity.md).

## References

- Beręsewicz, M., & Pawlukiewicz, K. (2020). Estimation of the number of
  irregular foreigners in Poland using non-linear count regression
  models. arXiv preprint arXiv:2008.09407.
  <https://arxiv.org/abs/2008.09407>

- Greifer, N. (2026). fwb: An R package for the fractional weighted
  bootstrap. *Journal of Open Source Software*, 11(124), 8880.
  <https://doi.org/10.21105/joss.08880>

- Lam, H. (2022). A cheap bootstrap method for fast inference. arXiv
  preprint. <https://doi.org/10.48550/arXiv.2202.00090>

- Santos Silva, J. M. C. and Tenreyro, S. (2006). The log of gravity.
  *The Review of Economics and Statistics*, 88(4), 641–658.
  <https://doi.org/10.1162/rest.88.4.641>

- Xu, L., Gotwalt, C., Hong, Y., King, C. B., and Meeker, W. Q. (2020).
  Applications of the fractional-random-weight bootstrap. *The American
  Statistician*, 74(4), 345–358.
  <https://doi.org/10.1080/00031305.2020.1731599>

- Zhang, L.-C. (2008). *Developing methods for determining the number of
  unauthorized foreigners in Norway* (Documents 2008/11). Statistics
  Norway.
  <https://www.ssb.no/a/english/publikasjoner/pdf/doc_200811_en/doc_200811_en.pdf>
