# Estimate the Size of an Unauthorized Migrant Population

Fits a power-law model relating observed counts of unauthorized migrants
to reference population sizes and auxiliary registration counts. The
model is estimated using OLS, NLS, Poisson MLE, or Negative Binomial
MLE, and supports covariate-varying parameters through interaction
formulas.

## Usage

``` r
estimate_hidden_pop(
  data,
  observed,
  auxiliary,
  reference_pop,
  method = c("poisson", "nb", "ols", "nls", "iols"),
  cov_alpha = NULL,
  cov_beta = NULL,
  gamma = "estimate",
  cov_gamma = NULL,
  gamma_bounds = c(1e-10, 0.5),
  theta_start = 1,
  link_rho = c("power", "cloglog", "logit", "probit"),
  estimator = c("mle", "gmm", "el"),
  vcov = "HC3",
  weights = NULL,
  constrained = FALSE,
  countries = NULL,
  cluster = NULL
)
```

## Arguments

- data:

  A data frame containing all variables.

- observed:

  One-sided formula for the observed count (e.g., `~ m`).

- auxiliary:

  One-sided formula for the auxiliary count (e.g., `~ n`).

- reference_pop:

  One-sided formula for the reference population (e.g., `~ N`).

- method:

  Estimation method: `"poisson"` (default), `"nb"`, `"ols"`, or `"nls"`.
  See Details.

- cov_alpha:

  Formula for covariates in alpha (e.g., `~ sex + year`). Default `NULL`
  means a single alpha (intercept only).

- cov_beta:

  Formula for covariates in beta (e.g., `~ sex`). Default `NULL` means a
  single beta (intercept only).

- gamma:

  Controls the gamma offset in the rate term:

  - `"estimate"` (default): gamma is estimated from data.

  - A numeric value: fixed gamma, uses \\\log(\gamma + n_i / N_i)\\.

  - `NULL`: no gamma, uses \\\log(n_i / N_i)\\ (requires \\n_i \> 0\\).

- cov_gamma:

  Optional formula for covariates in gamma (e.g., `~ sex` or `~ year`).
  When specified, gamma varies across observations via \\\gamma_i =
  \exp(\mathbf{X}\_{\gamma,i}' \boldsymbol{\delta})\\, where
  \\\boldsymbol{\delta}\\ are estimated coefficients on the log scale.
  Requires `gamma = "estimate"` and `method` `"poisson"` or `"nb"`.
  Default `NULL` means a single scalar gamma (intercept only).

- gamma_bounds:

  Numeric vector of length 2: lower and upper bounds for the estimated
  gamma parameter. Default `c(1e-10, 0.5)`.

- theta_start:

  Starting value for the NB dispersion parameter (used only when
  `method = "nb"`). Default 1.

- link_rho:

  Link for the detection component \\\rho_i = h(\beta_i \log(\gamma_i +
  n_i / N_i))\\. Supported values are `"power"` (default, \\\rho_i =
  \exp(\eta_i)\\), `"cloglog"` (\\\rho_i = 1 - \exp(-\exp(\eta_i))\\),
  and `"logit"` (\\\rho_i = 1 / (1 + \exp(-\eta_i))\\), and `"probit"`
  (\\\rho_i = \Phi(\eta_i)\\). The legacy alias `"logistic"` is accepted
  and normalized to `"logit"`. Bounded links are currently available for
  `method = "poisson"`, `"nb"`, and `"nls"`.

- estimator:

  Estimation criterion for count models: `"mle"` (default), `"gmm"`, or
  `"el"`. Non-MLE estimators are currently available only for
  `method = "poisson"` and `method = "nb"`.

- vcov:

  Controls the variance-covariance estimator. Can be:

  - A character string specifying the HC type: `"HC0"` through `"HC5"`,
    or `"HC4m"`. Default `"HC3"`. When `cluster` is also provided, this
    type is passed to
    [`sandwich::vcovCL()`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html).
    For count models fitted with `estimator = "gmm"` or
    `estimator = "el"`, the default robust covariance is `"HC1"` and
    HC2+ requests are downgraded to `"HC1"`.

  - A function that takes the fitted `uncounted` object and returns a
    variance-covariance matrix. For example,
    [`sandwich::vcovHC`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html)
    or `function(x) sandwich::vcovCL(x, cluster = data$country_code)`.

- weights:

  Optional numeric vector of observation weights.

- constrained:

  Logical. If `TRUE` and `method %in% c("poisson", "nb")`, the fitted
  per-observation parameters are constrained on the response scale via
  \\\alpha_i = \mathrm{logit}^{-1}(\mathbf{X}\_{\alpha,i}' \mathbf{a})\\
  and \\\beta_i = \exp(\mathbf{X}\_{\beta,i}' \mathbf{b})\\. This keeps
  \\\alpha_i \in (0, 1)\\ and \\\beta_i \> 0\\. The reported coefficient
  vectors `alpha_coefs` and `beta_coefs` remain on the transformed
  scales (logit for alpha, log for beta); use
  [`summary()`](https://rdrr.io/r/base/summary.html) or `alpha_values` /
  `beta_values` for response-scale values. Default `FALSE`. See Details.

- countries:

  One-sided formula identifying a country or group variable used when
  reporting population size estimates (e.g., `~ country`). This does
  **not** enable cluster-robust variance; use `cluster` for that.

- cluster:

  Optional one-sided formula identifying a cluster variable for
  cluster-robust variance estimation (e.g., `~ country_code`). When
  provided and `vcov` is a character string, the variance is computed
  using
  [`sandwich::vcovCL()`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html)
  with the specified HC type. Note: clustered HC2 and HC3 are only
  applicable to standard linear and generalized linear models;
  [`sandwich::vcovCL()`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html)
  will emit a warning for non-GLM objects. Use `"HC0"` or `"HC1"` with
  clustering to avoid this. Ignored when `vcov` is a function.

## Value

An object of class `"uncounted"`, a list containing:

- `coefficients`:

  Named vector of all estimated coefficients.

- `alpha_coefs`, `beta_coefs`:

  Coefficient sub-vectors for the alpha and beta equations.

- `alpha_values`, `beta_values`:

  Per-observation alpha and beta on the response scale (after applying
  link inverse when constrained).

- `vcov`:

  HC-robust variance-covariance matrix.

- `vcov_model`:

  Model-based (homoscedastic) variance-covariance matrix, used for bias
  correction. For Poisson and NB count models, this is the Fisher-style
  inverse information for the mean-model parameters, evaluated at the
  fitted coefficients even when `estimator = "gmm"` or
  `estimator = "el"`.

- `rho_values`:

  Per-observation fitted detection component \\\hat\rho_i\\.

- `fitted.values`:

  Fitted values \\\hat{m}\_i\\.

- `residuals`:

  Raw residuals \\m_i - \hat{m}\_i\\.

- `gamma`:

  Estimated or fixed gamma value (or `NULL`).

- `theta`:

  NB dispersion parameter (only for `method = "nb"`).

- `loglik`:

  Log-likelihood (for MLE count fits only).

- `method`:

  The model family used.

- `estimator`:

  The estimation criterion used.

- `link_rho`:

  The detection-link transformation used.

## Details

**Model specification.** For observation \\i\\, the expected observed
count \\m_i\\ is modelled as:

\$\$E(m_i) = \xi_i \rho_i\$\$

where \\\xi_i = E(M_i \mid N_i)\\ is the theoretical unauthorized
population size and \\\rho_i = E(p_i \mid N_i, n_i)\\ is the theoretical
detection rate. In the baseline empirical specification, \\\xi_i =
N_i^{\alpha_i}\\ and, under the default `link_rho = "power"`, \\\rho_i =
(\gamma_i + n_i / N_i)^{\beta_i}\\. More generally, the package writes
the detection component as \$\$\eta_i = \beta_i \log(\gamma_i + n_i /
N_i), \qquad \rho_i = h(\eta_i)\$\$ where \\h\\ is chosen by `link_rho`.
Supported links are:

- `"power"`:

  \\\rho_i = \exp(\eta_i)\\

- `"cloglog"`:

  \\\rho_i = 1 - \exp(-\exp(\eta_i))\\

- `"logit"`:

  \\\rho_i = 1 / (1 + \exp(-\eta_i))\\

- `"probit"`:

  \\\rho_i = \Phi(\eta_i)\\

The `"power"` link reproduces the paper's baseline power-law detection
model exactly; the bounded links are package extensions. On the log
scale,

\$\$\log E(m_i) = \alpha_i \log N_i + \log(\rho_i)\$\$

so for every supported `link_rho`, \\\xi_i = E(m_i) / \rho_i =
N_i^{\alpha_i}\\.

**Parameter interpretation.**

- \\\alpha\\:

  Elasticity of the observed count with respect to the reference
  population. When \\\alpha \< 1\\, the total (unobserved) population
  \\\xi = \sum_i N_i^\alpha\\ is smaller than \\\sum_i N_i\\, reflecting
  incomplete coverage. Values of \\\alpha\\ near 0 imply weak dependence
  on population size; values near 1 imply near-proportional scaling.

- \\\beta\\:

  Elasticity with respect to the registration rate \\\gamma + n_i /
  N_i\\. A positive \\\beta\\ means higher auxiliary rates are
  associated with more observed unauthorized migrants.

- \\\gamma\\:

  Baseline registration-rate offset, ensuring that the rate term is
  positive even when \\n_i = 0\\. Typically small (close to zero). When
  `gamma = "estimate"`, it is profiled out (OLS) or jointly optimized
  (MLE methods).

**Estimation methods.**

- OLS:

  Ordinary least squares on the log-linearized model. Fast and
  transparent but ignores the count nature of \\m_i\\ and may be
  inefficient under heteroscedasticity. Uses \\\log(m_i)\\ as the
  response (or \\\log(m_i + 1)\\ when zeros are present). Note: when the
  \\\log(m_i + 1)\\ transformation is used, fitted values are
  \\\exp(\hat{\mu})\\ where \\\hat{\mu}\\ was estimated on the shifted
  scale, so response-scale residuals and diagnostics are approximate.

- NLS:

  Nonlinear least squares on the original scale \\m_i = N_i^{\alpha}
  (\gamma + n_i/N_i)^{\beta} + \varepsilon_i\\. Avoids the
  log-transformation bias of OLS but still treats \\m_i\\ as continuous.

- Poisson:

  Poisson pseudo-maximum likelihood (PPML). Consistent under
  heteroscedasticity as long as the conditional mean is correctly
  specified (Santos Silva and Tenreyro, 2006). Recommended as the
  default.

- NB:

  Negative Binomial MLE. Adds a dispersion parameter \\\theta\\ to
  accommodate overdispersion beyond what the Poisson allows. Standard
  errors for the regression coefficients are computed conditional on
  \\\hat{\theta}\\, consistent with the approach used in
  [`MASS::glm.nb`](https://rdrr.io/pkg/MASS/man/glm.nb.html). This means
  coefficient SEs do not account for uncertainty in \\\theta\\
  estimation. For `method = "poisson"` and `method = "nb"`, the same
  moment conditions can also be fit with `estimator = "gmm"` or
  `estimator = "el"` using momentfit.

**Constrained vs. unconstrained estimation.** When `constrained = FALSE`
(default), \\\alpha\\ and \\\beta\\ are estimated on the real line
without restrictions. When `constrained = TRUE`, the linear predictors
are \\\eta\_{\alpha,i} = \mathbf{X}\_{\alpha,i}' \mathbf{a}\\ and
\\\eta\_{\beta,i} = \mathbf{X}\_{\beta,i}' \mathbf{b}\\, with
response-scale parameters \\\alpha_i =
\mathrm{logit}^{-1}(\eta\_{\alpha,i})\\ and \\\beta_i =
\exp(\eta\_{\beta,i})\\. This constrains the fitted \\\alpha_i\\ values
to \\(0, 1)\\ and the fitted \\\beta_i\\ values to be strictly positive.
The stored coefficient vectors `alpha_coefs` and `beta_coefs` remain on
the transformed scales (logit for alpha, log for beta); use
[`summary()`](https://rdrr.io/r/base/summary.html) or the `alpha_values`
/ `beta_values` elements of the returned object for response-scale
values. Constrained estimation is available for the Poisson and NB
methods only.

**Population size estimation.** The total unauthorized population is
estimated as \\\hat{\xi} = \sum_i N_i^{\hat{\alpha}\_i}\\. See
[`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
for bias correction and confidence intervals.

## References

Zhang, L.-C. (2008). Developing methods for determining the number of
unauthorized foreigners in Norway. *Documents* 2008/11, Statistics
Norway.
<https://www.ssb.no/a/english/publikasjoner/pdf/doc_200811_en/doc_200811_en.pdf>

Beręsewicz, M. and Pawlukiewicz, K. (2020). Estimation of the number of
irregular foreigners in Poland using non-linear count regression models.
*arXiv preprint* arXiv:2008.09407.

Santos Silva, J. M. C. and Tenreyro, S. (2006). The log of gravity. *The
Review of Economics and Statistics*, 88(4), 641–658.

## Examples

``` r
# Simulate data: 50 groups with known population structure
set.seed(42)
sim_data <- data.frame(
  N = sample(1000:50000, 50, replace = TRUE)
)
sim_data$n <- rpois(50, lambda = sim_data$N * 0.05)
alpha_true <- 0.6
beta_true <- 1.2
gamma_true <- 0.01
mu <- sim_data$N^alpha_true * (gamma_true + sim_data$n / sim_data$N)^beta_true
sim_data$m <- rpois(50, lambda = mu)

# Poisson with estimated gamma (default)
fit_pois <- estimate_hidden_pop(
  data = sim_data, observed = ~ m,
  auxiliary = ~ n, reference_pop = ~ N
)
summary(fit_pois)
#> Unauthorized population estimation
#> Method: POISSON | estimator: MLE | link_rho: power | vcov: HC3 
#> N obs: 50 
#> Gamma: 0.337767 (estimated) 
#> Log-likelihood: -134.08 
#> AIC: 274.16  BIC: 279.89 
#> Deviance: 45.82 
#> 
#> Coefficients:
#>         Estimate Std. Error z value Pr(>|z|)    
#> alpha   0.559800   0.066707  8.3919   <2e-16 ***
#> beta    3.113112 221.225021  0.0141   0.9888    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>       Observed Estimate Estimate (BC) CI lower CI upper
#> (all)      777   14,853        12,134    3,155   46,669

# OLS with fixed gamma
fit_ols <- estimate_hidden_pop(
  data = sim_data, observed = ~ m,
  auxiliary = ~ n, reference_pop = ~ N,
  method = "ols", gamma = 0.01
)
summary(fit_ols)
#> Unauthorized population estimation
#> Method: OLS | estimator: MLE | link_rho: power | vcov: HC3 
#> N obs: 50 
#> Gamma: 0.01 (fixed) 
#> Deviance: 4.13 
#> 
#> Coefficients:
#>       Estimate Std. Error t value  Pr(>|t|)    
#> alpha 0.647108   0.082283  7.8644 3.532e-10 ***
#> beta  1.374506   0.300282  4.5774 3.344e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>       Observed Estimate Estimate (BC) CI lower CI upper
#> (all)      777   36,563        33,187    6,271  175,626

# Constrained Poisson (alpha in (0,1), beta > 0)
fit_constr <- estimate_hidden_pop(
  data = sim_data, observed = ~ m,
  auxiliary = ~ n, reference_pop = ~ N,
  constrained = TRUE
)
summary(fit_constr)
#> Unauthorized population estimation
#> Method: POISSON | estimator: MLE | link_rho: power | vcov: HC3 
#> N obs: 50 
#> Gamma: 0 (estimated) 
#> Log-likelihood: -133.7 
#> AIC: 273.39  BIC: 279.13 
#> Deviance: 45.06 
#> 
#> Coefficients (link scale: logit for alpha, log for beta):
#>        Estimate Std. Error z value Pr(>|z|)
#> alpha  0.229427   0.263416  0.8710   0.3838
#> beta  -0.025829   0.230379 -0.1121   0.9107
#> 
#> Response-scale fitted parameters (alpha_values in (0,1), beta_values > 0):
#>   Alpha (response scale):
#>        alpha SE(alpha)
#> (all) 0.5571     0.065
#>   Beta (response scale):
#>     beta SE(beta)
#> 1 0.9745   0.2245
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>       Observed Estimate Estimate (BC) CI lower CI upper
#> (all)      777   14,447        11,858    3,192   44,055
```
