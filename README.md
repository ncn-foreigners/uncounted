
# The {uncounted} package <img src="man/figures/logo.png" align="right" width="150" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/ncn-foreigners/uncounted/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ncn-foreigners/uncounted/actions/workflows/R-CMD-check.yaml)
[![CRAN-check](https://github.com/ncn-foreigners/uncounted/actions/workflows/cran-check.yaml/badge.svg)](https://github.com/ncn-foreigners/uncounted/actions/workflows/cran-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/ncn-foreigners/uncounted/graph/badge.svg)](https://app.codecov.io/gh/ncn-foreigners/uncounted)
<!-- badges: end -->

R package for estimating the size of unauthorized migrant populations
using a power-law model that relates observed counts to reference
populations and auxiliary detection data. Supports OLS, NLS, Poisson,
Negative Binomial, and iOLS estimation, alternative bounded detection
links, optional gamma offsets, moment-based count estimators via the
`momentfit` package, formal population-size hypothesis tests, and
Bayesian Poisson/NB workflows through `brms` and Stan.

## Model

The paper’s theoretical model factors the expected observed count into a
latent population component and a detection component:

$$
\mu_i = E(m_i \mid N_i, n_i) = \xi_i \rho_i,
$$

where $\xi_i = E(M_i \mid N_i)$ is the theoretical unauthorized
population size and $\rho_i = E(p_i \mid N_i, n_i)$ is the theoretical
detection rate.

In the baseline empirical specification implemented by default,

$$
\xi_i = N_i^{\alpha_i}, \qquad
\rho_i = \left(\gamma_i + \frac{n_i}{N_i}\right)^{\beta_i},
$$

where $N_i$ is the reference (total registered) population, $n_i$ is an
auxiliary count (e.g. police records), and $\gamma_i \geq 0$ is a
baseline detection offset. This gives

$$
\mu_i = N_i^{\alpha_i}\left(\gamma_i + \frac{n_i}{N_i}\right)^{\beta_i}.
$$

The package also supports bounded alternatives for the detection
component by writing

$$
\eta_i = \beta_i \log\left(\gamma_i + \frac{n_i}{N_i}\right), \qquad
\rho_i = h(\eta_i),
$$

with `link_rho = "power"` (the paper’s baseline specification),
`"cloglog"`, `"logit"`, or `"probit"`. On the log scale, the mean
structure remains linear in the power-link case and is handled directly
for the nonlinear estimators.

## Installation

Install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("ncn-foreigners/uncounted")
```

## Quick start

``` r
library(uncounted)
data(irregular_migration)
irregular_migration$year <- as.factor(irregular_migration$year)
irregular_migration$ukr <- as.integer(irregular_migration$country_code == "UKR")
## Fit Poisson PML with year x UKR interaction in alpha
fit <- estimate_hidden_pop(
  data = irregular_migration,
  observed = ~ m,
  auxiliary = ~ n,
  reference_pop = ~ N,
  method = "poisson",
  cov_alpha = ~ year * ukr + sex,
  cov_beta = ~ year,
  countries = ~ country_code
)

summary(fit)
#> Unauthorized population estimation
#> Method: POISSON | vcov: HC3 
#> N obs: 1382 
#> Gamma: 0.006954 (estimated) 
#> Log-likelihood: -5738.09 
#> AIC: 11516.19  BIC: 11620.81 
#> Deviance: 8873.39 
#> 
#> Coefficients:
#>                       Estimate  Std. Error z value  Pr(>|z|)    
#> alpha:(Intercept)   8.2203e-01  1.4949e-01  5.4989 3.822e-08 ***
#> alpha:year2020     -1.4497e-02  1.7895e-01 -0.0810   0.93543    
#> alpha:year2021     -1.7807e-02  1.5228e-01 -0.1169   0.90691    
#> alpha:year2022     -1.4721e-02  1.4893e-01 -0.0988   0.92126    
#> alpha:year2023      3.3403e-02  1.5129e-01  0.2208   0.82525    
#> alpha:year2024      1.1142e-02  1.6500e-01  0.0675   0.94616    
#> alpha:ukr           7.5122e-03  7.5863e-02  0.0990   0.92112    
#> alpha:sexMale       1.0979e-02  3.0053e-02  0.3653   0.71487    
#> alpha:year2020:ukr  9.0019e-03  8.7711e-02  0.1026   0.91826    
#> alpha:year2021:ukr  6.0173e-05  7.7873e-02  0.0008   0.99938    
#> alpha:year2022:ukr -7.6404e-02  8.5391e-02 -0.8947   0.37092    
#> alpha:year2023:ukr -1.5379e-01  8.7395e-02 -1.7597   0.07846 .  
#> alpha:year2024:ukr -1.7573e-01  1.0221e-01 -1.7192   0.08558 .  
#> beta:(Intercept)    7.1261e-01  2.9983e-01  2.3767   0.01747 *  
#> beta:year2020       1.5887e-01  3.5610e-01  0.4461   0.65550    
#> beta:year2021       2.2779e-01  3.0389e-01  0.7496   0.45350    
#> beta:year2022       2.0279e-01  2.9758e-01  0.6815   0.49558    
#> beta:year2023       1.8773e-01  2.9954e-01  0.6267   0.53083    
#> beta:year2024      -1.1125e-02  3.2368e-01 -0.0344   0.97258    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>                              Observed Estimate Estimate (BC) CI lower CI upper
#> year=2019, ukr=0, sex=Female      426    9,754         9,735    1,137   83,323
#> year=2019, ukr=1, sex=Female    1,109   21,953        21,893    1,723  278,181
#> year=2019, ukr=0, sex=Male      2,177   24,431        24,382    2,297  258,819
#> year=2019, ukr=1, sex=Male      2,892   40,670        40,572    2,545  646,714
#> year=2020, ukr=0, sex=Female      193    9,703         9,671    1,821   51,370
#> year=2020, ukr=1, sex=Female      505   22,809        22,702    2,857  180,413
#> year=2020, ukr=0, sex=Male      1,158   24,242        24,149    4,242  137,490
#> year=2020, ukr=1, sex=Male      1,542   41,068        40,875    5,500  303,783
#> year=2021, ukr=0, sex=Female      144   11,648        11,614    4,350   31,008
#> year=2021, ukr=1, sex=Female      339   22,486        22,387    5,936   84,433
#> year=2021, ukr=0, sex=Male      1,337   29,892        29,789   12,296   72,168
#> year=2021, ukr=1, sex=Male      1,285   40,132        39,953   13,980  114,183
#> year=2022, ukr=0, sex=Female      211   14,366        14,339    7,110   28,918
#> year=2022, ukr=1, sex=Female      106   14,062        14,011    2,983   65,817
#> year=2022, ukr=0, sex=Male      1,903   38,823        38,735   21,717   69,088
#> year=2022, ukr=1, sex=Male        729   16,379        16,323    3,975   67,040
#> year=2023, ukr=0, sex=Female      444   24,162        24,128    9,999   58,218
#> year=2023, ukr=1, sex=Female       79    9,632         9,605    2,051   44,990
#> year=2023, ukr=0, sex=Male      3,283   70,186        70,070   30,299  162,047
#> year=2023, ukr=1, sex=Male        556   11,049        11,020    2,642   45,963
#> year=2024, ukr=0, sex=Female      887   22,442        22,416    5,744   87,473
#> year=2024, ukr=1, sex=Female       69    5,085         5,075      641   40,157
#> year=2024, ukr=0, sex=Male      5,041   62,216        62,137   14,156  272,745
#> year=2024, ukr=1, sex=Male        690    6,192         6,179      759   50,296
```

### Population size estimates

``` r
## By year
popsize(fit, by = ~ year)
#>   group observed  estimate estimate_bc     lower     upper share_pct
#> 1  2019     6604  96808.06    96581.24  8228.355 1133633.1  16.31464
#> 2  2020     3398  97821.08    97397.04 15271.644  621163.2  16.48536
#> 3  2021     3105 104157.85   103743.46 38794.612  277427.9  17.55327
#> 4  2022     2949  83630.32    83407.98 37693.079  184566.8  14.09385
#> 5  2023     4362 115029.17   114822.64 49241.000  267749.2  19.38537
#> 6  2024     6687  95935.01    95806.63 23177.438  396027.8  16.16751
```

### Hypothesis tests for population size

`hypotheses_popsize()` tests threshold and contrast statements about
grouped population-size estimates. The expression uses explicit
`xi[...]` filters, so the target is unambiguous.

``` r
ps_year <- popsize(fit, by = ~ year)

## H0: xi[2024] <= 200,000
## H1: xi[2024] > 200,000
hypotheses_popsize(
  ps_year,
  "xi[year == 2024] > 200000",
  estimate = "estimate_bc"
)

## H0: xi[2024] - xi[2019] <= 0
## H1: xi[2024] - xi[2019] > 0
hypotheses_popsize(
  ps_year,
  "xi[year == 2024] - xi[year == 2019] > 0",
  estimate = "estimate_bc"
)
```

The same syntax works with fractional weighted bootstrap results:

``` r
boot_year <- bootstrap_popsize(
  fit,
  R = 499,
  cluster = ~ country_code,
  by = ~ year,
  verbose = FALSE
)

hypotheses_popsize(
  boot_year,
  "xi[year == 2024] > 200000",
  estimate = "boot_median"
)
```

For frequentist `popsize()` and `bootstrap_popsize()` objects, the
helper reports Wald-style statistics and p values. Use
`exceedance_popsize()` when you want an empirical bootstrap tail area
instead.

### Bootstrap confidence intervals

``` r
## Cluster FWB bootstrap (by country)
boot <- bootstrap_popsize(fit, R = 499,
                          cluster = ~ country_code,
                          by = ~ year)
boot
```

### Bayesian brms/Stan workflow

Bayesian models use `brms` as the modelling interface and Stan as the
sampler backend. The default backend is `rstan`; `cmdstanr` can also be
used if it is installed.

``` r
fit_bayes <- estimate_hidden_pop_bayes(
  data = irregular_migration,
  observed = ~ m,
  auxiliary = ~ n,
  reference_pop = ~ N,
  method = "poisson",
  cov_alpha = ~ year * ukr + sex,
  cov_beta = ~ year,
  gamma = 0.005,
  backend = "rstan",
  chains = 4,
  cores = min(4, parallel::detectCores()),
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed = 20260425
)

summary(fit_bayes, total = TRUE)

ps_bayes_year <- popsize(fit_bayes, by = ~ year, total = TRUE)

hypotheses_popsize(
  ps_bayes_year,
  "xi[year == 2024] > 200000"
)
```

For Bayesian `popsize()` objects, directional hypotheses are posterior
propositions. The output reports quantities such as `Pr(H1 | data)`,
`Pr(H0 | data)`, posterior odds, and credible intervals, not frequentist
p values.

The underlying `brmsfit` is available for diagnostics:

``` r
bfit <- as_brmsfit(fit_bayes)

plot(bfit)
brms::pp_check(bfit, type = "rootogram", ndraws = 100)
brms::loo(bfit)
```

### Dependence and threshold diagnostics

``` r
## Identification bounds for residual dependence between xi and detection
bounds <- dependence_bounds(fit, Gamma = c(1, 1.1, 1.25, 1.5))
bounds

## Moving-xi profile under a parametric dependence offset
prof_dep <- profile_dependence(
  fit,
  delta_grid = seq(-0.5, 0.5, length.out = 9),
  plot = FALSE
)
head(prof_dep)

## Smallest dependence strength needed to cross a target
robustness_dependence(prof_dep, threshold = 500000)

## Bootstrap tail area for threshold questions
boot_total <- bootstrap_popsize(
  fit,
  R = 199,
  cluster = ~ country_code,
  total = TRUE,
  verbose = FALSE
)
exceedance_popsize(boot_total, threshold = 500000)

## Omitted-frailty sensitivity for grouped Xi
frailty_sensitivity(
  fit,
  by = ~ year,
  r2_d = c(0, 0.05, 0.10),
  r2_y = c(0, 0.05, 0.10),
  plot = FALSE
)
```

### Diagnostics

``` r
## Leave-one-out influence
db <- dfbeta(fit)
dp <- dfpopsize(fit)

## Which countries have the largest influence on population size?
head(sort(abs(dp), decreasing = TRUE))
```

## Features

- **Estimation methods**: OLS, NLS, Poisson (MLE, GMM, EL), Negative
  Binomial (MLE, GMM, EL), iOLS
- **Bayesian workflows**: Poisson/NB models via `brms` and Stan, with
  posterior population-size summaries and chain diagnostics
- **Covariate-varying parameters**: $\alpha$ and $\beta$ via formula
  interface
- **Detection links**: `power`, `cloglog`, `logit`, and `probit`
- **Gamma offset**: estimated, fixed, or excluded
- **Constrained estimation**: for Poisson/NB fits with
  `constrained = TRUE`, $\alpha_i = \mathrm{logit}^{-1}(X_{\alpha,i} a)$
  and $\beta_i = \exp(X_{\beta,i} b)$, so fitted $\alpha_i \in (0,1)$
  and fitted $\beta_i > 0$
- **Robust inference**: HC0–HC5, cluster-robust via `sandwich`,
  fractional weighted bootstrap via `fwb`
- **Population size**: bias-corrected point estimates with delta-method
  or bootstrap CIs
- **Hypotheses for $\xi$**: threshold, change-over-time, grouped,
  bootstrap, and Bayesian posterior-proposition tests with explicit
  `xi[...]` syntax
- **Diagnostics**: `dfbeta()`, `dfpopsize()`, `loo()`, `rootogram()`,
  `profile_gamma()`, `dependence_bounds()`, `profile_dependence()`,
  `robustness_dependence()`, `exceedance_popsize()`, residual plots
- **Model comparison**: AIC/BIC, pseudo-$R^2$, likelihood ratio tests
- **Interactive app**: `run_app()` launches a Shiny dashboard

## Funding

This work is supported by the National Science Centre, OPUS 20 grant no.
2021/43/B/HS4/00469.

## References

- Zhang, L.-C. (2008). *Developing methods for determining the number of
  unauthorized foreigners in Norway* (Documents 2008/11). Statistics
  Norway.
  <https://www.ssb.no/a/english/publikasjoner/pdf/doc_200811_en/doc_200811_en.pdf>
- Beręsewicz, M., & Pawlukiewicz, K. (2020). Estimation of the number of
  irregular foreigners in Poland using non-linear count regression
  models. arXiv preprint arXiv:2008.09407.
- Santos Silva, J.M.C. and Tenreyro, S. (2006). The log of gravity. *The
  Review of Economics and Statistics*, 88(4), 641–658.
