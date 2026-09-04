
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
#> Method: POISSON | estimator: MLE | link_rho: power | vcov: HC3 
#> N obs: 1382 
#> Gamma: 0.006979 (estimated) 
#> Log-likelihood: -5738.09 
#> AIC: 11516.18  BIC: 11620.81 
#> Deviance: 8873.39 
#> 
#> Coefficients:
#>                       Estimate  Std. Error z value  Pr(>|z|)    
#> alpha:(Intercept)   0.82280902  0.14938744  5.5079 3.632e-08 ***
#> alpha:year2020     -0.01556569  0.17904385 -0.0869   0.93072    
#> alpha:year2021     -0.01852846  0.15233174 -0.1216   0.90319    
#> alpha:year2022     -0.01554454  0.14896985 -0.1043   0.91689    
#> alpha:year2023      0.03271782  0.15129788  0.2162   0.82879    
#> alpha:year2024      0.01050503  0.16501614  0.0637   0.94924    
#> alpha:ukr           0.00728555  0.07588936  0.0960   0.92352    
#> alpha:sexMale       0.01099129  0.03008355  0.3654   0.71484    
#> alpha:year2020:ukr  0.00930099  0.08772054  0.1060   0.91556    
#> alpha:year2021:ukr  0.00026938  0.07788864  0.0035   0.99724    
#> alpha:year2022:ukr -0.07619193  0.08544340 -0.8917   0.37254    
#> alpha:year2023:ukr -0.15357545  0.08741910 -1.7568   0.07896 .  
#> alpha:year2024:ukr -0.17554573  0.10225718 -1.7167   0.08603 .  
#> beta:(Intercept)    0.71468425  0.30011634  2.3814   0.01725 *  
#> beta:year2020       0.15638361  0.35655031  0.4386   0.66095    
#> beta:year2021       0.22615449  0.30443173  0.7429   0.45756    
#> beta:year2022       0.20083922  0.29809353  0.6737   0.50047    
#> beta:year2023       0.18621054  0.29998544  0.6207   0.53478    
#> beta:year2024      -0.01259269  0.32410443 -0.0389   0.96901    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>                              Observed Estimate Estimate (BC) CI lower CI upper
#> year=2019, ukr=0, sex=Female      426    9,810         9,791    1,145   83,732
#> year=2019, ukr=1, sex=Female    1,109   22,099        22,038    1,729  280,903
#> year=2019, ukr=0, sex=Male      2,177   24,588        24,538    2,311  260,523
#> year=2019, ukr=1, sex=Male      2,892   40,960        40,861    2,550  654,746
#> year=2020, ukr=0, sex=Female      193    9,682         9,650    1,813   51,366
#> year=2020, ukr=1, sex=Female      505   22,748        22,642    2,854  179,606
#> year=2020, ukr=0, sex=Male      1,158   24,187        24,094    4,226  137,383
#> year=2020, ukr=1, sex=Male      1,542   40,960        40,767    5,507  301,815
#> year=2021, ukr=0, sex=Female      144   11,653        11,619    4,348   31,051
#> year=2021, ukr=1, sex=Female      339   22,496        22,397    5,937   84,492
#> year=2021, ukr=0, sex=Male      1,337   29,909        29,806   12,295   72,258
#> year=2021, ukr=1, sex=Male      1,285   40,158        39,979   14,000  114,160
#> year=2022, ukr=0, sex=Female      211   14,361        14,334    7,098   28,947
#> year=2022, ukr=1, sex=Female      106   14,051        14,000    2,975   65,887
#> year=2022, ukr=0, sex=Male      1,903   38,811        38,723   21,673   69,183
#> year=2022, ukr=1, sex=Male        729   16,369        16,313    3,964   67,128
#> year=2023, ukr=0, sex=Female      444   24,180        24,146   10,001   58,299
#> year=2023, ukr=1, sex=Female       79    9,642         9,615    2,052   45,054
#> year=2023, ukr=0, sex=Male      3,283   70,254        70,138   30,312  162,289
#> year=2023, ukr=1, sex=Male        556   11,062        11,032    2,644   46,024
#> year=2024, ukr=0, sex=Female      887   22,468        22,442    5,748   87,613
#> year=2024, ukr=1, sex=Female       69    5,092         5,081      641   40,247
#> year=2024, ukr=0, sex=Male      5,041   62,303        62,223   14,169  273,259
#> year=2024, ukr=1, sex=Male        690    6,200         6,188      759   50,422
```

### Population size estimates

``` r
## By year
popsize(fit, by = ~ year)
#>   group observed  estimate estimate_bc     lower     upper share_pct
#> 1  2019     6604  97456.70    97227.91  8261.575 1144245.1  16.40572
#> 2  2020     3398  97575.98    97153.06 15251.770  618860.5  16.42580
#> 3  2021     3105 104216.31   103801.46 38812.409  277610.8  17.54362
#> 4  2022     2949  83591.57    83369.29 37622.445  184741.8  14.07169
#> 5  2023     4362 115137.79   114930.93 49264.304  268127.6  19.38213
#> 6  2024     6687  96062.56    95933.92 23197.622  396735.4  16.17103
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

This work is supported by the National Science Centre, OPUS 27 grant no.
2024/53/B/HS4/01580.

## References

- Zhang, L.-C. (2008). *Developing methods for determining the number of
  unauthorized foreigners in Norway* (Documents 2008/11). Statistics
  Norway.
  <https://www.ssb.no/a/english/publikasjoner/pdf/doc_200811_en/doc_200811_en.pdf>
- Beręsewicz, M., & Pawlukiewicz, K. (2020). Estimation of the number of
  irregular foreigners in Poland using non-linear count regression
  models. arXiv preprint arXiv:2008.09407.
  <https://arxiv.org/abs/2008.09407>
- Santos Silva, J.M.C. and Tenreyro, S. (2006). The log of gravity. *The
  Review of Economics and Statistics*, 88(4), 641–658.
  <https://doi.org/10.1162/rest.88.4.641>
