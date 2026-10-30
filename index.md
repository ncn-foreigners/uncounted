# The {uncounted} package

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

``` math
\mu_i = E(m_i \mid N_i, n_i) = \xi_i \rho_i,
```

where $`\xi_i = E(M_i \mid N_i)`$ is the theoretical unauthorized
population size and $`\rho_i = E(p_i \mid N_i, n_i)`$ is the theoretical
detection rate.

In the baseline empirical specification implemented by default,

``` math
\xi_i = N_i^{\alpha_i}, \qquad
\rho_i = \left(\gamma_i + \frac{n_i}{N_i}\right)^{\beta_i},
```

where $`N_i`$ is the reference (total registered) population, $`n_i`$ is
an auxiliary count (e.g. police records), and $`\gamma_i \geq 0`$ is a
baseline detection offset. This gives

``` math
\mu_i = N_i^{\alpha_i}\left(\gamma_i + \frac{n_i}{N_i}\right)^{\beta_i}.
```

The package also supports bounded alternatives for the detection
component by writing

``` math
\eta_i = \beta_i \log\left(\gamma_i + \frac{n_i}{N_i}\right), \qquad
\rho_i = h(\eta_i),
```

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
#> Gamma: 0.00697 (estimated) 
#> Log-likelihood: -5738.09 
#> AIC: 11516.18  BIC: 11620.81 
#> Deviance: 8873.38 
#> 
#> Coefficients:
#>                       Estimate  Std. Error z value  Pr(>|z|)    
#> alpha:(Intercept)   0.82247146  0.14943954  5.5037 3.719e-08 ***
#> alpha:year2020     -0.01497293  0.17900376 -0.0836   0.93334    
#> alpha:year2021     -0.01825640  0.15231671 -0.1199   0.90460    
#> alpha:year2022     -0.01519229  0.14896112 -0.1020   0.91877    
#> alpha:year2023      0.03299293  0.15130742  0.2181   0.82739    
#> alpha:year2024      0.01078376  0.16501697  0.0653   0.94790    
#> alpha:ukr           0.00738549  0.07588025  0.0973   0.92246    
#> alpha:sexMale       0.01098736  0.03007214  0.3654   0.71484    
#> alpha:year2020:ukr  0.00913864  0.08771871  0.1042   0.91703    
#> alpha:year2021:ukr  0.00018515  0.07788140  0.0024   0.99810    
#> alpha:year2022:ukr -0.07628503  0.08542764 -0.8930   0.37187    
#> alpha:year2023:ukr -0.15366171  0.08741086 -1.7579   0.07876 .  
#> alpha:year2024:ukr -0.17562869  0.10224462 -1.7177   0.08585 .  
#> beta:(Intercept)    0.71379989  0.30000749  2.3793   0.01735 *  
#> beta:year2020       0.15778586  0.35639170  0.4427   0.65796    
#> beta:year2021       0.22676150  0.30420688  0.7454   0.45602    
#> beta:year2022       0.20168159  0.29788829  0.6770   0.49838    
#> beta:year2023       0.18681604  0.29982712  0.6231   0.53323    
#> beta:year2024      -0.01194996  0.32393841 -0.0369   0.97057    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>                              Observed Estimate Estimate (BC) CI lower CI upper
#> year=2019, ukr=0, sex=Female      426    9,786         9,767    1,142   83,563
#> year=2019, ukr=1, sex=Female    1,109   22,036        21,975    1,726  279,769
#> year=2019, ukr=0, sex=Male      2,177   24,520        24,470    2,305  259,806
#> year=2019, ukr=1, sex=Male      2,892   40,835        40,736    2,548  651,375
#> year=2020, ukr=0, sex=Female      193    9,700         9,668    1,818   51,404
#> year=2020, ukr=1, sex=Female      505   22,801        22,695    2,857  180,277
#> year=2020, ukr=0, sex=Male      1,158   24,237        24,144    4,238  137,550
#> year=2020, ukr=1, sex=Male      1,542   41,058        40,865    5,506  303,314
#> year=2021, ukr=0, sex=Female      144   11,647        11,613    4,347   31,028
#> year=2021, ukr=1, sex=Female      339   22,483        22,383    5,935   84,415
#> year=2021, ukr=0, sex=Male      1,337   29,891        29,789   12,289   72,208
#> year=2021, ukr=1, sex=Male      1,285   40,130        39,951   13,992  114,074
#> year=2022, ukr=0, sex=Female      211   14,362        14,335    7,102   28,937
#> year=2022, ukr=1, sex=Female      106   14,055        14,004    2,977   65,874
#> year=2022, ukr=0, sex=Male      1,903   38,815        38,726   21,687   69,152
#> year=2022, ukr=1, sex=Male        729   16,372        16,317    3,967   67,112
#> year=2023, ukr=0, sex=Female      444   24,167        24,134    9,994   58,280
#> year=2023, ukr=1, sex=Female       79    9,636         9,608    2,050   45,026
#> year=2023, ukr=0, sex=Male      3,283   70,211        70,095   30,282  162,249
#> year=2023, ukr=1, sex=Male        556   11,054        11,025    2,642   45,997
#> year=2024, ukr=0, sex=Female      887   22,457        22,431    5,746   87,560
#> year=2024, ukr=1, sex=Female       69    5,089         5,078      641   40,215
#> year=2024, ukr=0, sex=Male      5,041   62,267        62,188   14,162  273,066
#> year=2024, ukr=1, sex=Male        690    6,197         6,184      759   50,379
```

### Population size estimates

``` r

## By year
popsize(fit, by = ~ year)
#>   group observed  estimate estimate_bc     lower     upper share_pct
#> 1  2019     6604  97176.45    96948.51  8246.165 1139804.1  16.36501
#> 2  2020     3398  97796.22    97372.18 15272.198  620823.6  16.46939
#> 3  2021     3105 104151.16   103736.69 38793.630  277398.7  17.53959
#> 4  2022     2949  83604.67    83382.37 37643.431  184696.8  14.07946
#> 5  2023     4362 115068.21   114861.52 49218.327  268054.0  19.37808
#> 6  2024     6687  96009.35    95880.82 23187.477  396469.7  16.16847
```

### Hypothesis tests for population size

[`hypotheses_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/hypotheses_popsize.md)
tests threshold and contrast statements about grouped population-size
estimates. The expression uses explicit `xi[...]` filters, so the target
is unambiguous.

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

For frequentist
[`popsize()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
and
[`bootstrap_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/bootstrap_popsize.md)
objects, the helper reports Wald-style statistics and p values. Use
[`exceedance_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/exceedance_popsize.md)
when you want an empirical bootstrap tail area instead.

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

For Bayesian
[`popsize()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
objects, directional hypotheses are posterior propositions. The output
reports quantities such as `Pr(H1 | data)`, `Pr(H0 | data)`, posterior
odds, and credible intervals, not frequentist p values.

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
- **Covariate-varying parameters**: $`\alpha`$ and $`\beta`$ via formula
  interface
- **Detection links**: `power`, `cloglog`, `logit`, and `probit`
- **Gamma offset**: estimated, fixed, or excluded
- **Constrained estimation**: for Poisson/NB fits with
  `constrained = TRUE`,
  $`\alpha_i = \mathrm{logit}^{-1}(X_{\alpha,i} a)`$ and
  $`\beta_i = \exp(X_{\beta,i} b)`$, so fitted $`\alpha_i \in (0,1)`$
  and fitted $`\beta_i > 0`$
- **Robust inference**: HC0–HC5, cluster-robust via `sandwich`,
  fractional weighted bootstrap via `fwb`
- **Population size**: bias-corrected point estimates with delta-method
  or bootstrap CIs
- **Hypotheses for $`\xi`$**: threshold, change-over-time, grouped,
  bootstrap, and Bayesian posterior-proposition tests with explicit
  `xi[...]` syntax
- **Diagnostics**:
  [`dfbeta()`](https://rdrr.io/r/stats/influence.measures.html),
  [`dfpopsize()`](https://ncn-foreigners.github.io/uncounted/reference/dfpopsize.md),
  [`loo()`](https://ncn-foreigners.github.io/uncounted/reference/loo.md),
  [`rootogram()`](https://ncn-foreigners.github.io/uncounted/reference/rootogram.md),
  [`profile_gamma()`](https://ncn-foreigners.github.io/uncounted/reference/profile_gamma.md),
  [`dependence_bounds()`](https://ncn-foreigners.github.io/uncounted/reference/dependence_bounds.md),
  [`profile_dependence()`](https://ncn-foreigners.github.io/uncounted/reference/profile_dependence.md),
  [`robustness_dependence()`](https://ncn-foreigners.github.io/uncounted/reference/robustness_dependence.md),
  [`exceedance_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/exceedance_popsize.md),
  residual plots
- **Model comparison**: AIC/BIC, pseudo-$`R^2`$, likelihood ratio tests
- **Interactive app**:
  [`run_app()`](https://ncn-foreigners.github.io/uncounted/reference/run_app.md)
  launches a Shiny dashboard

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
