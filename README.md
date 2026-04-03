
# The {uncounted} package <img src="man/figures/logo.png" align="right" width="150" />

<!-- badges: start -->

<!-- badges: end -->

R package for estimating the size of unauthorised (hidden) populations
using a power-law model that relates observed counts to reference
populations and auxiliary detection data.

## Model

For observation $i$, the expected observed count $m_i$ is modelled as:

$$E(m_i) = N_i^{\alpha_i} \cdot (\gamma + n_i / N_i)^{\beta_i}$$

where $N_i$ is the reference (total registered) population, $n_i$ is an
auxiliary count (e.g. new registrations), and $\gamma \geq 0$ is an
intercept-like offset. On the log scale the model is linear, enabling
estimation via OLS, NLS, Poisson PML, or Negative Binomial MLE.

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
## Fit Poisson PML with year x UKR interaction in alpha
fit <- estimate_hidden_pop(
  data = irregular_migration,
  observed = ~ m,
  auxiliary = ~ n,
  reference_pop = ~ N,
  method = "poisson",
  cov_alpha = ~ year + sex,
  cov_beta = ~ year,
  countries = ~ country_code
)

summary(fit)
#> Unauthorized population estimation
#> Method: POISSON | vcov: HC3 
#> N obs: 1382 
#> Gamma: 0.00725 (estimated) 
#> Log-likelihood: -8217.24 
#> AIC: 16462.49  BIC: 16535.73 
#> Deviance: 13831.69 
#> 
#> Coefficients:
#>                     Estimate Std. Error z value  Pr(>|z|)    
#> alpha:(Intercept)  0.7887405  0.0902787  8.7367 < 2.2e-16 ***
#> alpha:year2020     0.0014903  0.0958996  0.0155 0.9876014    
#> alpha:year2021    -0.0171693  0.0871702 -0.1970 0.8438568    
#> alpha:year2022    -0.0933541  0.1092056 -0.8548 0.3926354    
#> alpha:year2023    -0.1134401  0.1772828 -0.6399 0.5222492    
#> alpha:year2024    -0.1690320  0.1611266 -1.0491 0.2941488    
#> alpha:sexMale      0.0419908  0.0458827  0.9152 0.3600994    
#> beta:(Intercept)   0.6748933  0.1945975  3.4682 0.0005241 ***
#> beta:year2020      0.1895758  0.2273090  0.8340 0.4042807    
#> beta:year2021      0.2339935  0.2123637  1.1019 0.2705257    
#> beta:year2022      0.0896797  0.2485925  0.3607 0.7182865    
#> beta:year2023     -0.0516332  0.3783151 -0.1365 0.8914403    
#> beta:year2024     -0.3222270  0.3182389 -1.0125 0.3112839    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>                       Observed Estimate Estimate (BC) CI lower CI upper
#> year=2019, sex=Female    1,535   21,081        21,050    3,403  130,201
#> year=2019, sex=Male      5,069   59,923        59,840   10,229  350,054
#> year=2020, sex=Female      698   23,651        23,586    6,023   92,366
#> year=2020, sex=Male      2,700   66,767        66,576   22,922  193,364
#> year=2021, sex=Female      483   22,780        22,715    7,267   71,000
#> year=2021, sex=Male      2,622   65,115        64,921   32,333  130,353
#> year=2022, sex=Female      317   14,137        14,101    2,478   80,247
#> year=2022, sex=Male      2,632   32,956        32,878    7,367  146,727
#> year=2023, sex=Female      523   12,072        12,052      533  272,295
#> year=2023, sex=Male      3,839   28,853        28,809    1,219  680,842
#> year=2024, sex=Female      956    7,107         7,100      577   87,433
#> year=2024, sex=Male      5,731   17,461        17,447    1,177  258,640
```

### Population size estimates

``` r
## By year
popsize(fit, by = ~ year)
#>   group observed estimate estimate_bc     lower    upper share_pct
#> 1  2019     6604 81004.56    80890.35 14302.040 457504.6 21.781061
#> 2  2020     3398 90418.18    90161.89 30764.863 264235.4 24.312259
#> 3  2021     3105 87895.64    87636.36 42384.688 181200.6 23.633982
#> 4  2022     2949 47092.35    46978.74 10365.446 212919.1 12.662513
#> 5  2023     4362 40925.35    40860.84  1802.407 926321.7 11.004288
#> 6  2024     6687 24567.57    24547.14  1798.948 334952.4  6.605897
```

### Bootstrap confidence intervals

``` r
## Cluster FWB bootstrap (by country)
boot <- bootstrap_popsize(fit, R = 499,
                          cluster = ~ country_code,
                          by = ~ year)
boot
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

- **Estimation methods**: OLS, NLS, Poisson PML, Negative Binomial MLE
- **Covariate-varying parameters**: $\alpha$ and $\beta$ via formula
  interface
- **Gamma offset**: estimated, fixed, or excluded
- **Constrained estimation**: $\alpha \in (0,1)$ via logit, $\beta > 0$
  via log link
- **Robust inference**: HC0–HC5, cluster-robust via `sandwich`,
  fractional weighted bootstrap via `fwb`
- **Population size**: bias-corrected point estimates with delta-method
  or bootstrap CIs
- **Diagnostics**: `dfbeta()`, `dfpopsize()`, `loo()`, `rootogram()`,
  `profile_gamma()`, residual plots
- **Model comparison**: AIC/BIC, pseudo-$R^2$, likelihood ratio tests
- **Interactive app**: `run_app()` launches a Shiny dashboard

## Funding

This work is supported by the National Science Centre, OPUS 22 grant no.
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
