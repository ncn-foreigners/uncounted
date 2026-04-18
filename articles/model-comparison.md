# Model Comparison and Reporting

## Introduction

This vignette demonstrates how to compare estimation methods and present
results using the **uncounted** package’s plotting and reporting tools.

``` r
library(uncounted)
data(irregular_migration)
d <- irregular_migration
```

## Fitting multiple models

The package supports five estimation methods. We fit Poisson, NB, and
iOLS to the same data for comparison:

``` r
fit_po <- estimate_hidden_pop(d, ~ m, ~ n, ~ N,
  method = "poisson", cov_alpha = ~ year + sex,
  gamma = "estimate", countries = ~ country_code)

fit_nb <- estimate_hidden_pop(d, ~ m, ~ n, ~ N,
  method = "nb", cov_alpha = ~ year + sex,
  gamma = fit_po$gamma, countries = ~ country_code)

fit_io <- estimate_hidden_pop(d, ~ m, ~ n, ~ N,
  method = "iols", cov_alpha = ~ year + sex,
  gamma = fit_po$gamma, countries = ~ country_code)
```

## Coefficient tables with modelsummary

The package provides
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) methods
compatible with the **modelsummary** package:

``` r
library(modelsummary)
modelsummary(
  list(Poisson = fit_po, NB = fit_nb, iOLS = fit_io),
  stars = TRUE
)
```

|                                                              | Poisson      | NB           | iOLS         |
|--------------------------------------------------------------|--------------|--------------|--------------|
| alpha × (Intercept)                                          | 0.772\*\*\*  | 0.867\*\*\*  | 0.764\*\*\*  |
|                                                              | (0.050)      | (0.026)      | (0.065)      |
| alpha × year2020                                             | -0.062\*\*   | -0.112\*\*\* | -0.116\*\*\* |
|                                                              | (0.022)      | (0.016)      | (0.024)      |
| alpha × year2021                                             | -0.093\*\*\* | -0.148\*\*\* | -0.169\*\*\* |
|                                                              | (0.023)      | (0.019)      | (0.026)      |
| alpha × year2022                                             | -0.122\*\*\* | -0.127\*\*\* | -0.141\*\*\* |
|                                                              | (0.025)      | (0.018)      | (0.028)      |
| alpha × year2023                                             | -0.097\*     | -0.077\*\*\* | -0.092\*\*\* |
|                                                              | (0.046)      | (0.016)      | (0.021)      |
| alpha × year2024                                             | -0.060       | 0.008        | -0.008       |
|                                                              | (0.058)      | (0.016)      | (0.020)      |
| alpha × sexMale                                              | 0.042        | 0.017        | 0.049\*      |
|                                                              | (0.038)      | (0.012)      | (0.023)      |
| beta                                                         | 0.613\*\*\*  | 0.783\*\*\*  | 0.638\*\*\*  |
|                                                              | (0.122)      | (0.032)      | (0.071)      |
| gamma                                                        | 0.006        |              |              |
|                                                              | (0.005)      |              |              |
| theta                                                        |              | 1.276        |              |
|                                                              |              | (0.084)      |              |
| Num.Obs.                                                     | 1382         | 1382         | 1382         |
| AIC                                                          | 17569.0      | 5303.9       | 925616.0     |
| BIC                                                          | 17616.1      | 5351.0       | 925657.9     |
| Log.Lik.                                                     | -8775.519    | -2642.969    | -462800.005  |
| method                                                       | POISSON      | NB           | IOLS         |
| link_rho                                                     | power        | power        | power        |
| \+ p \< 0.1, \* p \< 0.05, \*\* p \< 0.01, \*\*\* p \< 0.001 |              |              |              |

Individual coefficient tables:

``` r
tidy(fit_po, conf.int = TRUE)
#>                term     estimate   std.error statistic      p.value
#> 1 alpha:(Intercept)  0.771782271 0.050387228 15.317022 5.885290e-53
#> 2    alpha:year2020 -0.061924082 0.021868386 -2.831671 4.630541e-03
#> 3    alpha:year2021 -0.093015705 0.022994522 -4.045125 5.229532e-05
#> 4    alpha:year2022 -0.122303174 0.024719166 -4.947706 7.509313e-07
#> 5    alpha:year2023 -0.096578133 0.045874341 -2.105276 3.526732e-02
#> 6    alpha:year2024 -0.059814566 0.057933918 -1.032462 3.018558e-01
#> 7     alpha:sexMale  0.041791940 0.037663368  1.109618 2.671638e-01
#> 8              beta  0.613267108 0.122280964  5.015230 5.297018e-07
#> 9             gamma  0.006084249 0.005101257        NA           NA
#>       conf.low    conf.high
#> 1  0.673025119  0.870539423
#> 2 -0.104785331 -0.019062833
#> 3 -0.138084140 -0.047947270
#> 4 -0.170751850 -0.073854498
#> 5 -0.186490190 -0.006666075
#> 6 -0.173362959  0.053733826
#> 7 -0.032026904  0.115610784
#> 8  0.373600822  0.852933394
#> 9  0.001176328  0.031469188
glance(fit_po)
#>    method estimator link_rho nobs    logLik      AIC      BIC deviance
#> 1 POISSON       MLE    power 1382 -8775.519 17569.04 17616.12 14948.24
#>   df.residual
#> 1        1373
```

## Plotting population size estimates

The
[`popsize()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
function returns an object that can be plotted directly:

``` r
# By year
ps <- popsize(fit_po, by = ~ year)
plot(ps)
```

![](model-comparison_files/figure-html/plot-popsize-1.png)

``` r

# Compare plug-in vs bias-corrected
plot(ps, type = "compare")
```

![](model-comparison_files/figure-html/plot-popsize-2.png)

## Comparing models visually

[`compare_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/compare_popsize.md)
produces side-by-side population size estimates:

``` r
comp <- compare_popsize(fit_po, fit_nb, fit_io,
  labels = c("Poisson", "NB", "iOLS"),
  by = ~ year)

print(comp)
#> Population size comparison: Poisson vs NB vs iOLS 
#> 
#>      model group  estimate estimate_bc      lower     upper
#> 1  Poisson  2019  67536.84    67504.29  29387.549 155059.83
#> 2  Poisson  2020  38787.40    38766.50  18790.590  79978.43
#> 3  Poisson  2021  33243.00    33224.95  16150.800  68349.37
#> 4  Poisson  2022  29249.61    29233.11  13739.451  62198.63
#> 5  Poisson  2023  40825.92    40804.83  11715.867 142117.88
#> 6  Poisson  2024  62218.64    62188.79  11951.535 323594.03
#> 7       NB  2019 156232.46   152016.63  94807.807 243746.34
#> 8       NB  2020  51814.02    50505.31  33151.834  76942.55
#> 9       NB  2021  42132.35    41172.73  25999.216  65201.72
#> 10      NB  2022  63811.39    62276.28  39962.923  97048.33
#> 11      NB  2023 116535.56   113818.31  71444.473 181324.14
#> 12      NB  2024 305601.61   298353.14 182612.609 487450.45
#> 13    iOLS  2019  65745.64    64612.50  22417.086 186231.84
#> 14    iOLS  2020  21998.78    21671.81   8247.433  56947.11
#> 15    iOLS  2021  15444.44    15246.53   5642.191  41199.73
#> 16    iOLS  2022  23655.48    23334.48   8887.972  61262.33
#> 17    iOLS  2023  41622.40    41061.39  15118.867 111518.78
#> 18    iOLS  2024 104625.34   103172.17  37066.214 287175.20
plot(comp)
```

![](model-comparison_files/figure-html/compare-popsize-1.png)

## Predictions

The [`predict()`](https://rdrr.io/r/stats/predict.html) method supports
new data:

``` r
# Fitted values
head(predict(fit_po))
#> [1]  3.2837247  0.3538453  2.3653567  2.2142433 21.0701268  0.8119226

# Log-scale (linear predictor)
head(predict(fit_po, type = "link"))
#> [1]  1.1889783 -1.0388956  0.8609288  0.7949107  3.0478562 -0.2083502

# Predictions for new data
new_d <- d[1:10, ]
predict(fit_po, newdata = new_d)
#>  [1]  3.2837247  0.3538453  2.3653567  2.2142433 21.0701268  0.8119226
#>  [7]  3.9496912  2.1480204  0.4921352 84.9170456
```

## Model selection guidance

| Situation                      | Recommended method             |
|--------------------------------|--------------------------------|
| Default analysis               | Poisson (robust, well-studied) |
| Overdispersion suspected       | NB, then LR test vs Poisson    |
| Sensitivity to large countries | iOLS alongside Poisson         |
| Quick exploration              | OLS with fixed gamma           |

When Poisson and iOLS agree, the result is robust to the choice of
weighting. When they diverge, report both and discuss which observations
drive the difference.

## References

- Santos Silva, J. M. C. & Tenreyro, S. (2006). The Log of Gravity.
  *Review of Economics and Statistics*, 88(4), 641–658.
- Benatia, D., Bellego, C. & Pape, L.-D. (2024). Dealing with Logs and
  Zeros in Regression Models. arXiv:2203.11820v3.
- Zhang, L.-C. (2008). *Developing methods for determining the number of
  unauthorized foreigners in Norway* (Documents 2008/11). Statistics
  Norway.
- Beresewicz, M. & Pawlukiewicz, K. (2020). Estimation of the number of
  irregular foreigners in Poland using non-linear count regression
  models. arXiv:2008.09407.
