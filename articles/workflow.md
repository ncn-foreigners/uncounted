# Complete Analysis Workflow

## Introduction

This vignette walks through a complete analysis estimating the
unauthorized foreign population in Poland using the
`irregular_migration` dataset shipped with the **uncounted** package.
The workflow covers every step from loading the data through fitting
models, obtaining population size estimates with confidence intervals,
running bootstrap inference, selecting among competing specifications,
and performing sensitivity analysis. By the end you will have a
reproducible pipeline that can be adapted to other datasets.

``` r

library(uncounted)
```

## Data Exploration

The `irregular_migration` dataset is a country-level panel (2019–2024)
with three count variables per observation:

- **m** – foreigners apprehended by the Border Guard for unauthorized
  stay (the quantity whose hidden population we want to estimate).
- **n** – foreigners identified by the Police (an auxiliary, partially
  overlapping administrative source).
- **N** – foreigners registered in the Social Insurance Institution
  (ZUS), used as a proxy for the total known foreign population from a
  given country.

``` r

data(irregular_migration)
str(irregular_migration)
#> 'data.frame':    1382 obs. of  8 variables:
#>  $ year        : Factor w/ 6 levels "2019","2020",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ sex         : Factor w/ 2 levels "Female","Male": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ country_code: chr  "AFG" "AGO" "ALB" "ARG" ...
#>  $ country     : chr  "Afghanistan" "Angola" "Albania" "Argentina" ...
#>  $ continent   : chr  "Asia" "Africa" "Europe" "Americas" ...
#>  $ m           : int  0 2 0 0 8 2 2 1 0 55 ...
#>  $ n           : int  0 0 1 1 17 0 1 1 0 59 ...
#>  $ N           : int  269 15 64 52 1089 44 219 47 23 11023 ...
```

A large share of observations have zero apprehensions, reflecting the
many small countries of origin with no recorded unauthorized stay.

``` r

cat("Share m == 0:", round(mean(irregular_migration$m == 0) * 100, 1), "%\n")
#> Share m == 0: 48.6 %
cat("Share n == 0:", round(mean(irregular_migration$n == 0) * 100, 1), "%\n")
#> Share n == 0: 48.3 %
```

A log–log scatter of the apprehension count against the reference
population reveals the power-law relationship that motivates the model.
Observations with `m > 0` and `n > 0` are shown.

``` r

pos <- irregular_migration$m > 0 & irregular_migration$n > 0
plot(
  log(irregular_migration$N[pos]),
  log(irregular_migration$m[pos]),
  xlab = "log(N)", ylab = "log(m)",
  main = "Apprehensions vs reference population (log-log)",
  pch = 16, cex = 0.6, col = adjustcolor("steelblue", 0.6)
)
```

![](workflow_files/figure-html/explore-plot-1.png)

## Fitting Models

We fit four specifications, all using the same formula interface. The
core model is

``` math
E(m_i) = N_i^{\alpha}\,(\gamma + n_i / N_i)^{\beta},
```

where $`\alpha`$ governs the elasticity with respect to the reference
population, $`\beta`$ captures the relationship with the auxiliary
detection rate, and $`\gamma`$ is a baseline offset ensuring the rate
term is positive even when $`n_i = 0`$.

### Model 1: Poisson (unconstrained, estimated gamma)

This is the recommended default. Poisson pseudo-maximum likelihood
(PPML) is consistent for the conditional mean even under
heteroscedasticity, and gamma is estimated jointly.

``` r

fit_pois <- estimate_hidden_pop(
  data = irregular_migration,
  observed = ~ m,
  auxiliary = ~ n,
  reference_pop = ~ N,
  method = "poisson",
  gamma = "estimate",
  countries = ~ country
)
summary(fit_pois)
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
```

Key output to look for:

- **alpha** – the elasticity with respect to `N`. Values below 1
  indicate that the unauthorized population grows less than
  proportionally with the reference population.
- **beta** – the elasticity with respect to the detection rate. Positive
  values mean higher auxiliary rates are associated with more
  apprehensions.
- **gamma** – the estimated baseline offset (typically small).

### Model 2: Poisson constrained (alpha in (0,1), beta \> 0)

Constraining alpha to the unit interval and beta to be positive enforces
theoretical expectations. With `constrained = TRUE`, the fitted
observation-level parameters are
`alpha_values = plogis(X_alpha %*% alpha_coefs)` and
`beta_values = exp(X_beta %*% beta_coefs)`, so fitted alpha stays in
`(0, 1)` and fitted beta stays positive. The reported coefficient
vectors remain on the transformed scales (logit for alpha, log for
beta).

``` r

fit_pois_c <- estimate_hidden_pop(
  data = irregular_migration,
  observed = ~ m,
  auxiliary = ~ n,
  reference_pop = ~ N,
  method = "poisson",
  gamma = "estimate",
  constrained = TRUE,
  countries = ~ country
)
summary(fit_pois_c)
#> Unauthorized population estimation
#> Method: POISSON | estimator: MLE | link_rho: power | vcov: HC3 
#> N obs: 1382 
#> Gamma: 0.001831 (estimated) 
#> Log-likelihood: -11380.65 
#> AIC: 22767.3  BIC: 22782.99 
#> Deviance: 20158.5 
#> 
#> Coefficients (link scale: logit for alpha, log for beta):
#>       Estimate Std. Error z value  Pr(>|z|)    
#> alpha  1.02861    0.20752  4.9568 7.167e-07 ***
#> beta  -0.55331    0.15974 -3.4637 0.0005328 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Response-scale fitted parameters (alpha_values in (0,1), beta_values > 0):
#>   Alpha (response scale):
#>        alpha SE(alpha)
#> (all) 0.7366    0.0403
#>   Beta (response scale):
#>    beta SE(beta)
#> 1 0.575   0.0919
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>       Observed Estimate Estimate (BC) CI lower CI upper
#> (all)   27,105  290,597       290,520  127,736  660,752
```

The [`summary()`](https://rdrr.io/r/base/summary.html) output reports
both link-scale coefficients (with standard errors and p-values) and
response-scale alpha/beta values with delta-method standard errors.

### Model 3: Negative Binomial (estimated gamma)

The NB model adds a dispersion parameter theta to accommodate
overdispersion beyond what Poisson allows.

``` r

fit_nb <- estimate_hidden_pop(
  data = irregular_migration,
  observed = ~ m,
  auxiliary = ~ n,
  reference_pop = ~ N,
  method = "nb",
  gamma = "estimate",
  countries = ~ country
)
summary(fit_nb)
#> Unauthorized population estimation
#> Method: NB | estimator: MLE | link_rho: power | vcov: HC1 
#> N obs: 1382 
#> Gamma: 0.01901 (estimated) 
#> Theta (NB dispersion): 1.0617 
#> Log-likelihood: -2693.48 
#> AIC: 5394.97  BIC: 5415.89 
#> Deviance: 1272.2 
#> 
#> Coefficients:
#>       Estimate Std. Error z value  Pr(>|z|)    
#> alpha 0.871633   0.023597  36.939 < 2.2e-16 ***
#> beta  1.034268   0.080108  12.911 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>       Observed  Estimate Estimate (BC) CI lower  CI upper
#> (all)   27,105 1,259,529     1,231,371  731,343 2,073,274
```

### Model 4: NB two-stage (gamma from Poisson, fixed in NB)

A pragmatic two-stage approach: first estimate gamma from the Poisson
model (which is better identified), then fix it in the NB model. This
avoids potential numerical instability when estimating gamma and theta
simultaneously.

``` r

gamma_from_pois <- fit_pois$gamma

fit_nb2 <- estimate_hidden_pop(
  data = irregular_migration,
  observed = ~ m,
  auxiliary = ~ n,
  reference_pop = ~ N,
  method = "nb",
  gamma = gamma_from_pois,
  countries = ~ country
)
summary(fit_nb2)
#> Unauthorized population estimation
#> Method: NB | estimator: MLE | link_rho: power | vcov: HC1 
#> N obs: 1382 
#> Gamma: 0.001831 (fixed) 
#> Theta (NB dispersion): 0.9491 
#> Log-likelihood: -2727.07 
#> AIC: 5460.14  BIC: 5475.84 
#> Deviance: 1267.27 
#> 
#> Coefficients:
#>       Estimate Std. Error z value  Pr(>|z|)    
#> alpha 0.745560   0.015958  46.721 < 2.2e-16 ***
#> beta  0.597130   0.023991  24.890 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>       Observed Estimate Estimate (BC) CI lower CI upper
#> (all)   27,105  318,959       315,416  227,270  437,751
```

## Population Size Estimates

The
[`popsize()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
function extracts the estimated total unauthorized population
$`\hat{\xi} = \sum_i N_i^{\hat{\alpha}}`$ along with bias-corrected
estimates and confidence intervals.

``` r

ps_pois <- popsize(fit_pois)
ps_pois
#>   group observed estimate estimate_bc    lower  upper share_pct
#> 1 (all)    27105 290597.4    290498.3 127726.4 660703       100
```

The columns are:

- **group** – label for each alpha-covariate group (here `"(all)"` since
  we have no covariates on alpha).
- **observed** – total observed count $`\sum m_i`$.
- **estimate** – plug-in estimate $`\hat{\xi}`$.
- **estimate_bc** – analytical bias-corrected estimate. For
  unconstrained fits this is the multiplicative lognormal correction;
  constrained fits use a second-order Taylor approximation.
- **lower**, **upper** – confidence interval bounds. The CI is
  constructed by a delta-method/log-normal approximation on the
  population-size scale. When bias correction is requested, the bounds
  are rescaled by the same `estimate_bc / estimate` ratio as the
  reported point estimate.
- **share_pct** – group share as a percentage of the total.

Comparing across all four models:

``` r

models <- list(
  Poisson       = fit_pois,
  Poisson_constr = fit_pois_c,
  NB            = fit_nb,
  NB_twostage   = fit_nb2
)

ps_table <- do.call(rbind, lapply(names(models), function(nm) {
  ps <- popsize(models[[nm]])
  data.frame(
    Model    = nm,
    Estimate = sum(ps$estimate),
    BC       = sum(ps$estimate_bc),
    Lower    = sum(ps$lower),
    Upper    = sum(ps$upper)
  )
}))
ps_table
#>            Model  Estimate        BC    Lower     Upper
#> 1        Poisson  290597.4  290498.3 127726.4  660703.0
#> 2 Poisson_constr  290597.4  290520.0 127736.0  660752.5
#> 3             NB 1259529.5 1231370.7 731342.8 2073273.8
#> 4    NB_twostage  318959.1  315416.4 227269.8  437750.8
```

### Stratified population size

The `by` parameter allows computing population size by any variable in
the data, not just the `cov_alpha` groups. This is useful when alpha
varies by sex but you want totals by year, continent, or country.

``` r

# Population size by year (using the same fitted alpha)
popsize(fit_pois, by = ~ factor(year))
#>   group observed estimate estimate_bc    lower     upper share_pct
#> 1  2019     6604 34740.14    34728.95 15668.65  76975.38  11.95473
#> 2  2020     3398 37922.60    37910.18 16980.30  84638.15  13.04987
#> 3  2021     3105 44585.24    44570.30 19761.07 100526.51  15.34261
#> 4  2022     2949 54737.13    54717.54 23564.68 127054.96  18.83607
#> 5  2023     4362 58339.61    58319.10 25297.62 134444.16  20.07575
#> 6  2024     6687 60272.70    60252.18 26470.28 137147.20  20.74096

# By continent
popsize(fit_pois, by = ~ continent)
#>               group observed    estimate estimate_bc     lower       upper
#> 1            Africa     1304  13815.9904  13814.5986  8830.439  21611.9656
#> 2          Americas      931  11053.0193  11051.7084  6830.839  17880.7098
#> 3              Asia    10968  69969.5494  69956.2380 37521.505 130428.5449
#> 4            Europe    13826 192230.5202 192147.9773 74709.948 494189.1420
#> 5           Oceania       14    441.4146    441.3908   317.687    613.2635
#> 6 unknown/stateless       62   3086.9284   3086.3373  1635.241   5825.1204
#>   share_pct
#> 1  4.754340
#> 2  3.803550
#> 3 24.077829
#> 4 66.150112
#> 5  0.151899
#> 6  1.062270

# By country (one estimate per country)
popsize(fit_pois, by = ~ country)
#>                           group observed     estimate  estimate_bc        lower
#> 1                   Afghanistan      105 9.119919e+02 9.118977e+02   571.839576
#> 2                       Albania       48 6.449180e+02 6.448572e+02   414.172947
#> 3                       Algeria      134 9.401769e+02 9.400627e+02   569.178290
#> 4                       Andorra        0 2.250028e+01 2.250016e+01    20.287782
#> 5                        Angola       16 2.742566e+02 2.742407e+02   193.851914
#> 6                     Argentina       55 6.756970e+02 6.756310e+02   430.732894
#> 7                       Armenia      178 2.646683e+03 2.646260e+03  1480.598061
#> 8                         Aruba        0 2.246309e+00 2.246301e+00     2.059782
#> 9                     Australia       11 3.152166e+02 3.151966e+02   219.262332
#> 10                   Azerbaijan      272 2.410165e+03 2.409763e+03  1335.419422
#> 11                      Bahamas        0 2.776547e+00 2.776531e+00     2.488844
#> 12                   Bangladesh      186 2.392319e+03 2.391875e+03  1285.424916
#> 13                     Barbados        0 1.578149e+01 1.578138e+01    13.927294
#> 14                      Belarus     1415 3.134118e+04 3.133029e+04 13318.917849
#> 15                       Belize        0 2.350442e+01 2.350364e+01    18.031632
#> 16                        Benin        1 1.802780e+01 1.802768e+01    15.985784
#> 17                       Bhutan        0 4.492617e+00 4.492601e+00     4.119563
#> 18                      Bolivia        2 1.402283e+02 1.402229e+02   106.007404
#> 19         Bosnia & Herzegovina        4 2.399502e+02 2.399357e+02   168.320120
#> 20                     Botswana        0 1.592818e+01 1.592804e+01    13.953433
#> 21                       Brazil       54 1.359581e+03 1.359411e+03   814.084064
#> 22       British Virgin Islands        0 1.141252e+01 1.141244e+01    10.146022
#> 23                 Burkina Faso        7 1.481437e+01 1.481429e+01    13.365767
#> 24                      Burundi        2 6.584412e+01 6.584216e+01    51.339923
#> 25                     Cambodia        0 4.799067e+01 4.799013e+01    41.224287
#> 26                     Cameroon       37 4.463643e+02 4.463305e+02   299.909417
#> 27                       Canada       15 4.742112e+02 4.741748e+02   317.488408
#> 28                   Cape Verde        0 6.442844e+00 6.442826e+00     6.059840
#> 29     Central African Republic        0 1.519824e+01 1.519809e+01    13.302839
#> 30                         Chad        1 1.171201e+01 1.171196e+01    10.672976
#> 31                        Chile        7 2.501084e+02 2.500937e+02   176.777898
#> 32                        China      181 3.314421e+03 3.313851e+03  1813.217200
#> 33                     Colombia      634 2.650544e+03 2.650018e+03  1397.416679
#> 34                      Comoros        0 4.193127e+00 4.193080e+00     3.596265
#> 35          Congo - Brazzaville       15 4.130729e+02 4.130429e+02   279.750456
#> 36             Congo - Kinshasa        3 1.000000e+00 1.000000e+00           NA
#> 37                   Costa Rica        0 2.303829e+02 2.303710e+02   166.184414
#> 38                Cote d'Ivoire        3 1.116227e+02 1.116193e+02    87.087292
#> 39                         Cuba       28 5.112563e+02 5.112138e+02   337.146718
#> 40                       Cyprus        0 1.234401e+02 1.234355e+02    93.908933
#> 41           Dominican Republic        2 1.041464e+02 1.041434e+02    81.777883
#> 42                      Ecuador        3 1.915230e+02 1.915144e+02   141.139914
#> 43                        Egypt      120 1.353337e+03 1.353144e+03   785.371553
#> 44                  El Salvador        1 7.093653e+01 7.093498e+01    57.666408
#> 45                      Eritrea        1 2.021741e+01 2.021717e+01    17.329041
#> 46                     Ethiopia       51 7.068322e+02 7.067630e+02   449.735482
#> 47  French Southern Territories        0 1.377164e+01 1.377145e+01    11.669258
#> 48                        Gabon        0 2.950714e+01 2.950690e+01    26.075181
#> 49                       Gambia       15 8.410462e+01 8.410175e+01    64.710091
#> 50                      Georgia     4429 1.091515e+04 1.091224e+04  5162.549301
#> 51                        Ghana       36 3.472638e+02 3.472386e+02   236.157686
#> 52                      Grenada        0 2.666297e+00 2.666295e+00     2.576702
#> 53                    Guatemala        7 2.217124e+02 2.216993e+02   156.563353
#> 54                       Guinea       20 1.107083e+02 1.107046e+02    85.388102
#> 55                       Guyana        0 2.234074e+01 2.234059e+01    19.910265
#> 56                        Haiti        0 3.452667e+01 3.452615e+01    29.105424
#> 57                     Honduras        1 5.787925e+01 5.787799e+01    46.905954
#> 58          Hong Kong SAR China        0 2.246309e+00 2.246301e+00     2.059782
#> 59                        India      826 7.781252e+03 7.779331e+03  3790.380816
#> 60                    Indonesia      246 2.951829e+03 2.951298e+03  1598.031613
#> 61                         Iran       60 8.589063e+02 8.588175e+02   539.064953
#> 62                         Iraq      178 5.354914e+02 5.354431e+02   347.751179
#> 63                      Ireland        0 6.531085e+02 6.530445e+02   416.092442
#> 64                       Israel       13 5.091177e+02 5.090746e+02   334.790059
#> 65                      Jamaica        2 5.590924e+01 5.590836e+01    46.863783
#> 66                        Japan        5 6.574955e+02 6.574374e+02   427.033338
#> 67                       Jordan       24 4.299589e+02 4.299238e+02   285.143465
#> 68                   Kazakhstan      155 2.462590e+03 2.462203e+03  1384.897096
#> 69                        Kenya       36 6.650419e+02 6.649792e+02   426.687269
#> 70                       Kuwait        0 1.730153e+01 1.730134e+01    15.084623
#> 71                   Kyrgyzstan      171 1.194060e+03 1.193904e+03   709.127630
#> 72                         Laos        0 3.712930e+01 3.712870e+01    31.066775
#> 73                      Lebanon       18 4.729088e+02 4.728692e+02   311.893804
#> 74                      Liberia        1 8.185209e+00 8.185174e+00     7.512948
#> 75                        Libya       16 1.750429e+02 1.750330e+02   124.405852
#> 76              Macao SAR China        0 5.849754e+00 5.849655e+00     4.841283
#> 77                   Madagascar        0 6.216191e+01 6.216057e+01    50.498101
#> 78                       Malawi        0 9.392338e+00 9.392194e+00     7.885198
#> 79                     Malaysia        7 2.107511e+02 2.107413e+02   154.283115
#> 80                         Mali        8 1.673425e+02 1.673341e+02   121.590816
#> 81                   Mauritania        0 2.870659e+01 2.870633e+01    25.161600
#> 82                    Mauritius       14 1.039031e+02 1.038997e+02    80.357419
#> 83                       Mexico       20 9.402925e+02 9.401917e+02   585.099080
#> 84                      Moldova     1386 7.641567e+03 7.639799e+03  3800.245707
#> 85                     Mongolia       60 6.424446e+02 6.423889e+02   418.748884
#> 86                   Montenegro        2 9.889813e+01 9.889422e+01    74.423728
#> 87                      Morocco      115 9.787625e+02 9.786489e+02   598.620194
#> 88                   Mozambique        3 8.918056e+01 8.917803e+01    70.179498
#> 89              Myanmar (Burma)        1 6.036308e+01 6.036160e+01    48.405371
#> 90                      Namibia        1 1.595491e+01 1.595481e+01    14.318543
#> 91                        Nepal      228 3.241028e+03 3.240418e+03  1731.461397
#> 92                  New Zealand        3 1.186191e+02 1.186153e+02    91.808732
#> 93                    Nicaragua        2 7.124506e+01 7.124318e+01    56.916576
#> 94                        Niger        1 2.221460e+02 2.221339e+02   159.128445
#> 95                      Nigeria      207 1.325439e+03 1.325264e+03   784.357907
#> 96              North Macedonia       14 3.627555e+02 3.627306e+02   248.644037
#> 97                         Oman        0 1.575925e+01 1.575888e+01    12.629731
#> 98                     Pakistan      154 1.221848e+03 1.221692e+03   729.220810
#> 99      Palestinian Territories       14 2.628414e+02 2.628231e+02   179.858169
#> 100                      Panama        0 5.674112e+01 5.673996e+01    46.589976
#> 101            Papua New Guinea        0 5.332595e+00 5.332590e+00     5.153404
#> 102                    Paraguay       11 7.899192e+01 7.898902e+01    60.535272
#> 103                        Peru       18 5.090879e+02 5.090458e+02   336.293340
#> 104                 Philippines      285 4.713271e+03 4.712307e+03  2447.312333
#> 105            Pitcairn Islands        0 2.246309e+00 2.246301e+00     2.059782
#> 106                       Qatar        0 2.262079e+01 2.262006e+01    17.433074
#> 107                      Russia     1020 7.076379e+03 7.074810e+03  3569.047896
#> 108                      Rwanda       61 7.785386e+02 7.784587e+02   490.004340
#> 109         Sao Tome & Principe        0 5.332595e+00 5.332590e+00     5.153404
#> 110                Saudi Arabia       23 9.920814e+01 9.920516e+01    77.655061
#> 111                     Senegal       11 1.479975e+02 1.479911e+02   110.096118
#> 112                      Serbia       17 6.233232e+02 6.232426e+02   370.844052
#> 113                Sierra Leone        2 4.006803e+01 4.006746e+01    33.849191
#> 114                   Singapore        0 6.603344e+01 6.603232e+01    54.753804
#> 115                     Somalia       11 6.710474e+01 6.710306e+01    53.556015
#> 116                South Africa       28 4.797127e+02 4.796766e+02   322.127289
#> 117                 South Korea       33 6.827031e+02 6.826338e+02   431.795630
#> 118                   Sri Lanka       57 5.857001e+02 5.856378e+02   368.177612
#> 119                   St. Lucia        0 2.776547e+00 2.776531e+00     2.488844
#> 120                       Sudan       16 1.179993e+02 1.179948e+02    89.355972
#> 121                    Suriname        0 3.838724e+01 3.838683e+01    33.128144
#> 122                       Syria      228 7.655524e+02 7.654715e+02   479.115676
#> 123                      Taiwan        6 2.308852e+02 2.308740e+02   167.776816
#> 124                  Tajikistan      393 8.820419e+02 8.819377e+02   538.453544
#> 125                    Tanzania       17 2.635500e+02 2.635347e+02   186.114602
#> 126                    Thailand       34 9.185073e+02 9.184109e+02   574.154201
#> 127                        Togo        1 5.328275e+01 5.328140e+01    42.386593
#> 128           Trinidad & Tobago        0 1.652514e+01 1.652499e+01    14.457749
#> 129                     Tunisia       85 1.002552e+03 1.002432e+03   608.845897
#> 130                      Turkey      394 4.366647e+03 4.365740e+03  2258.470030
#> 131                Turkmenistan      385 1.271760e+03 1.271571e+03   731.401951
#> 132                      Uganda       20 3.124546e+02 3.124346e+02   216.974489
#> 133                     Ukraine     9901 1.412361e+05 1.411684e+05 51604.797919
#> 134        United Arab Emirates        1 3.743018e+00 3.742983e+00     3.249516
#> 135              United Kingdom       19 2.261041e+03 2.260677e+03  1265.107261
#> 136               United States       38 1.589667e+03 1.589455e+03   935.240653
#> 137           unknown/stateless       62 3.086928e+03 3.086337e+03  1635.241327
#> 138                     Uruguay        1 1.334334e+02 1.334289e+02   102.547257
#> 139                  Uzbekistan      925 3.110260e+03 3.109649e+03  1642.317589
#> 140                Vatican City        0 1.503438e+01 1.503431e+01    13.608365
#> 141                   Venezuela       30 5.065006e+02 5.064594e+02   335.498987
#> 142                     Vietnam      677 5.685574e+03 5.684390e+03  2929.983659
#> 143                       Yemen       16 2.272124e+02 2.271979e+02   158.129227
#> 144                      Zambia        8 1.585480e+02 1.585414e+02   118.370698
#> 145                    Zimbabwe      180 1.546966e+03 1.546751e+03   901.664646
#>            upper    share_pct
#> 1   1.454180e+03 3.138335e-01
#> 2   1.004027e+03 2.219283e-01
#> 3   1.552621e+03 3.235324e-01
#> 4   2.495379e+01 7.742767e-03
#> 5   3.879661e+02 9.437683e-02
#> 6   1.059769e+03 2.325199e-01
#> 7   4.729637e+03 9.107730e-01
#> 8   2.449709e+00 7.729967e-04
#> 9   4.531052e+02 1.084719e-01
#> 10  4.348414e+03 8.293828e-01
#> 11  3.097472e+00 9.554616e-04
#> 12  4.450721e+03 8.232418e-01
#> 13  1.788228e+01 5.430707e-03
#> 14  7.369872e+04 1.078509e+01
#> 15  3.063621e+01 8.088310e-03
#> 16  2.033038e+01 6.203704e-03
#> 17  4.899418e+00 1.545993e-03
#> 18  1.854819e+02 4.825516e-02
#> 19  3.420217e+02 8.257136e-02
#> 20  1.818207e+01 5.481185e-03
#> 21  2.270032e+03 4.678570e-01
#> 22  1.283692e+01 3.927260e-03
#> 23  1.641979e+01 5.097902e-03
#> 24  8.444090e+01 2.265819e-02
#> 25  5.586640e+01 1.651449e-02
#> 26  6.642370e+02 1.536023e-01
#> 27  7.081887e+02 1.631849e-01
#> 28  6.850017e+00 2.217103e-03
#> 29  1.736336e+01 5.229997e-03
#> 30  1.285208e+01 4.030321e-03
#> 31  3.538162e+02 8.606698e-02
#> 32  6.056421e+03 1.140554e+00
#> 33  5.025412e+03 9.121018e-01
#> 34  4.888938e+00 1.442933e-03
#> 35  6.098450e+02 1.421461e-01
#> 36            NA 3.441187e-04
#> 37  3.193489e+02 7.927905e-02
#> 38  1.430619e+02 3.841145e-02
#> 39  7.751507e+02 1.759329e-01
#> 40  1.622458e+02 4.247805e-02
#> 41  1.326258e+02 3.583873e-02
#> 42  2.598682e+02 6.590664e-02
#> 43  2.331379e+03 4.657085e-01
#> 44  8.725655e+01 2.441058e-02
#> 45  2.358665e+01 6.957188e-03
#> 46  1.110684e+03 2.432342e-01
#> 47  1.625236e+01 4.739077e-03
#> 48  3.339027e+01 1.015396e-02
#> 49  1.093045e+02 2.894197e-02
#> 50  2.306554e+04 3.756108e+00
#> 51  5.105684e+02 1.195000e-01
#> 52  2.759003e+00 9.175227e-04
#> 53  3.139341e+02 7.629538e-02
#> 54  1.435270e+02 3.809678e-02
#> 55  2.506758e+01 7.687868e-03
#> 56  4.095645e+01 1.188127e-02
#> 57  7.141655e+01 1.991733e-02
#> 58  2.449709e+00 7.729967e-04
#> 59  1.596620e+04 2.677674e+00
#> 60  5.450554e+03 1.015780e+00
#> 61  1.368235e+03 2.955657e-01
#> 62  8.244382e+02 1.842726e-01
#> 63  1.024934e+03 2.247468e-01
#> 64  7.740881e+02 1.751969e-01
#> 65  6.669851e+01 1.923941e-02
#> 66  1.012155e+03 2.262565e-01
#> 67  6.482156e+02 1.479569e-01
#> 68  4.377542e+03 8.474234e-01
#> 69  1.036350e+03 2.288533e-01
#> 70  1.984382e+01 5.953781e-03
#> 71  2.010084e+03 4.108985e-01
#> 72  4.437347e+01 1.277689e-02
#> 73  7.169277e+02 1.627367e-01
#> 74  8.917547e+00 2.816683e-03
#> 75  2.462630e+02 6.023553e-02
#> 76  7.068058e+00 2.013010e-03
#> 77  7.651646e+01 2.139107e-02
#> 78  1.118720e+01 3.232079e-03
#> 79  2.878599e+02 7.252340e-02
#> 80  2.302864e+02 5.758570e-02
#> 81  3.275043e+01 9.878473e-03
#> 82  1.343392e+02 3.575501e-02
#> 83  1.510788e+03 3.235722e-01
#> 84  1.535862e+04 2.629606e+00
#> 85  9.854676e+02 2.210772e-01
#> 86  1.314106e+02 3.403269e-02
#> 87  1.599936e+03 3.368104e-01
#> 88  1.133197e+02 3.068870e-02
#> 89  7.527105e+01 2.077206e-02
#> 90  1.777807e+01 5.490382e-03
#> 91  6.064419e+03 1.115298e+00
#> 92  1.532488e+02 4.081904e-02
#> 93  8.917597e+01 2.451675e-02
#> 94  3.100857e+02 7.644457e-02
#> 95  2.239187e+03 4.561082e-01
#> 96  5.291642e+02 1.248309e-01
#> 97  1.966332e+01 5.423052e-03
#> 98  2.046747e+03 4.204606e-01
#> 99  3.840580e+02 9.044863e-02
#> 100 6.910120e+01 1.952568e-02
#> 101 5.518007e+00 1.835045e-03
#> 102 1.030683e+02 2.718260e-02
#> 103 7.705406e+02 1.751867e-01
#> 104 9.073561e+03 1.621925e+00
#> 105 2.449709e+00 7.729967e-04
#> 106 2.935037e+01 7.784236e-03
#> 107 1.402417e+04 2.435114e+00
#> 108 1.236719e+03 2.679097e-01
#> 109 5.518007e+00 1.835045e-03
#> 110 1.267356e+02 3.413937e-02
#> 111 1.989294e+02 5.092869e-02
#> 112 1.047425e+03 2.144972e-01
#> 113 4.742806e+01 1.378816e-02
#> 114 7.963404e+01 2.272334e-02
#> 115 8.407685e+01 2.309199e-02
#> 116 7.142817e+02 1.650781e-01
#> 117 1.079188e+03 2.349309e-01
#> 118 9.315385e+02 2.015504e-01
#> 119 3.097472e+00 9.554616e-04
#> 120 1.558124e+02 4.060576e-02
#> 121 4.448028e+01 1.320977e-02
#> 122 1.222975e+03 2.634409e-01
#> 123 3.177007e+02 7.945191e-02
#> 124 1.444533e+03 3.035271e-01
#> 125 3.731599e+02 9.069249e-02
#> 126 1.469080e+03 3.160755e-01
#> 127 6.697654e+01 1.833559e-02
#> 128 1.888781e+01 5.686608e-03
#> 129 1.650449e+03 3.449967e-01
#> 130 8.439204e+03 1.502645e+00
#> 131 2.210674e+03 4.376364e-01
#> 132 4.498932e+02 1.075215e-01
#> 133 3.861756e+05 4.860198e+01
#> 134 4.311386e+00 1.288042e-03
#> 135 4.039705e+03 7.780665e-01
#> 136 2.701300e+03 5.470342e-01
#> 137 5.825120e+03 1.062270e+00
#> 138 1.736104e+02 4.591691e-02
#> 139 5.887970e+03 1.070298e+00
#> 140 1.660968e+01 5.173612e-03
#> 141 7.645361e+02 1.742963e-01
#> 142 1.102815e+04 1.956512e+00
#> 143 3.264349e+02 7.818803e-02
#> 144 2.123445e+02 5.455932e-02
#> 145 2.653356e+03 5.323400e-01
```

The total across all `by`-groups equals the total from the default
grouping.

## Bootstrap Inference

The
[`bootstrap_popsize()`](https://ncn-foreigners.github.io/uncounted/reference/bootstrap_popsize.md)
function uses the fractional weighted bootstrap (FWB) of Xu et
al. (2020), via the `fwb` package (Greifer, 2026). Cluster bootstrap by
country is recommended to account for within-country correlation across
years. When refits are slow, `ci_type = "cheap"` (Lam, 2022) gives valid
intervals from very few replications.

``` r

set.seed(2025)
boot_pois <- bootstrap_popsize(
  fit_pois,
  R = 49,
  cluster = ~ country_code,
  ci_type = "perc",
  seed = 2025
)
boot_pois
#> Bootstrap population size estimation
#> R = 49 | CI type: perc | Point estimate: median | Converged: 49 / 49 
#> Cluster bootstrap
#> 95% CI
#> 
#>   Point estimate: bootstrap median (recommended) | CI: bootstrap percentile
#>        Plugin Plugin (BC) Boot median Boot mean CI lower CI upper
#> (all) 290,597     290,498     289,317   344,128  183,808  728,846
# #> Bootstrap population size estimation
# #> R = 49 | CI type: perc | Point estimate: median | Converged: 998 / 999
# #> Cluster bootstrap
# #> 95% CI
```

Access the detailed results:

``` r

# Main table (selected point estimate + CI)
boot_pois$popsize
#>   group estimate    lower  upper
#> 1 (all) 289317.3 183807.5 728846

# Full table with all point estimate types
boot_pois$popsize_full
#>   group   plugin plugin_bc boot_median boot_mean    lower  upper
#> 1 (all) 290597.4  290498.3    289317.3    344128 183807.5 728846

# Bootstrap distribution summary
summary(boot_pois)
#> Bootstrap population size estimation
#> R = 49 | CI type: perc | Point estimate: median | Converged: 49 / 49 
#> Cluster bootstrap
#> 95% CI
#> 
#>   Point estimate: bootstrap median (recommended) | CI: bootstrap percentile
#>        Plugin Plugin (BC) Boot median Boot mean CI lower CI upper
#> (all) 290,597     290,498     289,317   344,128  183,808  728,846
#> 
#> Bootstrap distribution summary:
#>         mean     sd q025.2.5% q50.50% q975.97.5%
#> (all) 344128 150262    183808  289317     728846
#> 
#> Parameter bootstrap summary:
#>           mean       sd q025.2.5%  q50.50% q975.97.5%
#> alpha 0.744933 0.036480  0.691930 0.736222   0.822479
#> beta  0.602218 0.098744  0.462785 0.588591   0.842640
#> gamma 0.002686 0.002476  0.000570 0.001978   0.008180
```

### Comparing bootstrap and analytical CIs

``` r

ps_analytical <- popsize(fit_pois)
data.frame(
  Type       = c("Analytical", "Bootstrap"),
  Estimate   = c(sum(ps_analytical$estimate_bc),
                 boot_pois$popsize$estimate),
  Lower      = c(sum(ps_analytical$lower),
                 boot_pois$popsize$lower),
  Upper      = c(sum(ps_analytical$upper),
                 boot_pois$popsize$upper)
)
#>         Type Estimate    Lower  Upper
#> 1 Analytical 290498.3 127726.4 660703
#> 2  Bootstrap 289317.3 183807.5 728846
```

The bootstrap median is typically smaller than the plug-in estimate.
This is expected: $`\xi(\alpha) = \sum N_i^\alpha`$ is convex in
$`\alpha`$ for $`N_i > 1`$, so by Jensen’s inequality
$`E[\hat{\xi}] \geq \xi`$. The bootstrap mean inherits this upward bias,
while the bootstrap median is more robust to the right skew of the
bootstrap distribution and is therefore the recommended point estimate.

## Model Selection

### Information criteria and fit statistics

``` r

comp <- compare_models(
  Poisson        = fit_pois,
  Poisson_constr = fit_pois_c,
  NB             = fit_nb,
  NB_twostage    = fit_nb2,
  sort_by = "AIC"
)
comp
#> Model comparison
#> ------------------------------------------------------------ 
#>           Model  Method Estimator  Link Constrained n_par    logLik      AIC
#>              NB      NB       MLE power       FALSE     4  -2693.48  5394.97
#>     NB_twostage      NB       MLE power       FALSE     3  -2727.07  5460.14
#>         Poisson POISSON       MLE power       FALSE     3 -11380.65 22767.30
#>  Poisson_constr POISSON       MLE power        TRUE     3 -11380.65 22767.30
#>       BIC Deviance Pearson_X2   RMSE R2_cor R2_D  R2_CW
#>   5415.89  1272.20    2943.87 139.85 0.5168    0 0.9432
#>   5475.84  1267.27    2600.28  82.27 0.5620    0 0.9442
#>  22782.99 20158.50   25366.63  81.47 0.5623    0 0.9749
#>  22782.99 20158.50   25366.63  81.47 0.5623    0 0.9749
```

The output includes three pseudo $`R^2`$ measures: `R2_cor` (squared
correlation), `R2_D` (explained deviance relative to a null model
without covariates, see Cameron & Trivedi, 2013), and `R2_CW`
(Cameron–Windmeijer, 1997, which accounts for the model’s variance
function).

**Note on `R2_D`**: The null model is the same specification (same
method, same gamma) but without `cov_alpha` / `cov_beta` — i.e., a
single $`\alpha`$ and single $`\beta`$. If your model already has no
covariates, then model = null and `R2_D = 0` by construction. `R2_D`
becomes informative when you compare models with covariates (e.g.,
`cov_alpha = ~ year + sex`) against the covariate-free baseline: it
tells you how much of the deviance is explained by group-varying
parameters beyond the basic power-law structure.

### Likelihood ratio test: Poisson vs NB

``` r

lr <- lrtest(fit_pois, fit_nb)
lr
#> Likelihood ratio test
#> ---------------------------------------- 
#> Model 1: POISSON   (logLik = -11380.65 )
#> Model 2: NB   (logLik = -2693.48 )
#> LR statistic: 17374.33 on 1 df
#> (Boundary-corrected: 0.5 * chi2_0 + 0.5 * chi2_1, Self & Liang 1987)
#> p-value: < 2.2e-16
```

When comparing Poisson (H0) against NB (H1), the dispersion parameter
theta lies on the boundary of its parameter space under H0
($`\theta \to \infty`$). The function automatically applies the Self &
Liang (1987) correction, yielding $`p = 0.5 \cdot \Pr(\chi^2_1 > LR)`$.

### Diagnostics

The four-panel diagnostic plot shows fitted vs observed, Anscombe
residuals, scale-location, and a normal Q-Q plot.

``` r

par(mfrow = c(2, 2))
plot(fit_pois, ask = FALSE)
```

![](workflow_files/figure-html/diagnostics-pois-1.png)

The rootogram compares observed and fitted count frequencies. Hanging
bars that touch the zero line indicate a good fit.

``` r

rootogram(fit_pois)
```

![](workflow_files/figure-html/rootogram-pois-1.png)

``` r

rootogram(fit_nb)
```

![](workflow_files/figure-html/rootogram-nb-1.png)

### Why Poisson as the main model?

Despite the NB often having a better AIC/BIC (reflecting genuine
overdispersion), the Poisson model is recommended as the primary
specification for several reasons:

1.  **Consistency.** Poisson PML is consistent for the conditional mean
    as long as $`E(m_i | X_i)`$ is correctly specified, regardless of
    the true variance function. The NB requires correct specification of
    both the mean and the variance.
2.  **Robustness.** HC-robust standard errors correct for overdispersion
    in the Poisson model without adding an extra parameter that must be
    estimated.
3.  **Stability.** The NB model’s joint estimation of gamma and theta
    can be numerically fragile when many observations are zero. The
    two-stage approach mitigates this but introduces a conditioning
    step.

The NB results serve as a useful robustness check on the population size
estimate.

### Exploratory log-log plots

The
[`plot_explore()`](https://ncn-foreigners.github.io/uncounted/reference/plot_explore.md)
function produces marginal log-log scatterplots that visualise the
power-law relationships.

``` r

plot_explore(fit_pois)
```

![](workflow_files/figure-html/explore-fit-1.png)

## Sensitivity Analysis

### Leave-one-out by country

LOO analysis refits the model dropping one country at a time and reports
the change in $`\hat{\xi}`$.

``` r

loo_pois <- loo(fit_pois, by = "country")
print(loo_pois)
#> Leave-one-out sensitivity analysis
#> Dropped by: country 
#> N iterations: 145 
#> Converged: 145 / 145 
#> 
#> Full model xi: 290597.4 
#> LOO xi range: 204772.5 to 900323.8 
#> 
#> Most influential (by |delta xi|):
#>         dropped      dxi pct_change
#>         Ukraine 609726.4     209.82
#>         Georgia -85824.9     -29.53
#>        Colombia  16718.5       5.75
#>      Uzbekistan  15906.1       5.47
#>         Vietnam  12464.1       4.29
#>           India  11926.5       4.10
#>  United Kingdom -10462.4      -3.60
#>          Russia   9177.5       3.16
#>         Belarus  -9061.5      -3.12
#>    Turkmenistan   7837.1       2.70
# #> Leave-one-out sensitivity analysis
# #> Dropped by: country
# #> ...
# #> Most influential (by |delta xi|):
# #>   dropped      dxi  pct_change
# #>   Ukraine  -12345.0       -8.50
# #>   ...
```

`dxi` holds the target fixed: it applies the refit’s coefficients to all
rows, so it measures how much the estimate moves because of that
country. The `dxi_subset` component also removes the country from the
total, so it includes the country’s own share of $`\hat\xi`$.

``` r

plot(loo_pois)
```

![](workflow_files/figure-html/loo-plot-1.png)

The coefficient stability summary shows how much alpha and beta shift:

``` r

summary(loo_pois)
#> Leave-one-out sensitivity analysis
#> Dropped by: country 
#> N iterations: 145 (converged: 145 )
#> 
#> === Coefficient stability ===
#>           full loo_mean   loo_sd  loo_min  loo_max
#> alpha 0.736646 0.737121 0.009246 0.702613 0.841615
#> beta  0.575046 0.576314 0.023531 0.500263 0.842034
#> 
#> === Xi stability ===
#> Full xi: 290597.4 
#> LOO mean: 294188.7 
#> LOO range: 204772.5 to 900323.8 
#> Max |%change|: 209.82 %
```

### Comparing LOO across models

``` r

loo_nb <- loo(fit_nb, by = "country")
comp_loo <- compare_loo(loo_pois, loo_nb, labels = c("Poisson", "NB"))
print(comp_loo)
#> LOO comparison: Poisson vs NB 
#> Dropped by: country 
#> Full xi -- Poisson : 290597 | NB : 1259529 
#> 
#> Top 15 most influential (by max |%change|):
#>  label             dxi_Poisson pct_Poisson dxi_NB     pct_NB max_abs_pct
#>  Ukraine           609726.42   209.82       301165.53  23.91 209.82     
#>  Georgia           -85824.94   -29.53      -103689.42  -8.23  29.53     
#>  unknown/stateless    745.69     0.26       240012.73  19.06  19.06     
#>  Belarus            -9061.47    -3.12       212862.67  16.90  16.90     
#>  Mauritius            708.16     0.24       144954.16  11.51  11.51     
#>  Burkina Faso         500.31     0.17       136739.46  10.86  10.86     
#>  Colombia           16718.47     5.75      -129251.63 -10.26  10.26     
#>  Guinea               727.23     0.25       114013.89   9.05   9.05     
#>  Congo - Kinshasa     320.97     0.11        88093.22   6.99   6.99     
#>  Uzbekistan         15906.09     5.47       -81297.90  -6.45   6.45     
#>  Tajikistan          7218.78     2.48        72733.23   5.77   5.77     
#>  Vietnam            12464.08     4.29       -66348.05  -5.27   5.27     
#>  Russia              9177.54     3.16       -63085.47  -5.01   5.01     
#>  Turkmenistan        7837.14     2.70       -52439.41  -4.16   4.16     
#>  India              11926.54     4.10       -28844.50  -2.29   4.10
plot(comp_loo, type = "scatter")
```

![](workflow_files/figure-html/compare-loo-1.png)

``` r

plot(comp_loo, type = "bar")
```

![](workflow_files/figure-html/compare-loo-2.png)

The scatter plot shows whether the same countries are influential under
both models. Points near the diagonal indicate agreement; points in
off-diagonal quadrants indicate divergent influence.

### Gamma profile

The
[`profile_gamma()`](https://ncn-foreigners.github.io/uncounted/reference/profile_gamma.md)
function refits the model across a grid of fixed gamma values and shows
how $`\hat{\xi}`$ and the log-likelihood depend on gamma. A flat xi
profile indicates robustness; a steep profile suggests sensitivity.

``` r

prof <- profile_gamma(fit_pois)
```

![](workflow_files/figure-html/gamma-profile-1.png)

The returned data frame can be inspected directly:

``` r

head(prof)
#>        gamma       xi    loglik
#> 1 0.00010000 217823.6 -11613.29
#> 2 0.02641053 409727.4 -12450.12
#> 3 0.05272105 388248.7 -13157.48
#> 4 0.07903158 366054.8 -13577.46
#> 5 0.10534211 349654.2 -13852.40
#> 6 0.13165263 337789.6 -14043.25
```

## Reporting Results

To prepare results for a paper, combine the analytical and bootstrap
estimates into a single summary table.

``` r

# Analytical results
ps_a <- popsize(fit_pois)

# Bootstrap results (assuming boot_pois was computed above)
ps_b <- boot_pois$popsize_full

results <- data.frame(
  Observed         = sum(ps_a$observed),
  Plugin           = round(sum(ps_a$estimate)),
  Plugin_BC        = round(sum(ps_a$estimate_bc)),
  Analytical_Lower = round(sum(ps_a$lower)),
  Analytical_Upper = round(sum(ps_a$upper)),
  Boot_Median      = round(sum(ps_b$boot_median)),
  Boot_Lower       = round(sum(ps_b$lower)),
  Boot_Upper       = round(sum(ps_b$upper))
)
results
#>   Observed Plugin Plugin_BC Analytical_Lower Analytical_Upper Boot_Median
#> 1    27105 290597    290498           127726           660703      289317
#>   Boot_Lower Boot_Upper
#> 1     183808     728846
```

Extract coefficients for reporting:

``` r

coefs <- coef(fit_pois)
se <- sqrt(diag(vcov(fit_pois)))
ci_lo <- coefs - 1.96 * se
ci_hi <- coefs + 1.96 * se

coef_table <- data.frame(
  Parameter = names(coefs),
  Estimate  = round(coefs, 4),
  SE        = round(se, 4),
  CI_lower  = round(ci_lo, 4),
  CI_upper  = round(ci_hi, 4),
  row.names = NULL
)
coef_table
#>   Parameter Estimate     SE CI_lower CI_upper
#> 1     alpha   0.7366 0.0403   0.6577   0.8156
#> 2      beta   0.5750 0.0919   0.3950   0.7551
```

The estimated gamma and model fit statistics:

``` r

cat("Gamma:", round(fit_pois$gamma, 6), "\n")
#> Gamma: 0.001831
cat("Log-likelihood:", round(fit_pois$loglik, 2), "\n")
#> Log-likelihood: -11380.65
cat("AIC:", round(AIC(fit_pois), 2), "\n")
#> AIC: 22767.3
cat("BIC:", round(BIC(fit_pois), 2), "\n")
#> BIC: 22782.99
cat("Deviance:", round(deviance(fit_pois), 2), "\n")
#> Deviance: 20158.5
```
