# Model Diagnostics and Comparison

This vignette demonstrates the diagnostic and model-comparison tools
available in the **uncounted** package. We cover goodness-of-fit
statistics, residual analysis, rootograms, exploratory log-log plots,
leave-one-out sensitivity analysis, and gamma profiling. Together these
tools help you decide which model specification is most appropriate and
whether the population-size estimates are robust.

## 1. Setup: fitting Poisson and Negative Binomial models

We use the `irregular_migration` dataset shipped with the package. Two
models are fitted—Poisson and Negative Binomial—with year and sex
covariates on alpha and year covariates on beta.

``` r
library(uncounted)
data(irregular_migration)

fit_po <- estimate_hidden_pop(
  data = irregular_migration,
  observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
  method = "poisson",
  cov_alpha = ~ factor(year) + sex, cov_beta = ~ factor(year),
  countries = ~ country_code
)

fit_nb <- estimate_hidden_pop(
  data = irregular_migration,
  observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
  method = "nb",
  cov_alpha = ~ factor(year) + sex, cov_beta = ~ factor(year),
  countries = ~ country_code
)
```

Quick summaries:

``` r
summary(fit_po)
#> Unauthorized population estimation
#> Method: POISSON | estimator: MLE | link_rho: power | vcov: HC3 
#> N obs: 1382 
#> Gamma: 0.00732 (estimated) 
#> Log-likelihood: -8217.25 
#> AIC: 16462.5  BIC: 16535.74 
#> Deviance: 13831.7 
#> 
#> Coefficients:
#>                          Estimate Std. Error z value  Pr(>|z|)    
#> alpha:(Intercept)       0.7886799  0.0903224  8.7318 < 2.2e-16 ***
#> alpha:factor(year)2020  0.0017375  0.0960815  0.0181 0.9855721    
#> alpha:factor(year)2021 -0.0170002  0.0872630 -0.1948 0.8455373    
#> alpha:factor(year)2022 -0.0934353  0.1092597 -0.8552 0.3924587    
#> alpha:factor(year)2023 -0.1134532  0.1773590 -0.6397 0.5223800    
#> alpha:factor(year)2024 -0.1689754  0.1613130 -1.0475 0.2948690    
#> alpha:sexMale           0.0421440  0.0458616  0.9189 0.3581279    
#> beta:(Intercept)        0.6756783  0.1950207  3.4646 0.0005309 ***
#> beta:factor(year)2020   0.1905329  0.2277210  0.8367 0.4027642    
#> beta:factor(year)2021   0.2347380  0.2126388  1.1039 0.2696242    
#> beta:factor(year)2022   0.0895142  0.2488296  0.3597 0.7190408    
#> beta:factor(year)2023  -0.0517138  0.3787723 -0.1365 0.8914023    
#> beta:factor(year)2024  -0.3223434  0.3189134 -1.0108 0.3121336    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>                       Observed Estimate Estimate (BC) CI lower CI upper
#> year=2019, sex=Female    1,535   21,068        21,037    3,398  130,226
#> year=2019, sex=Male      5,069   59,983        59,900   10,222  351,029
#> year=2020, sex=Female      698   23,697        23,632    6,016   92,828
#> year=2020, sex=Male      2,700   67,015        66,823   22,884  195,132
#> year=2021, sex=Female      483   22,806        22,741    7,274   71,092
#> year=2021, sex=Male      2,622   65,302        65,107   32,364  130,974
#> year=2022, sex=Female      317   14,115        14,079    2,476   80,059
#> year=2022, sex=Male      2,632   32,960        32,882    7,369  146,728
#> year=2023, sex=Female      523   12,063        12,043      533  272,157
#> year=2023, sex=Male      3,839   28,877        28,833    1,219  681,962
#> year=2024, sex=Female      956    7,106         7,100      575   87,640
#> year=2024, sex=Male      5,731   17,487        17,473    1,175  259,837
summary(fit_nb)
#> Unauthorized population estimation
#> Method: NB | estimator: MLE | link_rho: power | vcov: HC1 
#> N obs: 1382 
#> Gamma: 0.020902 (estimated) 
#> Theta (NB dispersion): 1.3552 
#> Log-likelihood: -2621.12 
#> AIC: 5272.23  BIC: 5350.7 
#> Deviance: 1275.92 
#> 
#> Coefficients:
#>                         Estimate Std. Error z value  Pr(>|z|)    
#> alpha:(Intercept)       0.815858   0.034303 23.7838 < 2.2e-16 ***
#> alpha:factor(year)2020 -0.041962   0.057019 -0.7359 0.4617745    
#> alpha:factor(year)2021  0.075292   0.059612  1.2630 0.2065778    
#> alpha:factor(year)2022 -0.030061   0.073391 -0.4096 0.6821026    
#> alpha:factor(year)2023  0.015592   0.055919  0.2788 0.7803731    
#> alpha:factor(year)2024  0.069599   0.047748  1.4576 0.1449455    
#> alpha:sexMale           0.044285   0.013138  3.3708 0.0007495 ***
#> beta:(Intercept)        0.873376   0.092225  9.4700 < 2.2e-16 ***
#> beta:factor(year)2020   0.140459   0.116675  1.2038 0.2286488    
#> beta:factor(year)2021   0.467946   0.121527  3.8505 0.0001179 ***
#> beta:factor(year)2022   0.185009   0.151482  1.2213 0.2219630    
#> beta:factor(year)2023   0.183925   0.110426  1.6656 0.0957949 .  
#> beta:factor(year)2024   0.132684   0.092752  1.4305 0.1525665    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>                       Observed Estimate Estimate (BC) CI lower CI upper
#> year=2019, sex=Female    1,535   27,938        25,619   12,674   51,788
#> year=2019, sex=Male      5,069   82,589        75,347   35,613  159,414
#> year=2020, sex=Female      698   19,965        18,070    6,727   48,539
#> year=2020, sex=Male      2,700   57,322        51,551   18,770  141,581
#> year=2021, sex=Female      483   82,785        72,665   24,233  217,890
#> year=2021, sex=Male      2,622  256,825       223,774   71,832  697,112
#> year=2022, sex=Female      317   38,077        34,486    7,929  149,989
#> year=2022, sex=Male      2,632   89,763        81,957   20,874  321,787
#> year=2023, sex=Female      523   66,694        60,766   20,291  181,981
#> year=2023, sex=Male      3,839  159,171       146,421   52,747  406,450
#> year=2024, sex=Female      956  120,792       111,067   46,511  265,228
#> year=2024, sex=Male      5,731  302,634       280,353  121,017  649,474
```

## 2. Model comparison

### Comparison table

[`compare_models()`](https://ncn-foreigners.github.io/uncounted/reference/compare_models.md)
produces a side-by-side table of log-likelihood, AIC, BIC, deviance,
Pearson chi-squared, RMSE, and three pseudo $R^{2}$ measures:

- **`R2_cor`**: ${cor}\left( m,\widehat{\mu} \right)^{2}$ — squared
  correlation between observed and fitted values.
- **`R2_D`**: Explained deviance
  $1 - D\left( \text{model} \right)/D\left( \text{null} \right)$, where
  the null model is the same specification without covariates (single
  $\alpha$, single $\beta$). Measures how much covariates improve the
  fit beyond the baseline power-law structure.
- **`R2_CW`**: Cameron–Windmeijer (1996) pseudo $R^{2}$, which uses the
  model-implied variance function $V(\mu)$ — specifically designed for
  count data regression.

Models are sorted by AIC (the default). Note that `R2_D` compares the
fitted model against a null model (same method, no covariates). If the
model has no `cov_alpha` / `cov_beta`, then `R2_D = 0` by construction —
it only becomes informative when covariates are present.

``` r
comp <- compare_models(Poisson = fit_po, NB = fit_nb, sort_by = "AIC")
comp
#> Model comparison
#> ------------------------------------------------------------ 
#>    Model  Method Estimator  Link Constrained n_par   logLik      AIC      BIC
#>       NB      NB       MLE power       FALSE    15 -2621.12  5272.23  5350.70
#>  Poisson POISSON       MLE power       FALSE    14 -8217.25 16462.50 16535.74
#>  Deviance Pearson_X2   RMSE R2_cor    R2_D  R2_CW
#>   1275.92    3090.18 206.02 0.4056 -0.0029 0.9526
#>  13831.70   17499.93  52.53 0.8280  0.3139 0.9827
```

Lower AIC and BIC values indicate a better trade-off between fit and
complexity. The Pearson chi-squared statistic divided by its degrees of
freedom (approximately `n_obs - n_par`) gives a quick dispersion check:
values much greater than 1 suggest overdispersion, which favours the NB
model over Poisson.

### Likelihood ratio test

Because the Poisson model is nested within the Negative Binomial (the
Poisson is the NB with theta -\> infinity), we can use a likelihood
ratio test.
[`lrtest()`](https://ncn-foreigners.github.io/uncounted/reference/lrtest.md)
applies the boundary correction of Self & Liang (1987) automatically
when comparing Poisson against NB, because the dispersion parameter lies
on the boundary of the parameter space under H0.

``` r
lrtest(fit_po, fit_nb)
#> Likelihood ratio test
#> ---------------------------------------- 
#> Model 1: POISSON   (logLik = -8217.25 )
#> Model 2: NB   (logLik = -2621.12 )
#> LR statistic: 11192.27 on 1 df
#> (Boundary-corrected: 0.5 * P(chi2 > LR), Self & Liang 1987)
#> p-value: < 2.2e-16
```

A small p-value indicates significant overdispersion, favouring the NB
model.

## 3. Residual diagnostics

### Residual types

The [`residuals()`](https://rdrr.io/r/stats/residuals.html) method
supports four types:

- **response**: raw residuals, $m_{i} - {\widehat{\mu}}_{i}$.
- **pearson**: standardised by the variance function,
  $\left( m_{i} - {\widehat{\mu}}_{i} \right)/\sqrt{V\left( {\widehat{\mu}}_{i} \right)}$.
- **deviance**: signed square root of individual deviance contributions.
- **anscombe**: variance-stabilising residuals (approximately standard
  normal for a well-specified model).

``` r
r_resp <- residuals(fit_nb, type = "response")
r_pear <- residuals(fit_nb, type = "pearson")
r_dev  <- residuals(fit_nb, type = "deviance")
r_ansc <- residuals(fit_nb, type = "anscombe")

summary(r_pear)
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -1.15132 -0.59481 -0.28876  0.03627  0.11109 22.91161
summary(r_ansc)
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -24.1664  -1.0757  -0.4675  -0.3348   0.1767   9.4214
```

### Four-panel diagnostic plot

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a fitted
model produces the Zhang (2008) four-panel display:

1.  **Fitted vs Observed (sqrt scale)** – points should cluster around
    the 45-degree line.
2.  **Anscombe residuals vs fitted** – the running mean should be flat
    near zero (no systematic misfit).
3.  **Scale-Location** – absolute Anscombe residuals vs fitted; an
    increasing trend signals under-modelled variance.
4.  **Normal Q-Q** – Anscombe residuals should follow the reference line
    if the distributional assumption holds.

``` r
op <- par(mfrow = c(2, 2))
plot(fit_nb, ask = FALSE)
```

![](diagnostics_files/figure-html/plot-diagnostics-1.png)

``` r
par(op)
```

**What to look for:**

- Curvature in Panel 1 suggests the mean function is misspecified.
- A non-flat smoother in Panel 2 indicates systematic bias at certain
  fitted-value ranges.
- A fan shape in Panel 3 points to overdispersion or heteroscedasticity.
- Heavy tails or S-shapes in Panel 4 indicate departure from the assumed
  distribution.

Compare with the Poisson fit to see whether switching to NB resolves any
patterns:

``` r
op <- par(mfrow = c(2, 2))
plot(fit_po, ask = FALSE)
```

![](diagnostics_files/figure-html/plot-diagnostics-po-1.png)

``` r
par(op)
```

## 4. Rootograms

Rootograms (Kleiber & Zeileis, 2016) compare the observed and fitted
count-frequency distributions directly. Three display styles are
available:

- **hanging** (default): observed bars are hung from the fitted curve.
  If the model fits well, the bottom of each bar touches zero.
- **standing**: observed bars stand on the x-axis with the fitted curve
  overlaid. Harder to judge discrepancies visually.
- **suspended**: plots the difference
  $\sqrt{f_{\text{obs}}} - \sqrt{f_{\text{exp}}}$ directly. Bars
  crossing zero indicate over- or under-prediction at that count.

``` r
rootogram(fit_nb, style = "hanging")
```

![](diagnostics_files/figure-html/rootogram-nb-1.png)

``` r
rootogram(fit_po, style = "hanging")
```

![](diagnostics_files/figure-html/rootogram-po-1.png)

Compare the two: if the Poisson rootogram shows bars consistently
hanging below the zero line at small counts and above at moderate
counts, this is a classic overdispersion signature that the NB model
should correct.

``` r
rootogram(fit_nb, style = "standing")
```

![](diagnostics_files/figure-html/rootogram-styles-1.png)

``` r
rootogram(fit_nb, style = "suspended")
```

![](diagnostics_files/figure-html/rootogram-styles-2.png)

## 5. Exploratory log-log plots

Before or after fitting,
[`plot_explore()`](https://ncn-foreigners.github.io/uncounted/reference/plot_explore.md)
produces two log-log scatter plots that reveal the marginal
relationships underlying the power-law model:

1.  **log(m/N) vs log(N)**: the slope approximates $\alpha - 1$. A slope
    near zero means the apprehension rate is roughly independent of
    population size ($\alpha \approx 1$). A negative slope indicates
    decreasing returns to scale ($\alpha < 1$).
2.  **log(m/N) vs log(n/N)**: the slope approximates $\beta$. A positive
    slope confirms that the auxiliary rate is predictive of the
    apprehension rate.

Observations with $m = 0$ or $n = 0$ are excluded (undefined on the log
scale).

``` r
plot_explore(fit_nb)
```

![](diagnostics_files/figure-html/plot-explore-1.png)

These plots are exploratory visual checks. The OLS lines overlaid are
rough guides—the formal model accounts for covariates and the gamma
offset.

## 6. Leave-one-out sensitivity analysis

LOO refits the model repeatedly, dropping one observation (or one
country) at a time, to assess how much each unit influences the total
population-size estimate.

### LOO by country

Dropping all rows for a country measures its structural contribution.
This is computationally intensive (one refit per country).

``` r
loo_po_ctry <- loo(fit_po, by = "country", verbose = TRUE)
loo_nb_ctry <- loo(fit_nb, by = "country", verbose = TRUE)
```

``` r
print(loo_nb_ctry)
summary(loo_nb_ctry)
```

The [`print()`](https://rdrr.io/r/base/print.html) method shows the top
10 most influential countries ranked by $|\Delta\xi|$, and
[`summary()`](https://rdrr.io/r/base/summary.html) reports coefficient
and xi stability across all LOO refits.

``` r
plot(loo_nb_ctry, type = "xi")
plot(loo_nb_ctry, type = "coef")
```

The xi bar plot shows which countries pull the total estimate up (red
bars, dropping the country decreases the estimate) or down (blue bars).

### LOO by observation

Dropping individual rows identifies specific data points that drive the
estimate. Useful for outlier detection.

``` r
loo_nb_obs <- loo(fit_nb, by = "obs", verbose = TRUE)
print(loo_nb_obs)
plot(loo_nb_obs)
```

### Comparing LOO across models

[`compare_loo()`](https://ncn-foreigners.github.io/uncounted/reference/compare_loo.md)
brings together LOO results from two models and produces a scatter plot
where each axis shows the percentage change in xi when dropping a given
unit:

``` r
comp_loo <- compare_loo(
  loo_po_ctry, loo_nb_ctry,
  labels = c("Poisson", "NB")
)
print(comp_loo)
```

**Scatter plot interpretation (4 quadrants):**

``` r
plot(comp_loo, type = "scatter")
```

- **Top-right (+, +):** dropping this unit increases the estimate under
  both models. The unit was pulling the estimate down in both.
- **Bottom-left (-, -):** dropping this unit decreases the estimate
  under both models. The unit was pulling the estimate up in both.
- **Off-diagonal (top-left or bottom-right):** divergent influence. The
  unit affects the two models in opposite directions—a sign that model
  choice matters for that unit.
- **Near the origin:** non-influential under both models.

Points near the diagonal indicate consistent influence; points far from
the diagonal warrant investigation.

``` r
plot(comp_loo, type = "bar")
```

The side-by-side bar plot ranks countries by maximum absolute percentage
change across both models, making it easy to spot the most consequential
units.

## 7. Gamma profile

[`profile_gamma()`](https://ncn-foreigners.github.io/uncounted/reference/profile_gamma.md)
refits the model across a grid of fixed gamma values and plots two
curves:

1.  **$\widehat{\xi}(\gamma)$**: how the total population-size estimate
    changes with gamma.
2.  **$\ell(\gamma)$**: the profile log-likelihood as a function of
    gamma.

``` r
prof <- profile_gamma(fit_nb, gamma_grid = seq(1e-4, 0.5, length.out = 10))
```

![](diagnostics_files/figure-html/profile-gamma-1.png)

**Interpretation:**

- A **sharp peak** in the log-likelihood curve means gamma is
  well-identified by the data. The MLE is precise and the
  population-size estimate is robust to small perturbations in gamma.
- A **flat** log-likelihood curve means gamma is weakly identified. The
  data cannot distinguish between a range of gamma values, and the
  population-size estimate may be sensitive to the gamma assumption. In
  this case, reporting results for several gamma values (or fixing gamma
  based on external information) is advisable.

If the original model estimated gamma, its point estimate is marked with
a dashed red vertical line. The xi profile shows whether the
population-size estimate is stable (flat) or sensitive (steep) across
the gamma range.

The profiling results are returned invisibly as a data frame:

``` r
prof <- profile_gamma(fit_nb, plot = FALSE)
head(prof)
#>        gamma        xi    loglik
#> 1 0.00010000  105417.3 -2719.707
#> 2 0.02641053 1423540.0 -2621.493
#> 3 0.05272105 1744889.5 -2626.657
#> 4 0.07903158 1894473.3 -2631.701
#> 5 0.10534211 1979041.5 -2635.659
#> 6 0.13165263 2037460.9 -2638.698
```

## References

- Beręsewicz, M., & Pawlukiewicz, K. (2020). Estimation of the number of
  irregular foreigners in Poland using non-linear count regression
  models. arXiv preprint arXiv:2008.09407.
- Kleiber, C. and Zeileis, A. (2016). Visualizing Count Data Regressions
  Using Rootograms. *The American Statistician*, 70(3), 296–303.
- McCullagh, P. and Nelder, J. A. (1989). *Generalized Linear Models*,
  2nd ed. Chapman & Hall.
- Self, S. G. and Liang, K.-Y. (1987). Asymptotic properties of maximum
  likelihood estimators and likelihood ratio tests under nonstandard
  conditions. *JASA*, 82(398), 605–610.
- Zhang, L.-C. (2008). *Developing methods for determining the number of
  unauthorized foreigners in Norway* (Documents 2008/11). Statistics
  Norway.
  <https://www.ssb.no/a/english/publikasjoner/pdf/doc_200811_en/doc_200811_en.pdf>
- Beręsewicz, M., & Pawlukiewicz, K. (2020). Estimation of the number of
  irregular foreigners in Poland using non-linear count regression
  models. arXiv preprint arXiv:2008.09407.
