# Bootstrap Confidence Intervals for Population Size

Computes bootstrap confidence intervals for the estimated population
size using the fractional weighted bootstrap (FWB) of Xu et al. (2020).
Supports cluster bootstrap to account for within-country correlation
across years and sex.

## Usage

``` r
bootstrap_popsize(
  object,
  R = 199,
  cluster = NULL,
  level = 0.95,
  ci_type = c("perc", "bc", "cheap"),
  point_estimate = c("median", "plugin", "mean"),
  seed = NULL,
  verbose = TRUE,
  by = NULL,
  total = FALSE
)
```

## Arguments

- object:

  An `"uncounted"` object (fitted model).

- R:

  Number of bootstrap replications (default 199).

- cluster:

  One-sided formula for cluster variable (e.g., `~ country_code`). When
  provided, FWB weights are generated at the cluster level (all
  observations in a cluster get the same weight). When `NULL`,
  observation-level FWB is used.

- level:

  Confidence level (default 0.95).

- ci_type:

  Type of bootstrap CI: `"perc"` (percentile), `"bc"` (bias-corrected
  percentile), or `"cheap"` (Lam, 2022), which keeps nominal coverage
  with very few replications because it needs only their spread around
  the plug-in estimate: \\\exp\\\log\hat\xi \pm t\_{R, 1-a/2} S\\\\ with
  \\S^2 = R^{-1}\sum_r (\log\hat\xi^\*\_r - \log\hat\xi)^2\\, computed
  on the log scale so the interval stays positive. Useful when refits
  are slow; `R` can be as small as 1, and `"perc"`/`"bc"` need at least
  10 converged replications. Default `"perc"`.

- point_estimate:

  Which point estimate to report: `"median"` (default, recommended,
  except with `ci_type = "cheap"` where the default is `"plugin"`) –
  bootstrap median, robust to skewness; `"plugin"` – plug-in \\\hat{\xi}
  = \sum N^{\hat{\alpha}}\\; `"mean"` – bootstrap mean.

- seed:

  Random seed for reproducibility.

- verbose:

  Print progress bar? Default TRUE.

- by:

  Optional formula for stratified population size estimation (e.g.,
  `~ year`, `~ country`). When provided, bootstrap population sizes are
  computed per `by`-group instead of per `cov_alpha` group. See
  [`popsize`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
  for details.

- total:

  Logical; if `TRUE` and multiple groups exist, include a total row
  computed by summing per-replicate group estimates before taking
  quantiles. Default `FALSE`.

## Value

An object of class `"uncounted_boot"` with components:

- t:

  Matrix (`R` x `n_groups`): bootstrap population size estimates per
  group. Each row is one replicate.

- t0:

  Numeric vector: plug-in point estimates from original fit.

- t0_bc:

  Numeric vector: analytical bias-corrected point estimates from
  [`popsize()`](https://ncn-foreigners.github.io/uncounted/reference/popsize.md)
  (multiplicative lognormal for unconstrained fits, Taylor approximation
  for constrained fits).

- popsize:

  Data frame with columns: `group`, `estimate`, `lower`, `upper`, where
  `estimate` uses the chosen `point_estimate` type.

- popsize_full:

  Data frame with all point estimate types: `plugin`, `plugin_bc`,
  `boot_median`, `boot_mean`, `lower`, `upper`.

- total:

  List with total across groups (`plugin`, `plugin_bc`, `median`,
  `mean`, `lower`, `upper`), or `NULL` if only one group.

- R:

  Number of replications requested.

- ci_type:

  CI type used (`"perc"`, `"bc"` or `"cheap"`).

- point_estimate:

  Point estimate type used.

- level:

  Confidence level.

- n_converged:

  Number of bootstrap replications where the model converged.

- cluster:

  Logical; whether cluster bootstrap was used.

## Details

**Fractional weighted bootstrap (FWB).** Instead of resampling rows, FWB
generates random non-negative weights from a Dirichlet distribution and
refits the model with weighted likelihood (Xu et al., 2020). This avoids
duplicate-row issues with discrete data and is better suited to GLMs.
The fwb package provides the implementation.

**Cluster bootstrap.** When `cluster` is specified, all observations
within the same cluster (e.g., country) receive the same FWB weight.
This preserves within-cluster correlation structure (e.g., the same
country observed across multiple years and sex categories). Cluster
bootstrap is recommended whenever observations are not independent.

**Why bootstrap median as point estimate.** The transformation \\\xi =
\sum N_i^{\alpha}\\ is convex in \\\alpha\\ for \\N \> 1\\, so by
Jensen's inequality \\E\[\hat{\xi}\] \geq \xi\\. The bootstrap
distribution of \\\hat{\xi}\\ is therefore typically right-skewed, and
the bootstrap mean overestimates \\\xi\\. The bootstrap median is a more
robust central-tendency summary that is less sensitive to this skewness,
and is the default point estimate.

**Percentile CI (`ci_type = "perc"`).** The interval is
\\\[\hat{\xi}^{\*}\_{\alpha/2}, \\ \hat{\xi}^{\*}\_{1-\alpha/2}\]\\,
i.e., the \\\alpha/2\\ and \\1-\alpha/2\\ quantiles of the bootstrap
distribution.

**Bias-corrected percentile CI (`ci_type = "bc"`).** Adjusts the
quantile probabilities for median bias. Let \\z_0 =
\Phi^{-1}(\mathrm{Pr}(\hat{\xi}^\* \< \hat{\xi}\_0))\\. The adjusted
quantile probabilities are \\p_L = \Phi(2z_0 + z\_{\alpha/2})\\ and
\\p_U = \Phi(2z_0 - z\_{\alpha/2})\\, yielding
\\\[\hat{\xi}^{\*}\_{p_L}, \\ \hat{\xi}^{\*}\_{p_U}\]\\.

**Total across groups.** When multiple groups exist, the total is
computed by summing per-replicate group estimates before taking
quantiles. This correctly accounts for between-group correlation within
each bootstrap replicate, unlike summing marginal quantiles.

## References

Xu, L., Gotwalt, C., Hong, Y., King, C. B., & Meeker, W. Q. (2020).
Applications of the fractional-random-weight bootstrap. *The American
Statistician*, 74(4), 345–358.
[doi:10.1080/00031305.2020.1731599](https://doi.org/10.1080/00031305.2020.1731599)

Greifer, N. (2026). fwb: An R package for the fractional weighted
bootstrap. *Journal of Open Source Software*, 11(124), 8880.
[doi:10.21105/joss.08880](https://doi.org/10.21105/joss.08880)

Lam, H. (2022). A cheap bootstrap method for fast inference. *arXiv
preprint*.
[doi:10.48550/arXiv.2202.00090](https://doi.org/10.48550/arXiv.2202.00090)

## Examples

``` r
# Simulate synthetic data
set.seed(42)
n_obs <- 15
sim_data <- data.frame(
  country = rep(paste0("C", 1:5), each = 3),
  year    = rep(2018:2020, 5),
  N       = rpois(n_obs, lambda = 500000),
  n       = rpois(n_obs, lambda = 1000),
  m       = rpois(n_obs, lambda = 50)
)

fit <- estimate_hidden_pop(
  data = sim_data, observed = ~m, auxiliary = ~n,
  reference_pop = ~N, method = "poisson",
  countries = ~country
)
#> Warning: Some alpha values < 0 (min = -0.046). Consider using constrained = TRUE.
#> Estimated gamma is at its lower bound (1e-10); the optimum may lie outside gamma_bounds = c(1e-10, 0.5).

# \donttest{
# Cluster bootstrap (recommended when data has repeated country obs)
boot <- bootstrap_popsize(fit, R = 99, cluster = ~country, seed = 123)
#> Warning: Some alpha values < 0 (min = -0.046). Consider using constrained = TRUE.
#> 
  |                                                  | 0 % ~calculating  
  |+                                                 | 1 % ~01s          
  |++                                                | 2 % ~01s          
#> Warning: Some alpha values < 0 (min = -0.204). Consider using constrained = TRUE.
#> 
  |++                                                | 3 % ~01s          
#> Warning: Some alpha values < 0 (min = -0.437). Consider using constrained = TRUE.
#> 
  |+++                                               | 4 % ~01s          
  |+++                                               | 5 % ~01s          
  |++++                                              | 6 % ~01s          
#> Warning: Some alpha values < 0 (min = -0.183). Consider using constrained = TRUE.
#> 
  |++++                                              | 7 % ~01s          
  |+++++                                             | 8 % ~01s          
  |+++++                                             | 9 % ~00s          
#> Warning: Some alpha values < 0 (min = -0.083). Consider using constrained = TRUE.
#> 
  |++++++                                            | 10% ~01s          
#> Warning: Some alpha values < 0 (min = -0.445). Consider using constrained = TRUE.
#> 
  |++++++                                            | 11% ~01s          
  |+++++++                                           | 12% ~01s          
  |+++++++                                           | 13% ~01s          
#> Warning: Some alpha values < 0 (min = -0.417). Consider using constrained = TRUE.
#> 
  |++++++++                                          | 14% ~00s          
#> Warning: Some alpha values < 0 (min = -0.069). Consider using constrained = TRUE.
#> 
  |++++++++                                          | 15% ~01s          
#> Warning: Some alpha values < 0 (min = -0.184). Consider using constrained = TRUE.
#> 
  |+++++++++                                         | 16% ~00s          
  |+++++++++                                         | 17% ~00s          
  |++++++++++                                        | 18% ~01s          
  |++++++++++                                        | 19% ~01s          
  |+++++++++++                                       | 20% ~00s          
  |+++++++++++                                       | 21% ~00s          
#> Warning: Some alpha values < 0 (min = -0.048). Consider using constrained = TRUE.
#> 
  |++++++++++++                                      | 22% ~00s          
#> Warning: Some alpha values < 0 (min = -0.176). Consider using constrained = TRUE.
#> 
  |++++++++++++                                      | 23% ~00s          
#> Warning: Some alpha values < 0 (min = -0.331). Consider using constrained = TRUE.
#> 
  |+++++++++++++                                     | 24% ~00s          
  |+++++++++++++                                     | 25% ~00s          
#> Warning: Some alpha values < 0 (min = -0.093). Consider using constrained = TRUE.
#> 
  |++++++++++++++                                    | 26% ~00s          
  |++++++++++++++                                    | 27% ~00s          
  |+++++++++++++++                                   | 28% ~00s          
  |+++++++++++++++                                   | 29% ~00s          
#> Warning: Some alpha values < 0 (min = -0.083). Consider using constrained = TRUE.
#> 
  |++++++++++++++++                                  | 30% ~00s          
#> Warning: Some alpha values < 0 (min = -0.214). Consider using constrained = TRUE.
#> 
  |++++++++++++++++                                  | 31% ~00s          
  |+++++++++++++++++                                 | 32% ~00s          
#> Warning: Some alpha values < 0 (min = -0.287). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++                                 | 33% ~00s          
#> Warning: Some alpha values < 0 (min = -0.388). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++                                | 34% ~00s          
  |++++++++++++++++++                                | 35% ~00s          
  |+++++++++++++++++++                               | 36% ~00s          
#> Warning: Some alpha values < 0 (min = -0.328). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++                               | 37% ~00s          
  |++++++++++++++++++++                              | 38% ~00s          
  |++++++++++++++++++++                              | 39% ~00s          
#> Warning: Some alpha values < 0 (min = -0.11). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++                             | 40% ~00s          
#> Warning: Some alpha values < 0 (min = -0.155). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++                             | 41% ~00s          
#> Warning: Some alpha values < 0 (min = -0.431). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++                            | 42% ~00s          
  |++++++++++++++++++++++                            | 43% ~00s          
  |+++++++++++++++++++++++                           | 44% ~00s          
  |+++++++++++++++++++++++                           | 45% ~00s          
  |++++++++++++++++++++++++                          | 46% ~00s          
#> Warning: Some alpha values < 0 (min = -0.418). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++                          | 47% ~00s          
#> Warning: Some alpha values < 0 (min = -0.116). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++                         | 48% ~00s          
#> Warning: Some alpha values < 0 (min = -0.061). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++                         | 49% ~00s          
#> Warning: Some alpha values < 0 (min = -0.312). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++                        | 51% ~00s          
  |++++++++++++++++++++++++++                        | 52% ~00s          
#> Warning: Some alpha values < 0 (min = -0.257). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++                       | 53% ~00s          
  |+++++++++++++++++++++++++++                       | 54% ~00s          
  |++++++++++++++++++++++++++++                      | 55% ~00s          
#> Warning: Some alpha values < 0 (min = -0.112). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++                      | 56% ~00s          
#> Warning: Some alpha values > 1 (max = 1.935). Population size estimates may be unreliable. Consider using constrained = TRUE or simplifying cov_alpha.
#> 
  |+++++++++++++++++++++++++++++                     | 57% ~00s          
#> Warning: Some alpha values < 0 (min = -0.412). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++                     | 58% ~00s          
  |++++++++++++++++++++++++++++++                    | 59% ~00s          
  |++++++++++++++++++++++++++++++                    | 60% ~00s          
#> Warning: Some alpha values < 0 (min = -0.286). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++                   | 61% ~00s          
  |+++++++++++++++++++++++++++++++                   | 62% ~00s          
  |++++++++++++++++++++++++++++++++                  | 63% ~00s          
#> Warning: Some alpha values < 0 (min = -0.064). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++                  | 64% ~00s          
#> Warning: Some alpha values < 0 (min = -0.272). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++                 | 65% ~00s          
#> Warning: Some alpha values < 0 (min = -0.599). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++                 | 66% ~00s          
#> Warning: Some alpha values < 0 (min = -0.701). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++                | 67% ~00s          
  |++++++++++++++++++++++++++++++++++                | 68% ~00s          
  |+++++++++++++++++++++++++++++++++++               | 69% ~00s          
#> Warning: Some alpha values < 0 (min = -0.103). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++               | 70% ~00s          
  |++++++++++++++++++++++++++++++++++++              | 71% ~00s          
#> Warning: Some alpha values < 0 (min = -0.114). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++              | 72% ~00s          
  |+++++++++++++++++++++++++++++++++++++             | 73% ~00s          
  |+++++++++++++++++++++++++++++++++++++             | 74% ~00s          
#> Warning: Some alpha values < 0 (min = -0.076). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++            | 75% ~00s          
  |++++++++++++++++++++++++++++++++++++++            | 76% ~00s          
#> Warning: Some alpha values < 0 (min = -0.429). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++++++           | 77% ~00s          
  |+++++++++++++++++++++++++++++++++++++++           | 78% ~00s          
#> Warning: Some alpha values < 0 (min = -0.037). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++++          | 79% ~00s          
  |++++++++++++++++++++++++++++++++++++++++          | 80% ~00s          
#> Warning: Some alpha values < 0 (min = -0.572). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++++++++         | 81% ~00s          
  |+++++++++++++++++++++++++++++++++++++++++         | 82% ~00s          
  |++++++++++++++++++++++++++++++++++++++++++        | 83% ~00s          
#> Warning: Some alpha values < 0 (min = -0.281). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++++++        | 84% ~00s          
  |+++++++++++++++++++++++++++++++++++++++++++       | 85% ~00s          
#> Warning: Some alpha values < 0 (min = -0.156). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++++++++++       | 86% ~00s          
  |++++++++++++++++++++++++++++++++++++++++++++      | 87% ~00s          
#> Warning: Some alpha values < 0 (min = -0.013). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++++++++      | 88% ~00s          
#> Warning: Some alpha values < 0 (min = -0.034). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++++++++++++     | 89% ~00s          
  |+++++++++++++++++++++++++++++++++++++++++++++     | 90% ~00s          
  |++++++++++++++++++++++++++++++++++++++++++++++    | 91% ~00s          
#> Warning: Some alpha values < 0 (min = -0.312). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++++++++++    | 92% ~00s          
#> Warning: Some alpha values < 0 (min = -0.064). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++++++++++++++   | 93% ~00s          
#> Warning: Some alpha values < 0 (min = -0.248). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++++++++++++++   | 94% ~00s          
#> Warning: Some alpha values < 0 (min = -0.319). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++++++++++++  | 95% ~00s          
#> Warning: Some alpha values < 0 (min = -0.217). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++++++++++++  | 96% ~00s          
#> Warning: Some alpha values < 0 (min = -0.504). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++++++++++++++++ | 97% ~00s          
  |+++++++++++++++++++++++++++++++++++++++++++++++++ | 98% ~00s          
#> Warning: Some alpha values < 0 (min = -0.2). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 99% ~00s          
#> Warning: Some alpha values < 0 (min = -0.021). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  
boot                # prints table with median as point estimate
#> Bootstrap population size estimation
#> R = 99 | CI type: perc | Point estimate: median | Converged: 99 / 99 
#> Cluster bootstrap
#> 95% CI
#> 
#>   Point estimate: bootstrap median (recommended) | CI: bootstrap percentile
#>       Plugin Plugin (BC) Boot median      Boot mean CI lower CI upper
#> (all)      8           0          13 16,164,990,743        0   26,038
boot$popsize        # data.frame for further use
#>   group estimate     lower    upper
#> 1 (all) 12.70394 0.0135202 26038.06
boot$popsize_full   # all estimate types side by side
#>   group  plugin    plugin_bc boot_median   boot_mean     lower    upper
#> 1 (all) 8.22019 5.988907e-11    12.70394 16164990743 0.0135202 26038.06

# Bias-corrected percentile CI
boot_bc <- bootstrap_popsize(fit, R = 99, ci_type = "bc", seed = 123)
#> Warning: Some alpha values < 0 (min = -0.046). Consider using constrained = TRUE.
#> 
  |                                                  | 0 % ~calculating  
  |+                                                 | 1 % ~01s          
  |++                                                | 2 % ~01s          
#> Warning: Some alpha values < 0 (min = -0.315). Consider using constrained = TRUE.
#> 
  |++                                                | 3 % ~01s          
  |+++                                               | 4 % ~01s          
#> Warning: Some alpha values < 0 (min = -0.158). Consider using constrained = TRUE.
#> 
  |+++                                               | 5 % ~01s          
#> Warning: Some alpha values < 0 (min = -0.583). Consider using constrained = TRUE.
#> 
  |++++                                              | 6 % ~01s          
#> Warning: Some alpha values < 0 (min = -0.104). Consider using constrained = TRUE.
#> 
  |++++                                              | 7 % ~01s          
#> Warning: Some alpha values < 0 (min = -0.722). Consider using constrained = TRUE.
#> 
  |+++++                                             | 8 % ~01s          
#> Warning: Some alpha values < 0 (min = -0.102). Consider using constrained = TRUE.
#> 
  |+++++                                             | 9 % ~01s          
#> Warning: Some alpha values < 0 (min = -0.045). Consider using constrained = TRUE.
#> 
  |++++++                                            | 10% ~01s          
  |++++++                                            | 11% ~01s          
  |+++++++                                           | 12% ~01s          
  |+++++++                                           | 13% ~01s          
#> Warning: Some alpha values < 0 (min = -0.098). Consider using constrained = TRUE.
#> 
  |++++++++                                          | 14% ~01s          
#> Warning: Some alpha values < 0 (min = -0.5). Consider using constrained = TRUE.
#> 
  |++++++++                                          | 15% ~01s          
#> Warning: Some alpha values < 0 (min = -0.343). Consider using constrained = TRUE.
#> 
  |+++++++++                                         | 16% ~01s          
#> Warning: Some alpha values < 0 (min = -0.322). Consider using constrained = TRUE.
#> 
  |+++++++++                                         | 17% ~01s          
  |++++++++++                                        | 18% ~01s          
#> Warning: Some alpha values < 0 (min = -0.18). Consider using constrained = TRUE.
#> 
  |++++++++++                                        | 19% ~01s          
#> Warning: Some alpha values < 0 (min = -0.031). Consider using constrained = TRUE.
#> 
  |+++++++++++                                       | 20% ~01s          
#> Warning: Some alpha values < 0 (min = -0.197). Consider using constrained = TRUE.
#> 
  |+++++++++++                                       | 21% ~01s          
  |++++++++++++                                      | 22% ~01s          
#> Warning: Some alpha values < 0 (min = -0.486). Consider using constrained = TRUE.
#> 
  |++++++++++++                                      | 23% ~01s          
#> Warning: Some alpha values < 0 (min = -0.742). Consider using constrained = TRUE.
#> 
  |+++++++++++++                                     | 24% ~01s          
  |+++++++++++++                                     | 25% ~01s          
#> Warning: Some alpha values < 0 (min = -0.561). Consider using constrained = TRUE.
#> 
  |++++++++++++++                                    | 26% ~01s          
  |++++++++++++++                                    | 27% ~01s          
  |+++++++++++++++                                   | 28% ~01s          
#> Warning: Some alpha values < 0 (min = -0.129). Consider using constrained = TRUE.
#> 
  |+++++++++++++++                                   | 29% ~00s          
  |++++++++++++++++                                  | 30% ~00s          
  |++++++++++++++++                                  | 31% ~00s          
#> Warning: Some alpha values < 0 (min = -0.258). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++                                 | 32% ~00s          
  |+++++++++++++++++                                 | 33% ~00s          
  |++++++++++++++++++                                | 34% ~00s          
#> Warning: Some alpha values < 0 (min = -0.097). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++                                | 35% ~00s          
#> Warning: Some alpha values < 0 (min = -0.253). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++                               | 36% ~00s          
  |+++++++++++++++++++                               | 37% ~00s          
  |++++++++++++++++++++                              | 38% ~00s          
#> Warning: Some alpha values < 0 (min = -0.281). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++                              | 39% ~00s          
#> Warning: Some alpha values < 0 (min = -0.222). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++                             | 40% ~00s          
#> Warning: Some alpha values < 0 (min = -0.335). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++                             | 41% ~00s          
  |++++++++++++++++++++++                            | 42% ~00s          
#> Warning: Some alpha values < 0 (min = -0.121). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++                            | 43% ~00s          
#> Warning: Some alpha values < 0 (min = -0.396). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++                           | 44% ~00s          
  |+++++++++++++++++++++++                           | 45% ~00s          
  |++++++++++++++++++++++++                          | 46% ~00s          
#> Warning: Some alpha values < 0 (min = -0.427). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++                          | 47% ~00s          
  |+++++++++++++++++++++++++                         | 48% ~00s          
#> Warning: Some alpha values < 0 (min = -0.47). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++                         | 49% ~00s          
#> Warning: Some alpha values < 0 (min = -0.092). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++                        | 51% ~00s          
  |++++++++++++++++++++++++++                        | 52% ~00s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++                       | 53% ~00s          
#> Warning: Some alpha values < 0 (min = -0.144). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++                       | 54% ~00s          
#> Warning: Some alpha values < 0 (min = -0.102). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++                      | 55% ~00s          
#> Warning: Some alpha values < 0 (min = -0.574). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++                      | 56% ~00s          
  |+++++++++++++++++++++++++++++                     | 57% ~00s          
#> Warning: Some alpha values < 0 (min = -0.328). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++                     | 58% ~00s          
#> Warning: Some alpha values < 0 (min = -0.396). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++                    | 59% ~00s          
  |++++++++++++++++++++++++++++++                    | 60% ~00s          
  |+++++++++++++++++++++++++++++++                   | 61% ~00s          
  |+++++++++++++++++++++++++++++++                   | 62% ~00s          
  |++++++++++++++++++++++++++++++++                  | 63% ~00s          
  |++++++++++++++++++++++++++++++++                  | 64% ~00s          
  |+++++++++++++++++++++++++++++++++                 | 65% ~00s          
  |+++++++++++++++++++++++++++++++++                 | 66% ~00s          
#> Warning: Some alpha values < 0 (min = -0.599). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++                | 67% ~00s          
  |++++++++++++++++++++++++++++++++++                | 68% ~00s          
#> Warning: Some alpha values < 0 (min = -0.408). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++               | 69% ~00s          
  |+++++++++++++++++++++++++++++++++++               | 70% ~00s          
#> Warning: Some alpha values < 0 (min = -0.087). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++              | 71% ~00s          
  |++++++++++++++++++++++++++++++++++++              | 72% ~00s          
  |+++++++++++++++++++++++++++++++++++++             | 73% ~00s          
  |+++++++++++++++++++++++++++++++++++++             | 74% ~00s          
#> Warning: Some alpha values < 0 (min = -0.199). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++            | 75% ~00s          
#> Warning: Some alpha values < 0 (min = -0.483). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++            | 76% ~00s          
#> Warning: Some alpha values < 0 (min = -0.716). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++++++           | 77% ~00s          
#> Warning: Some alpha values < 0 (min = -0.095). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++++++           | 78% ~00s          
  |++++++++++++++++++++++++++++++++++++++++          | 79% ~00s          
  |++++++++++++++++++++++++++++++++++++++++          | 80% ~00s          
#> Warning: Some alpha values < 0 (min = -0.393). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++++++++         | 81% ~00s          
#> Warning: Some alpha values < 0 (min = -0.058). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++++++++         | 82% ~00s          
#> Warning: Some alpha values < 0 (min = -0.018). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++++++        | 83% ~00s          
  |++++++++++++++++++++++++++++++++++++++++++        | 84% ~00s          
#> Warning: Some alpha values < 0 (min = -0.758). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++++++++++       | 85% ~00s          
#> Warning: Some alpha values < 0 (min = -0.601). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++++++++++       | 86% ~00s          
#> Warning: Some alpha values < 0 (min = -0.325). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++++++++      | 87% ~00s          
#> Warning: Some alpha values < 0 (min = -0.728). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++++++++      | 88% ~00s          
#> Warning: Some alpha values < 0 (min = -0.427). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++++++++++++     | 89% ~00s          
  |+++++++++++++++++++++++++++++++++++++++++++++     | 90% ~00s          
#> Warning: Some alpha values < 0 (min = -0.161). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++++++++++    | 91% ~00s          
#> Warning: Some alpha values < 0 (min = -0.187). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++++++++++    | 92% ~00s          
#> Warning: Some alpha values < 0 (min = -0.251). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++++++++++++++   | 93% ~00s          
#> Warning: Some alpha values < 0 (min = -0.339). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++++++++++++++   | 94% ~00s          
  |++++++++++++++++++++++++++++++++++++++++++++++++  | 95% ~00s          
#> Warning: Some alpha values < 0 (min = -0.05). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++++++++++++  | 96% ~00s          
#> Warning: Some alpha values < 0 (min = -0.73). Consider using constrained = TRUE.
#> 
  |+++++++++++++++++++++++++++++++++++++++++++++++++ | 97% ~00s          
  |+++++++++++++++++++++++++++++++++++++++++++++++++ | 98% ~00s          
#> Warning: Some alpha values < 0 (min = -0.011). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 99% ~00s          
#> Warning: Some alpha values < 0 (min = -0.114). Consider using constrained = TRUE.
#> 
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  
# }
```
