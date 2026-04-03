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
  ci_type = c("perc", "bc"),
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

  Type of bootstrap CI: `"perc"` (percentile) or `"bc"` (bias-corrected
  percentile). Default `"perc"`.

- point_estimate:

  Which point estimate to report: `"median"` (default, recommended) –
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

## Value

An object of class `"uncounted_boot"` with components:

- t:

  Matrix (`R` x `n_groups`): bootstrap population size estimates per
  group. Each row is one replicate.

- t0:

  Numeric vector: plug-in point estimates from original fit.

- t0_bc:

  Numeric vector: bias-corrected point estimates (Taylor expansion).

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

  CI type used (`"perc"` or `"bc"`).

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
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.

# \donttest{
# Cluster bootstrap (recommended when data has repeated country obs)
boot <- bootstrap_popsize(fit, R = 99, cluster = ~country, seed = 123)
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |                                                  | 0 % ~calculating  
#> Warning: Some alpha values < 0 (min = -0.081). Consider using constrained = TRUE.
#>   |+                                                 | 1 % ~00s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |++                                                | 2 % ~00s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#>   |++                                                | 3 % ~00s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |+++                                               | 4 % ~02s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++                                               | 5 % ~03s          
#> Warning: Some alpha values < 0 (min = -0.077). Consider using constrained = TRUE.
#>   |++++                                              | 6 % ~03s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++                                              | 7 % ~02s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |+++++                                             | 8 % ~02s          
#> Warning: Some alpha values < 0 (min = -0.471). Consider using constrained = TRUE.
#>   |+++++                                             | 9 % ~02s          
#> Warning: Some alpha values < 0 (min = -2.087). Consider using constrained = TRUE.
#>   |++++++                                            | 10% ~02s          
#> Warning: Some alpha values < 0 (min = -0.076). Consider using constrained = TRUE.
#>   |++++++                                            | 11% ~02s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++                                           | 12% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++                                           | 13% ~01s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |++++++++                                          | 14% ~01s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#>   |++++++++                                          | 15% ~01s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |+++++++++                                         | 16% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++++                                         | 17% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++                                        | 18% ~01s          
#> Warning: Some alpha values < 0 (min = -0.588). Consider using constrained = TRUE.
#>   |++++++++++                                        | 19% ~01s          
#> Warning: Some alpha values < 0 (min = -0.465). Consider using constrained = TRUE.
#>   |+++++++++++                                       | 20% ~01s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |+++++++++++                                       | 21% ~01s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |++++++++++++                                      | 22% ~01s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |++++++++++++                                      | 23% ~01s          
#> Warning: Some alpha values < 0 (min = -0.078). Consider using constrained = TRUE.
#>   |+++++++++++++                                     | 24% ~01s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |+++++++++++++                                     | 25% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++                                    | 26% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++                                    | 27% ~01s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |+++++++++++++++                                   | 28% ~01s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |+++++++++++++++                                   | 29% ~01s          
#> Warning: Some alpha values < 0 (min = -0.465). Consider using constrained = TRUE.
#>   |++++++++++++++++                                  | 30% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++++                                  | 31% ~01s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |+++++++++++++++++                                 | 32% ~01s          
#> Warning: Some alpha values < 0 (min = -0.081). Consider using constrained = TRUE.
#>   |+++++++++++++++++                                 | 33% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++++++                                | 34% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++++++                                | 35% ~01s          
#> Warning: Some alpha values < 0 (min = -0.47). Consider using constrained = TRUE.
#>   |+++++++++++++++++++                               | 36% ~01s          
#> Warning: Some alpha values < 0 (min = -0.471). Consider using constrained = TRUE.
#>   |+++++++++++++++++++                               | 37% ~01s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#>   |++++++++++++++++++++                              | 38% ~01s          
#> Warning: Some alpha values < 0 (min = -0.6). Consider using constrained = TRUE.
#>   |++++++++++++++++++++                              | 39% ~01s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++                             | 40% ~01s          
#> Warning: Some alpha values < 0 (min = -0.077). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++                             | 41% ~01s          
#> Warning: Some alpha values < 0 (min = -0.596). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++                            | 42% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++                            | 43% ~01s          
#> Warning: Some alpha values < 0 (min = -0.076). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++                           | 44% ~01s          
#> Warning: Some alpha values < 0 (min = -2.562). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++                           | 45% ~01s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++                          | 46% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++                          | 47% ~01s          
#> Warning: Some alpha values < 0 (min = -0.465). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++                         | 48% ~01s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++                         | 49% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++                        | 51% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++                        | 52% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++                       | 53% ~01s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++                       | 54% ~00s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++                      | 55% ~00s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++                      | 56% ~00s          
#> Warning: Some alpha values < 0 (min = -2.356). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++                     | 57% ~00s          
#> Warning: Some alpha values < 0 (min = -0.291). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++                     | 58% ~00s          
#> Warning: Some alpha values < 0 (min = -0.47). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++                    | 59% ~00s          
#> Warning: Some alpha values < 0 (min = -0.078). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++                    | 60% ~00s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++                   | 61% ~00s          
#> Warning: Some alpha values < 0 (min = -0.465). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++                   | 62% ~00s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++                  | 63% ~00s          
#> Warning: Some alpha values < 0 (min = -0.464). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++                  | 64% ~00s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++                 | 65% ~00s          
#> Warning: Some alpha values < 0 (min = -0.08). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++                 | 66% ~00s          
#> Warning: Some alpha values < 0 (min = -2.694). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++                | 67% ~00s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++                | 68% ~00s          
#> Warning: Some alpha values < 0 (min = -2.539). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++               | 69% ~00s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++               | 70% ~00s          
#> Warning: Some alpha values < 0 (min = -0.462). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++              | 71% ~00s          
#> Warning: Some alpha values < 0 (min = -0.465). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++              | 72% ~00s          
#> Warning: Some alpha values < 0 (min = -0.465). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++             | 73% ~00s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++             | 74% ~00s          
#> Warning: Some alpha values < 0 (min = -0.47). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++            | 75% ~00s          
#> Warning: Some alpha values < 0 (min = -2.008). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++            | 76% ~00s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++           | 77% ~00s          
#> Warning: Some alpha values < 0 (min = -0.47). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++           | 78% ~00s          
#> Warning: Some alpha values < 0 (min = -0.498). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++          | 79% ~00s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++          | 80% ~00s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++         | 81% ~00s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++         | 82% ~00s          
#> Warning: Some alpha values < 0 (min = -0.465). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++        | 83% ~00s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++        | 84% ~00s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++++       | 85% ~00s          
#> Warning: Some alpha values < 0 (min = -1.934). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++++       | 86% ~00s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++++      | 87% ~00s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++++      | 88% ~00s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++++++     | 89% ~00s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++++++     | 90% ~00s          
#> Warning: Some alpha values < 0 (min = -2.205). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++++++    | 91% ~00s          
#> Warning: Some alpha values < 0 (min = -0.077). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++++++    | 92% ~00s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++++++++   | 93% ~00s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++++++++   | 94% ~00s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++++++++  | 95% ~00s          
#> Warning: Some alpha values < 0 (min = -0.47). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++++++++  | 96% ~00s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++++++++++ | 97% ~00s          
#> Warning: Some alpha values < 0 (min = -2.033). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++++++++++ | 98% ~00s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++++++++++| 99% ~00s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  
boot                # prints table with median as point estimate
#> Bootstrap population size estimation
#> R = 99 | CI type: perc | Point estimate: median | Converged: 99 / 99 
#> Cluster bootstrap
#> 95% CI
#> 
#>   Point estimate: bootstrap median (recommended) | CI: bootstrap percentile
#>       Plugin Plugin (BC) Boot median Boot mean CI lower CI upper
#> (all)      0      -1,978           0         1        0        5
boot$popsize        # data.frame for further use
#>   group   estimate        lower    upper
#> 1 (all) 0.03266777 2.817751e-13 5.450315
boot$popsize_full   # all estimate types side by side
#>   group     plugin plugin_bc boot_median boot_mean        lower    upper
#> 1 (all) 0.03262174  -1977.66  0.03266777 0.5721444 2.817751e-13 5.450315

# Bias-corrected percentile CI
boot_bc <- bootstrap_popsize(fit, R = 99, ci_type = "bc", seed = 123)
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |                                                  | 0 % ~calculating  
#> Warning: Some alpha values < 0 (min = -0.111). Consider using constrained = TRUE.
#>   |+                                                 | 1 % ~09s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++                                                | 2 % ~04s          
#> Warning: Some alpha values < 0 (min = -0.47). Consider using constrained = TRUE.
#>   |++                                                | 3 % ~03s          
#> Warning: Some alpha values < 0 (min = -0.192). Consider using constrained = TRUE.
#>   |+++                                               | 4 % ~04s          
#> Warning: Some alpha values < 0 (min = -0.433). Consider using constrained = TRUE.
#>   |+++                                               | 5 % ~04s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |++++                                              | 6 % ~04s          
#> Warning: Some alpha values < 0 (min = -0.078). Consider using constrained = TRUE.
#>   |++++                                              | 7 % ~03s          
#> Warning: Some alpha values < 0 (min = -0.463). Consider using constrained = TRUE.
#>   |+++++                                             | 8 % ~03s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |+++++                                             | 9 % ~03s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++                                            | 10% ~02s          
#> Warning: Some alpha values < 0 (min = -1.007). Consider using constrained = TRUE.
#>   |++++++                                            | 11% ~02s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++                                           | 12% ~02s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |+++++++                                           | 13% ~02s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++                                          | 14% ~02s          
#> Warning: Some alpha values < 0 (min = -2.028). Consider using constrained = TRUE.
#>   |++++++++                                          | 15% ~02s          
#> Warning: Some alpha values < 0 (min = -0.086). Consider using constrained = TRUE.
#>   |+++++++++                                         | 16% ~02s          
#> Warning: Some alpha values < 0 (min = -2.355). Consider using constrained = TRUE.
#>   |+++++++++                                         | 17% ~02s          
#> Warning: Some alpha values < 0 (min = -2.141). Consider using constrained = TRUE.
#>   |++++++++++                                        | 18% ~02s          
#> Warning: Some alpha values < 0 (min = -1.789). Consider using constrained = TRUE.
#>   |++++++++++                                        | 19% ~02s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#>   |+++++++++++                                       | 20% ~02s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++++++                                       | 21% ~02s          
#> Warning: Some alpha values < 0 (min = -0.081). Consider using constrained = TRUE.
#>   |++++++++++++                                      | 22% ~01s          
#> Warning: Some alpha values < 0 (min = -0.412). Consider using constrained = TRUE.
#>   |++++++++++++                                      | 23% ~02s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++++++++                                     | 24% ~01s          
#> Warning: Some alpha values < 0 (min = -0.465). Consider using constrained = TRUE.
#>   |+++++++++++++                                     | 25% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++                                    | 26% ~01s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |++++++++++++++                                    | 27% ~01s          
#> Warning: Some alpha values < 0 (min = -1.08). Consider using constrained = TRUE.
#>   |+++++++++++++++                                   | 28% ~01s          
#> Warning: Some alpha values < 0 (min = -0.463). Consider using constrained = TRUE.
#>   |+++++++++++++++                                   | 29% ~01s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |++++++++++++++++                                  | 30% ~01s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |++++++++++++++++                                  | 31% ~01s          
#> Warning: Some alpha values < 0 (min = -1.771). Consider using constrained = TRUE.
#>   |+++++++++++++++++                                 | 32% ~01s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |+++++++++++++++++                                 | 33% ~01s          
#> Warning: Some alpha values < 0 (min = -0.47). Consider using constrained = TRUE.
#>   |++++++++++++++++++                                | 34% ~01s          
#> Warning: Some alpha values < 0 (min = -0.376). Consider using constrained = TRUE.
#>   |++++++++++++++++++                                | 35% ~01s          
#> Warning: Some alpha values < 0 (min = -0.464). Consider using constrained = TRUE.
#>   |+++++++++++++++++++                               | 36% ~01s          
#> Warning: Some alpha values < 0 (min = -0.079). Consider using constrained = TRUE.
#>   |+++++++++++++++++++                               | 37% ~01s          
#> Warning: Some alpha values < 0 (min = -2.344). Consider using constrained = TRUE.
#>   |++++++++++++++++++++                              | 38% ~01s          
#> Warning: Some alpha values < 0 (min = -0.47). Consider using constrained = TRUE.
#>   |++++++++++++++++++++                              | 39% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++                             | 40% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++                             | 41% ~01s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++                            | 42% ~01s          
#> Warning: Some alpha values < 0 (min = -0.14). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++                            | 43% ~01s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++                           | 44% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++                           | 45% ~01s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++                          | 46% ~01s          
#> Warning: Some alpha values < 0 (min = -0.078). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++                          | 47% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++                         | 48% ~01s          
#> Warning: Some alpha values < 0 (min = -0.47). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++                         | 49% ~01s          
#> Warning: Some alpha values < 0 (min = -0.47). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++                        | 51% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++                        | 52% ~01s          
#> Warning: Some alpha values < 0 (min = -0.862). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++                       | 53% ~01s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++                       | 54% ~01s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++                      | 55% ~01s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++                      | 56% ~01s          
#> Warning: Some alpha values < 0 (min = -0.47). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++                     | 57% ~01s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++                     | 58% ~01s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++                    | 59% ~01s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++                    | 60% ~01s          
#> Warning: Some alpha values < 0 (min = -0.47). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++                   | 61% ~01s          
#> Warning: Some alpha values < 0 (min = -2.044). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++                   | 62% ~01s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++                  | 63% ~00s          
#> Warning: Some alpha values < 0 (min = -0.233). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++                  | 64% ~01s          
#> Warning: Some alpha values < 0 (min = -1.79). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++                 | 65% ~00s          
#> Warning: Some alpha values < 0 (min = -0.463). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++                 | 66% ~00s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++                | 67% ~00s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++                | 68% ~00s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++               | 69% ~00s          
#> Warning: Some alpha values < 0 (min = -0.018). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++               | 70% ~00s          
#> Warning: Some alpha values < 0 (min = -0.462). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++              | 71% ~00s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++              | 72% ~00s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++             | 73% ~00s          
#> Warning: Some alpha values < 0 (min = -0.124). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++             | 74% ~00s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++            | 75% ~00s          
#> Warning: Some alpha values < 0 (min = -0.465). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++            | 76% ~00s          
#> Warning: Some alpha values < 0 (min = -0.465). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++           | 77% ~00s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++           | 78% ~00s          
#> Warning: Some alpha values < 0 (min = -0.47). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++          | 79% ~00s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++          | 80% ~00s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++         | 81% ~00s          
#> Warning: Some alpha values < 0 (min = -0.463). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++         | 82% ~00s          
#> Warning: Some alpha values < 0 (min = -2.795). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++        | 83% ~00s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++        | 84% ~00s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++++       | 85% ~00s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++++       | 86% ~00s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++++      | 87% ~00s          
#> Warning: Some alpha values < 0 (min = -0.23). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++++      | 88% ~00s          
#> Warning: Some alpha values < 0 (min = -0.266). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++++++     | 89% ~00s          
#> Warning: Some alpha values < 0 (min = -0.468). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++++++     | 90% ~00s          
#> Warning: Some alpha values < 0 (min = -0.078). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++++++    | 91% ~00s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++++++    | 92% ~00s          
#> Warning: Some alpha values < 0 (min = -0.469). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++++++++   | 93% ~00s          
#> Warning: Some alpha values < 0 (min = -0.467). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++++++++   | 94% ~00s          
#> Warning: Some alpha values < 0 (min = -0.466). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++++++++  | 95% ~00s          
#> Warning: Some alpha values < 0 (min = -0.465). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++++++++  | 96% ~00s          
#> Warning: Some alpha values < 0 (min = -2.177). Consider using constrained = TRUE.
#>   |+++++++++++++++++++++++++++++++++++++++++++++++++ | 97% ~00s            |+++++++++++++++++++++++++++++++++++++++++++++++++ | 98% ~00s          
#> Warning: Some alpha values < 0 (min = -0.47). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++++++++++| 99% ~00s          
#> Warning: Some alpha values < 0 (min = -0.31). Consider using constrained = TRUE.
#>   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  
# }
```
