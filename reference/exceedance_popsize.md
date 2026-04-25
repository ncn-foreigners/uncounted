# Bootstrap Exceedance Probability for Population Size

Computes the empirical bootstrap tail area for a threshold question such
as “What fraction of bootstrap replications imply \\\xi\\ above a given
number?”

## Usage

``` r
exceedance_popsize(
  object,
  threshold,
  group = NULL,
  direction = c("above", "below")
)

# S3 method for class 'uncounted_popsize_exceedance'
print(x, ...)
```

## Arguments

- object:

  An `"uncounted_boot"` object returned by
  [`bootstrap_popsize`](https://ncn-foreigners.github.io/uncounted/reference/bootstrap_popsize.md).

- threshold:

  A single finite numeric threshold.

- group:

  Optional group label. When `NULL` (default), the helper uses the total
  bootstrap distribution if available; otherwise it uses the only group
  when the bootstrap object contains a single group.

- direction:

  Which tail area to compute: `"above"` for `P^*(xi > c)` or `"below"`
  for `P^*(xi < c)`.

- x:

  Object to print.

- ...:

  Additional arguments passed to print methods.

## Value

An object of class `"uncounted_popsize_exceedance"` with components
`group`, `threshold`, `direction`, `n_boot`, `n_finite`, `estimate`,
`count`, and `distribution_summary`.

## Details

The reported value is an empirical bootstrap exceedance probability,
i.e. a sample proportion over finite bootstrap draws. It is useful for
threshold questions, but it is not a Bayesian posterior probability.

For multi-group bootstrap objects with `total = TRUE`,
`exceedance_popsize()` reconstructs the total bootstrap distribution by
summing the per-replicate group estimates in `object$t`. This preserves
the within-replicate dependence across groups.

## Examples

``` r
# \donttest{
set.seed(42)
n_obs <- 15
sim_data <- data.frame(
  country = rep(paste0("C", 1:5), each = 3),
  year    = rep(2018:2020, 5),
  N       = round(exp(rnorm(n_obs, mean = 13, sd = 0.2)))
)
sim_data$n <- rpois(n_obs, lambda = pmax(1, 0.003 * sim_data$N))
sim_data$m <- rpois(n_obs, lambda = sim_data$N^0.6 *
  (0.005 + sim_data$n / sim_data$N)^0.8)

fit <- estimate_hidden_pop(
  data = sim_data, observed = ~m, auxiliary = ~n,
  reference_pop = ~N, method = "poisson",
  countries = ~country
)

if (requireNamespace("fwb", quietly = TRUE)) {
  boot <- bootstrap_popsize(
    fit, R = 49, cluster = ~country, seed = 123,
    total = TRUE, verbose = FALSE
  )
  exceedance_popsize(boot, threshold = 2000)
}
#> Bootstrap exceedance probability
#> Group: (all) 
#> P*(xi > 2000) = 0.9796 (48/49 finite draws)
#> Bootstrap draws: 49 
#> Distribution summary:
#>      mean    median        sd      q025      q975 
#> 23608.745 13468.667 24755.238  2441.417 84668.828 
# }
```
