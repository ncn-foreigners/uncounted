The {uncounted} package
================

## Basics

The {uncounted} package implements various methods to estimate
population size such as: graph-based methods, non-linear (count)
regression models or dual/multiple system estimation.

To install the package from GitHub you can use the `pak` package

``` r
install.packages("pak")
pak::pkg_install("ncn-foreigners/uncounted")
```

Load the package

``` r
library(uncounted)
```

Below you can find the basic functionality of the package on the Polish
example from 2018

``` r
data(foreigners_pl)
model_data <- subset(foreigners_pl, year == 2018 & half == 1)
head(model_data)
```

    ##     year half iso3n_new    sex border police pesel     country
    ## 153 2018    1       004   male      2     12    24 Afghanistan
    ## 154 2018    1       008 female      2      1   118     Albania
    ## 155 2018    1       008   male      2     10    35     Albania
    ## 156 2018    1       031 female      1      2   199  Azerbaijan
    ## 157 2018    1       031   male      9     16    80  Azerbaijan
    ## 158 2018    1       032   male      1      3    25   Argentina

Fit the basic Zhang model as described
[here](https://www.ssb.no/en/befolkning/artikler-og-publikasjoner/developing-methods-for-determining-the-number-of-unauthorized-foreigners-in-norway--18906)
to obtain estimate of the irregular population size along with 95%
confidence interval.

``` r
result_zhang <- estimate_hidden_pop(
   data = model_data,
   observed = ~ border,
   auxiliary = ~ police,
   reference_pop = ~ pesel,
   method = "mle",
   family = "poisson"
)

c(N_hat=result_zhang$estimates$xi, 
  lower=unname(result_zhang$confint_xi[1]),
  upper=unname(result_zhang$confint_xi[2]))
```

    ##     N_hat     lower     upper 
    ## 10420.367  5485.292 20261.303

## Acknowledgements

- Work on this package is supported by the National Science Centre, OPUS
  20 grant no. 2020/39/B/HS4/00941 (Towards census-like statistics for
  foreign-born populations – quality, data integration and estimation)

- Name of the package was motivated by the name of the [Zult, D. B.
  (2024). Counting the uncounted: Methodological extensions in multiple
  systems estimation. Dissertation, Maastricht University,
  doi:10.33540/2611](https://www.cbs.nl/en-gb/background/2024/51/counting-the-uncounted)
