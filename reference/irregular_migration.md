# Irregular migration data for Poland (2019–2024)

Country-level panel data on irregular (unauthorized) migration in
Poland, combining administrative records from the Social Insurance
Institution (ZUS), Border Guard, and Police. Covers non-Schengen
countries of origin with at least one insured foreigner registered in
ZUS during 2019–2024.

## Usage

``` r
irregular_migration
```

## Format

A data frame with 1,382 rows and 8 variables:

- year:

  Year of observation (factor, 2019–2024).

- sex:

  Sex (factor: `"Female"`, `"Male"`).

- country_code:

  ISO 3166-1 alpha-3 country code (e.g. `"UKR"`).

- country:

  Country name in English.

- continent:

  Continent (`"Africa"`, `"Americas"`, `"Asia"`, `"Europe"`,
  `"Oceania"`).

- m:

  Observed count: foreigners apprehended by the Border Guard for
  unauthorized stay (the variable whose hidden population we estimate).

- n:

  Auxiliary count: foreigners identified by the Police (a second,
  partially overlapping administrative source).

- N:

  Reference population: foreigners registered in ZUS (social insurance).
  Used as a proxy for the total known foreign population from a given
  country.

## Source

- Social Insurance Institution (ZUS) — register of insured foreigners

- Border Guard — apprehensions for unauthorized stay

- Police — identification of foreigners

## Details

The dataset is used to estimate the size of the unauthorized foreign
population using the power-law model of Zhang (2008).

The key assumption of the estimation framework is that \\m\\ and \\n\\
are partial observations from a larger unauthorized population of size
\\\xi\\. The reference population \\N\\ (insured foreigners) serves as a
scaling anchor: countries with more insured foreigners are assumed to
also have more unauthorized residents, with the elasticity governed by
the parameter \\\alpha\\.

Observations with `N = 0` have been excluded. About 49\\ observations
have `m = 0` and 48\\ large number of small countries with no recorded
apprehensions.

## References

Zhang, L.-C. (2008). Developing methods for determining the number of
unauthorized foreigners in Norway. *Documents* 2008/11, Statistics
Norway.
<https://www.ssb.no/a/english/publikasjoner/pdf/doc_200811_en/doc_200811_en.pdf>

Beresewicz, M. and Pawlukiewicz, K. (2020). Estimation of the number of
irregular foreigners in Poland using non-linear count regression models.
*arXiv preprint* arXiv:2008.09407.

## Examples

``` r
data(irregular_migration)
head(irregular_migration)
#>   year    sex country_code     country continent m  n    N
#> 1 2019 Female          AFG Afghanistan      Asia 0  0  269
#> 2 2019 Female          AGO      Angola    Africa 2  0   15
#> 3 2019 Female          ALB     Albania    Europe 0  1   64
#> 4 2019 Female          ARG   Argentina  Americas 0  1   52
#> 5 2019 Female          ARM     Armenia      Asia 8 17 1089
#> 6 2019 Female          AUS   Australia   Oceania 2  0   44

# Basic model
fit <- estimate_hidden_pop(
  data = irregular_migration,
  observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
  method = "poisson"
)
summary(fit)
#> Unauthorized population estimation
#> Method: POISSON | vcov: HC3 
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
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>       Observed Estimate Estimate (BC) CI lower CI upper
#> (all)   27,105  290,597       290,498  127,726  660,703
```
