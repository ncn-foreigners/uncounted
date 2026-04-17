# iOLS: Iterated OLS for Gamma Pseudo-Maximum Likelihood

## Introduction

The `uncounted` package supports five estimation methods for the
power-law model
$E\left( m_{i} \right) = N_{i}^{\alpha_{i}}\left( \gamma + n_{i}/N_{i} \right)^{\beta_{i}}$.
This vignette describes the iOLS (iterated Ordinary Least Squares)
estimator, which targets the Gamma Pseudo-Maximum Likelihood (GPML)
solution via a sequence of linear regressions.

The method was proposed by Benatia, Bellego and Pape (2024) as a
computationally simple alternative to nonlinear GLM estimation. The key
advantage: each iteration is a standard OLS regression, so no gradient
or Hessian computation is needed.

## Model recap

On the log scale the model is linear:

$$\log E\left( m_{i} \right) = \alpha_{i}\log N_{i} + \beta_{i}\log\left( \gamma + n_{i}/N_{i} \right) = \mathbf{z}_{i}\prime{\mathbf{θ}}$$

where $\mathbf{z}_{i}$ is the design vector and
${\mathbf{θ}} = (\alpha,\beta)\prime$. This is an exponential
conditional mean model
$E\left( m_{i}|\mathbf{z}_{i} \right) = \exp\left( \mathbf{z}_{i}\prime{\mathbf{θ}} \right)$.

## GPML vs PPML

There are two leading PML estimators for this model:

**Poisson PML (PPML)** solves
$\sum_{i}\left( m_{i} - \mu_{i} \right)\mathbf{z}_{i} = \mathbf{0}$. It
weights observations by $\mu_{i}$ in the information matrix, giving more
influence to large-count observations.

**Gamma PML (GPML)** solves
$\sum_{i}\left( m_{i}/\mu_{i} - 1 \right)\mathbf{z}_{i} = \mathbf{0}$.
It gives equal weight to all observations regardless of $\mu_{i}$,
because the Gamma variance function $V(\mu) = \mu^{2}$ implies unit
working weights under a log link.

Both are consistent for $\mathbf{θ}$ under correct specification of the
conditional mean. They differ in efficiency and sensitivity to outliers:
PPML is more efficient when the data are Poisson-like, while GPML is
more robust to observations with unusually large counts relative to
$\mu$.

## The iOLS algorithm

iOLS solves the GPML score equations via two phases.

### Phase 1: Warm-up with increasing delta

Initialize ${\widehat{\mathbf{θ}}}_{0}$ from OLS on $\log(m + 1)$. For
each $\delta$ in an increasing sequence (default: 1, 10, 100, 1000):

1.  Compute
    ${\widehat{\mu}}_{i} = \exp\left( \mathbf{z}_{i}\prime{\widehat{\mathbf{θ}}}_{t} \right)$
2.  Set
    ${\widetilde{y}}_{i} = \log\left( m_{i} + \delta{\widehat{\mu}}_{i} \right) - \bar{c}$,
    where $\bar{c}$ is the weighted mean of
    $\log\left( m_{i} + \delta{\widehat{\mu}}_{i} \right)$ (empirical
    centering)
3.  Update ${\widehat{\mathbf{θ}}}_{t + 1}$ by regressing
    $\widetilde{y}$ on $\mathbf{Z}$ using OLS
4.  Repeat until convergence

This phase is guaranteed to converge globally (contraction mapping with
modulus $\kappa < 1$). Larger $\delta$ gives a better approximation to
GPML but slower convergence.

### Phase 2: Exact GPML limiting transform

After Phase 1 converges, switch to the exact debiasing step:

$${\widetilde{y}}_{i} = \log{\widehat{\mu}}_{i} + \frac{m_{i}/{\widehat{\mu}}_{i} - 1}{1 + \rho}$$

where $\rho > 0$ is a stability parameter (default: 1). This is the
theoretically correct form from Benatia et al. (2024) that targets the
exact GPML score. At convergence:

$$\sum\limits_{i}\left( m_{i}/{\widehat{\mu}}_{i} - 1 \right)\mathbf{z}_{i} = \mathbf{0}$$

which are the GPML first-order conditions.

## Usage

``` r
library(uncounted)
data(irregular_migration)
d <- irregular_migration
d$ukr <- as.integer(d$country_code == "UKR")

# iOLS requires a fixed gamma (estimated gamma not yet supported)
fit_iols <- estimate_hidden_pop(
  data = d, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
  method = "iols",
  cov_alpha = ~ year * ukr,
  cov_beta = ~ year,
  gamma = 0.005,
  countries = ~ country_code
)

summary(fit_iols)
#> Unauthorized population estimation
#> Method: IOLS | estimator: MLE | link_rho: power | vcov: HC3 
#> N obs: 1382 
#> Gamma: 0.005 (fixed) 
#> Log-likelihood: -462765.9 
#> AIC: 925567.8  BIC: 925662 
#> Deviance: 2605.75 
#> 
#> Coefficients:
#>                     Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)   0.772397   0.047736 16.1807   <2e-16 ***
#> alpha:year2020     -0.158556   0.120951 -1.3109   0.1899    
#> alpha:year2021     -0.062778   0.113575 -0.5527   0.5804    
#> alpha:year2022     -0.084999   0.117596 -0.7228   0.4698    
#> alpha:year2023     -0.117000   0.139903 -0.8363   0.4030    
#> alpha:year2024      0.071404   0.091504  0.7803   0.4352    
#> alpha:ukr           0.041534   0.028920  1.4362   0.1510    
#> alpha:year2020:ukr  0.080418   0.068406  1.1756   0.2398    
#> alpha:year2021:ukr  0.022852   0.064421  0.3547   0.7228    
#> alpha:year2022:ukr -0.064579   0.083448 -0.7739   0.4390    
#> alpha:year2023:ukr -0.094507   0.102574 -0.9214   0.3569    
#> alpha:year2024:ukr -0.223856   0.085450 -2.6197   0.0088 ** 
#> beta:(Intercept)    0.601292   0.067822  8.8657   <2e-16 ***
#> beta:year2020      -0.052293   0.160779 -0.3252   0.7450    
#> beta:year2021       0.164709   0.155589  1.0586   0.2898    
#> beta:year2022       0.064846   0.172962  0.3749   0.7077    
#> beta:year2023      -0.061620   0.186935 -0.3296   0.7417    
#> beta:year2024       0.101640   0.131161  0.7749   0.4384    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> -----------------------
#> Population size estimation results:
#>   (BC = bias-corrected using model-based variance)
#>                  Observed Estimate Estimate (BC) CI lower CI upper
#> year=2019, ukr=0    2,603   21,829        21,321   10,383   43,780
#> year=2019, ukr=1    4,001   47,258        36,424   19,260   68,886
#> year=2020, ukr=0    1,351    7,329         7,183    1,494   34,548
#> year=2020, ukr=1    2,047   19,407        14,973    4,305   52,079
#> year=2021, ukr=0    1,481   18,128        17,734    3,612   87,075
#> year=2021, ukr=1    1,624   35,513        27,396    8,466   88,653
#> year=2022, ukr=0    2,114   18,056        17,726    3,196   98,296
#> year=2022, ukr=1      835   10,840         8,398    1,258   56,058
#> year=2023, ukr=0    3,727   15,730        15,455    1,904  125,473
#> year=2023, ukr=1      635    4,834         3,745      446   31,421
#> year=2024, ukr=0    5,928   86,505        84,645   21,594  331,802
#> year=2024, ukr=1      759   10,000         7,744    1,115   53,797
popsize(fit_iols)
#>               group observed  estimate estimate_bc      lower     upper
#> 1  year=2019, ukr=0     2603 21828.966   21320.974 10383.3310  43780.16
#> 2  year=2019, ukr=1     4001 47257.597   36424.421 19259.8263  68886.32
#> 3  year=2020, ukr=0     1351  7329.221    7183.299  1493.5486  34548.45
#> 4  year=2020, ukr=1     2047 19407.135   14973.485  4305.1368  52078.54
#> 5  year=2021, ukr=0     1481 18127.555   17733.917  3611.7351  87074.99
#> 6  year=2021, ukr=1     1624 35512.915   27395.617  8465.8230  88652.91
#> 7  year=2022, ukr=0     2114 18055.980   17725.590  3196.4209  98296.36
#> 8  year=2022, ukr=1      835 10840.200    8397.878  1258.0551  56058.25
#> 9  year=2023, ukr=0     3727 15729.967   15454.610  1903.5493 125473.49
#> 10 year=2023, ukr=1      635  4833.874    3744.789   446.3124  31420.70
#> 11 year=2024, ukr=0     5928 86504.752   84645.305 21593.7109 331801.59
#> 12 year=2024, ukr=1      759 10000.132    7744.133  1114.7721  53797.18
#>    share_pct
#> 1   7.388922
#> 2  15.996300
#> 3   2.480880
#> 4   6.569152
#> 5   6.136025
#> 6  12.020824
#> 7   6.111798
#> 8   3.669317
#> 9   5.324462
#> 10  1.636226
#> 11 29.281133
#> 12  3.384961
```

## Comparison with Poisson

On simulated data with known parameters ($\alpha = 0.70$,
$\beta = 0.55$, $\gamma = 0.005$), all methods recover the truth:

| Method  | $\widehat{\alpha}$ | $\widehat{\beta}$ | $\widehat{\xi}$ |
|---------|--------------------|-------------------|-----------------|
| OLS     | 0.713              | 0.583             | 154,609         |
| Poisson | 0.722              | 0.602             | 168,778         |
| NB      | 0.722              | 0.603             | 169,244         |
| iOLS    | 0.714              | 0.587             | 156,714         |

iOLS estimates are closer to OLS than to Poisson/NB, reflecting the
different weighting: GPML gives equal weight to all observations, while
PPML upweights large-$\mu$ observations.

On real migration data with ~50% zeros and complex covariates, iOLS
typically gives **lower** total population estimates than Poisson (about
50–70% of the Poisson total), because GPML downweights the few large
countries that dominate the PPML fit.

## Variance estimation

iOLS uses a sandwich variance estimator with the GPML score residual
$e_{i} = m_{i}/{\widehat{\mu}}_{i} - 1$:

$$\widehat{\mathbf{V}} = (\mathbf{Z}\prime\mathbf{Z})^{- 1}\mathbf{Z}\prime\text{diag}\left( {\widehat{e}}_{i}^{2} \right)\mathbf{Z}(\mathbf{Z}\prime\mathbf{Z})^{- 1}$$

The bread $(\mathbf{Z}\prime\mathbf{Z})^{- 1}$ reflects the GPML Fisher
information (unit working weights from the Gamma variance function). HC0
through HC5 adjustments are available via the `vcov` argument.

Standard errors are generally smaller than Poisson because GPML is more
efficient when the conditional variance is proportional to $\mu^{2}$
(Gamma-like). On count data, Poisson SEs may be more conservative.

## Bias correction for population size

The plug-in estimator
${\widehat{\xi}}_{g} = \sum_{i \in g}N_{i}^{{\widehat{\alpha}}_{g}}$ is
upward-biased by Jensen’s inequality. The bias correction uses the
model-based variance $\mathbf{V}_{\text{model}}$, which differs by
method:

| Method                             | $\mathbf{V}_{\text{model}}$                                                                                                                              |
|------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------|
| Poisson / Poisson GMM / Poisson EL | $\left( \mathbf{Z}\prime\text{diag}\left( \widehat{\mu} \right)\mathbf{Z} \right)^{- 1}$                                                                 |
| NB / NB GMM / NB EL                | $\left( \mathbf{Z}\prime\text{diag}\left( \widehat{\mu}\widehat{\theta}/\left( \widehat{\theta} + \widehat{\mu} \right) \right)\mathbf{Z} \right)^{- 1}$ |
| iOLS                               | $(\mathbf{Z}\prime\mathbf{Z})^{- 1}$                                                                                                                     |

For iOLS, the model-based variance is the GPML Fisher information
inverse **without** dispersion scaling. This is because the Gamma
quasi-likelihood has variance function $V(\mu) = \mu^{2}$, giving unit
working weights.

For unconstrained models, including iOLS, the current implementation
uses the multiplicative lognormal correction

$${\widehat{\xi}}^{BC} = \sum\limits_{i}N_{i}^{{\widehat{\alpha}}_{i}}\exp\!\left( - \frac{1}{2}\left( \log N_{i} \right)^{2}\mathbf{x}_{i}^{\top}{\widehat{\mathbf{V}}}_{\text{model}}\mathbf{x}_{i} \right),$$

which is exact under the Gaussian approximation for the linear predictor
and guarantees positivity. Constrained count fits fall back to a
subtractive second-order Taylor correction because the corresponding
logistic-normal integral has no closed form.

## Convergence diagnostics

Check `fit$convergence`:

- `0`: converged (GPML score condition satisfied)
- `1`: did not converge (consider simpler covariates or checking data)

The normalized GPML score `max|Z'(m/mu - 1)| / n` should be small (\<
$10^{- 4}$) at convergence.

## Limitations

1.  **`gamma = "estimate"` not supported.** The paper’s theory covers
    iOLS on the mean model; profiled gamma is a package extension not
    yet validated. Use a fixed gamma from Poisson estimation.

2.  **HC2/HC3 not available for the NB theta-aware path.** This only
    affects NB, not iOLS.

3.  **AIC/BIC not comparable across methods.** iOLS uses a GPML pseudo
    log-likelihood, which is not on the same scale as Poisson/NB count
    likelihoods.
    [`compare_models()`](https://ncn-foreigners.github.io/uncounted/reference/compare_models.md)
    warns about this.

4.  **Bootstrap recommended for inference.** The analytical bias
    correction uses the model-based variance, which may underestimate
    the true Jensen’s bias on zero-heavy data. Bootstrap bias correction
    `2*xi_hat - mean(xi_star)` is more robust.

## References

- Benatia, D., Bellego, C. and Pape, L.-D. (2024). Dealing with Logs and
  Zeros in Regression Models. *arXiv preprint* arXiv:2203.11820v3.
- Santos Silva, J. M. C. and Tenreyro, S. (2006). The log of gravity.
  *The Review of Economics and Statistics*, 88(4), 641–658.
- Gourieroux, C., Monfort, A. and Trognon, A. (1984). Pseudo Maximum
  Likelihood Methods: Applications to Poisson Models. *Econometrica*,
  52(3), 701–720.
- Zhang, L.-C. (2008). *Developing methods for determining the number of
  unauthorized foreigners in Norway* (Documents 2008/11). Statistics
  Norway.
  <https://www.ssb.no/a/english/publikasjoner/pdf/doc_200811_en/doc_200811_en.pdf>
- Beresewicz, M. and Pawlukiewicz, K. (2020). Estimation of the number
  of irregular foreigners in Poland using non-linear count regression
  models. *arXiv preprint* arXiv:2008.09407.
