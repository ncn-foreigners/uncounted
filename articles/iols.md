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
popsize(fit_iols)
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

| Method  | $\mathbf{V}_{\text{model}}$                                                              |
|---------|------------------------------------------------------------------------------------------|
| Poisson | $\left( \mathbf{Z}\prime\text{diag}\left( \widehat{\mu} \right)\mathbf{Z} \right)^{- 1}$ |
| iOLS    | $(\mathbf{Z}\prime\mathbf{Z})^{- 1}$                                                     |

For iOLS, the model-based variance is the GPML Fisher information
inverse **without** dispersion scaling. This is because the Gamma
quasi-likelihood has variance function $V(\mu) = \mu^{2}$, giving unit
working weights.

The current implementation uses a subtractive correction
${\widehat{\xi}}^{BC} = \widehat{\xi} - \text{bias}$, clamped to be
positive. A multiplicative correction based on the exact Gaussian MGF
${\widehat{\xi}}^{BC} = \sum N_{i}^{\widehat{\alpha}}\exp\left( - \frac{1}{2}\left( \log N_{i} \right)^{2}\mathbf{x}_{i}\prime\widehat{\mathbf{V}}\mathbf{x}_{i} \right)$
is under development and would eliminate the negativity issue entirely.

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
