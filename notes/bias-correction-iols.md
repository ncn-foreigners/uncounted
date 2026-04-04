# Bias Correction for Population Size Estimates under GPML/iOLS

## Problem Statement

The `uncounted` package estimates unauthorized population size as
$$\hat{\xi}_g = \sum_{i \in g} N_i^{\hat{\alpha}_g}$$
where $\hat{\alpha}$ is estimated by Poisson PML, Negative Binomial MLE, OLS, or iOLS (GPML).

Because $\xi(\alpha) = \sum N_i^\alpha$ is a convex function of $\alpha$ (for $N_i > 1$), Jensen's inequality implies $E[\hat{\xi}] \geq \xi(\alpha_0)$, so the plug-in estimate is upward-biased. The package applies a second-order Taylor bias correction:

$$\text{Bias}(\hat{\xi}_g) \approx \frac{1}{2} \sum_{i \in g} N_i^{\alpha_g} (\log N_i)^2 \cdot \mathbf{x}_i' \mathbf{V} \mathbf{x}_i$$

where $\mathbf{V}$ is a variance-covariance matrix of $\hat{\alpha}$.

**The question: which $\mathbf{V}$ should be used?**

## Current Implementation

The code uses $\mathbf{V}_{\text{model}}$: the **model-based (homoscedastic) variance**, not the sandwich (HC-robust) variance. This is stored as `vcov_model` on the fitted object.

For each estimation method, `vcov_model` is computed as:

| Method | $\mathbf{V}_{\text{model}}$ formula | What is $\sigma^2$? |
|--------|--------------------------------------|---------------------|
| OLS | $\sigma^2 (\mathbf{Z}'\mathbf{Z})^{-1}$ | $\sum (\log m_i - \mathbf{z}_i'\hat{\beta})^2 / (n-p)$ |
| Poisson | $(\mathbf{Z}' \text{diag}(\hat{\mu}) \mathbf{Z})^{-1}$ | Not needed (Fisher information) |
| NB | $(\mathbf{Z}' \text{diag}(w_i^{NB}) \mathbf{Z})^{-1}$ | Not needed |
| **iOLS** | $\sigma^2 (\mathbf{Z}'\mathbf{Z})^{-1}$ | $\sum (m_i/\hat{\mu}_i - 1)^2 / (n-p)$ |

## Why iOLS Has the Problem

For iOLS (GPML), the "working residuals" are $e_i = m_i/\hat{\mu}_i - 1$ (the GPML score contributions). These residuals have very different behavior from Poisson/NB:

- When $m_i = 0$: $e_i = -1$ (bounded)
- When $m_i \gg \hat{\mu}_i$ (large underprediction): $e_i = m_i/\hat{\mu}_i - 1 \gg 1$ (unbounded above)

With ~50% zeros in migration data, $\sigma^2 = \sum e_i^2 / (n-p)$ is dominated by a few observations where $m_i/\hat{\mu}_i$ is large. This inflates $\mathbf{V}_{\text{model}} = \sigma^2 (\mathbf{Z}'\mathbf{Z})^{-1}$, which inflates the bias term beyond the point estimate for groups with high-variance interaction terms (like UKR with `year*ukr` covariates).

**Concrete example from the migration data:**
- Poisson `vcov_model` diagonal entries for alpha: O(0.001)
- iOLS `vcov_model` diagonal entries for alpha: O(1) — three orders of magnitude larger
- This produces bias > $\hat{\xi}_g$ for UKR groups → negative bias-corrected estimate

## Theoretical Analysis

### What the bias correction actually needs

The Taylor expansion requires $\text{Var}(\hat{\alpha})$ — the **true sampling variance** of the estimator. Under correct specification:

$$\sqrt{n}(\hat{\alpha} - \alpha_0) \xrightarrow{d} \mathcal{N}(0, \mathbf{\Omega})$$

The bias correction should use an estimate of $\mathbf{\Omega}/n$, which is the finite-sample variance of $\hat{\alpha}$.

### Three candidates for $\mathbf{V}$

**1. Model-based $\mathbf{V}_{\text{model}}$ (current for Poisson/OLS):**

For Poisson: $\mathbf{V}_{\text{model}} = (\mathbf{Z}'\text{diag}(\hat{\mu})\mathbf{Z})^{-1}$. This is the inverse Fisher information, valid under correct Poisson specification (equidispersion). It gives the right bias correction when $\text{Var}(m_i | \mathbf{x}_i) = \mu_i$.

For iOLS: $\mathbf{V}_{\text{model}} = \sigma^2 (\mathbf{Z}'\mathbf{Z})^{-1}$ where $\sigma^2$ is the average squared GPML residual. This is NOT the inverse Fisher information for GPML. The correct GPML Fisher information is $(\mathbf{Z}'\mathbf{Z})^{-1}$ (unit working weight), NOT scaled by $\sigma^2$.

**2. Sandwich $\mathbf{V}_{\text{HC}}$ (robust):**

$$\mathbf{V}_{\text{HC}} = (\mathbf{Z}'\mathbf{Z})^{-1} \mathbf{Z}' \text{diag}(e_i^2) \mathbf{Z} (\mathbf{Z}'\mathbf{Z})^{-1}$$

This is consistent for $\text{Var}(\hat{\alpha})$ under heteroskedasticity. However, for Poisson data it tends to **overestimate** the variance (especially with HC3 leverage corrections), leading to overcorrection of the bias.

The original package documentation notes: *"Model-based variance is used for bias correction rather than HC-robust variance, because HC3 can be inflated by high-leverage observations."*

**3. GPML-specific Fisher information $(\mathbf{Z}'\mathbf{Z})^{-1}$:**

For the Gamma family with log link, the variance function is $V(\mu) = \mu^2$, giving working weights $w_i = 1$ (i.e., all observations contribute equally). The expected Fisher information is simply $\mathbf{Z}'\mathbf{Z}$, and the model-based variance is:

$$\mathbf{V}_{\text{GPML}} = (\mathbf{Z}'\mathbf{Z})^{-1}$$

**without any $\sigma^2$ scaling.** This is because the Gamma quasi-likelihood automatically normalizes the variance — the dispersion parameter is implicit.

## Recommended Fix

### For iOLS specifically

The correct model-based variance for GPML is $(\mathbf{Z}'\mathbf{Z})^{-1}$, **not** $\sigma^2 (\mathbf{Z}'\mathbf{Z})^{-1}$. The $\sigma^2$ scaling is appropriate for OLS/NLS (Gaussian quasi-likelihood) but NOT for GPML (Gamma quasi-likelihood).

**Implementation:** In `.fit_iols()`, compute `vcov_model` as:
```r
V_model <- .solve_safe(crossprod(Z))  # (Z'Z)^{-1}, no sigma2 scaling
```
instead of:
```r
V_model <- .compute_model_vcov(Z, weights, sigma2 = sigma2)  # sigma2 * (Z'Z)^{-1}
```

### Justification

The GPML score is $s_i(\beta) = (m_i/\mu_i - 1) \mathbf{z}_i$. The Fisher information is:
$$\mathcal{I}(\beta) = E\left[\sum_i s_i s_i'\right] = \sum_i E[(m_i/\mu_i - 1)^2] \mathbf{z}_i \mathbf{z}_i'$$

Under the Gamma model, $\text{Var}(m_i/\mu_i) = 1$ (the coefficient of variation is constant), so:
$$\mathcal{I}(\beta) = \sum_i \mathbf{z}_i \mathbf{z}_i' = \mathbf{Z}'\mathbf{Z}$$

The inverse $(\mathbf{Z}'\mathbf{Z})^{-1}$ is the correct model-based variance. The empirical $\hat{\sigma}^2 = \sum (m_i/\hat{\mu}_i - 1)^2 / (n-p)$ can be much larger than 1 (especially with count data that doesn't follow a Gamma distribution), which is why scaling by $\hat{\sigma}^2$ inflates the variance for bias correction.

### Why not use the sandwich for bias correction?

The sandwich $\mathbf{V}_{\text{HC}}$ is designed for **inference** (CIs, p-values) under misspecification. The bias correction is a different problem: it requires $E[\hat{\alpha} - \alpha_0]$, not $\text{Var}(\hat{\alpha})$ under misspecification. The Taylor expansion:

$$E[\hat{\xi}] \approx \xi(\alpha_0) + \frac{1}{2} \text{tr}(H \cdot \text{Var}(\hat{\alpha}))$$

uses $\text{Var}(\hat{\alpha})$ under the **assumed model** (because the bias arises from the nonlinearity of $\xi$, not from model misspecification). Using the sandwich here would "correct" for both the nonlinearity bias AND apparent heteroskedasticity, which is double-counting — the sandwich variance includes variance from model misspecification that has nothing to do with the Jensen's-inequality bias.

### For Poisson and NB

No change needed. Poisson uses the inverse Fisher information $(\mathbf{Z}'\text{diag}(\hat{\mu})\mathbf{Z})^{-1}$, which is the correct model-based variance. NB similarly uses $(\mathbf{Z}'\text{diag}(w_i^{NB})\mathbf{Z})^{-1}$.

## Implementation Details

### File: `R/estimators.R`, `.fit_iols()` function

Change (around line 754):
```r
# Current (incorrect for GPML):
sigma2 <- sum(weights * gpml_resid^2) / (n_obs - p)
V_model <- .compute_model_vcov(Z, weights, sigma2 = sigma2)

# Correct (GPML Fisher information, no sigma2):
V_model <- .compute_model_vcov(Z, weights, sigma2 = NULL)
```

The `.compute_model_vcov()` function already handles `sigma2 = NULL` by returning `(Z'WZ)^{-1}` without scaling.

### Tests to verify

1. **iOLS bias correction is positive for UKR groups**: With the fix, `est_bc > 0` for all groups on the migration data with `year*ukr` covariates.

2. **iOLS bias correction magnitude comparable to Poisson**: The ratio `bias / est` should be similar order of magnitude for Poisson and iOLS on the same data.

3. **Regression test on simulated data**: With DGP alpha=0.70, verify that all methods produce `est_bc > 0` and `est_bc < est`.

## Summary

| Method | Correct $\mathbf{V}_{\text{model}}$ for bias correction |
|--------|--------------------------------------------------------|
| OLS | $\hat{\sigma}^2 (\mathbf{Z}'\mathbf{Z})^{-1}$ — Gaussian quasi-likelihood |
| Poisson | $(\mathbf{Z}'\text{diag}(\hat{\mu})\mathbf{Z})^{-1}$ — Poisson Fisher information |
| NB | $(\mathbf{Z}'\text{diag}(w_i^{NB})\mathbf{Z})^{-1}$ — NB Fisher information |
| **iOLS** | $(\mathbf{Z}'\mathbf{Z})^{-1}$ — **Gamma Fisher information (no $\sigma^2$ scaling)** |

The one-line fix: remove $\sigma^2$ from `V_model` for iOLS.
