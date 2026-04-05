# ---- Diagnostics tests (CRAN) ----

test_that("logLik returns correct class and attributes for Poisson", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(!is.null(attr(ll, "df")))
  expect_true(!is.null(attr(ll, "nobs")))
  expect_equal(attr(ll, "nobs"), fit$n_obs)
})

test_that("logLik returns correct class for NB", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005)
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(attr(ll, "df") > length(coef(fit)))
})

test_that("AIC and BIC work on uncounted objects", {
  d <- small_data()
  fit_po <- quick_fit(d, method = "poisson", gamma = 0.005)
  fit_nb <- quick_fit(d, method = "nb", gamma = 0.005)
  expect_true(is.finite(AIC(fit_po)))
  expect_true(is.finite(BIC(fit_po)))
  expect_true(is.finite(AIC(fit_nb)))
})

test_that("logLik pseudo-loglik for OLS is finite", {
  d <- positive_data(small_data())
  fit <- quick_fit(d, method = "ols", gamma = NULL)
  ll <- logLik(fit)
  expect_true(is.finite(as.numeric(ll)))
})

test_that("nobs returns correct count", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  expect_equal(nobs(fit), nrow(d))
})

test_that("deviance is positive for Poisson and NB", {
  d <- small_data()
  for (method in c("poisson", "nb")) {
    fit <- quick_fit(d, method = method, gamma = 0.005)
    expect_true(deviance(fit) > 0)
  }
})

test_that("Pearson residuals are finite", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  r_p <- residuals(fit, type = "pearson")
  expect_length(r_p, nrow(d))
  expect_true(all(is.finite(r_p)))
})

test_that("deviance residuals are finite", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  r_d <- residuals(fit, type = "deviance")
  expect_true(all(is.finite(r_d)))
})

test_that("Anscombe residuals are finite for Poisson and NB", {
  d <- small_data()
  for (method in c("poisson", "nb")) {
    fit <- quick_fit(d, method = method, gamma = 0.005)
    r_a <- residuals(fit, type = "anscombe")
    expect_true(all(is.finite(r_a)))
  }
})

test_that("response residuals match m - fitted", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  r <- residuals(fit, type = "response")
  expect_equal(r, fit$m - fit$fitted.values)
})

# ---- Bug regressions ----

test_that("NB deviance is non-negative even with many zeros", {
  d <- small_data()
  d$m[1:5] <- 0L
  fit <- quick_fit(d, method = "nb", gamma = 0.005)
  expect_true(deviance(fit) >= 0)
  r_d <- residuals(fit, type = "deviance")
  expect_true(all(is.finite(r_d)))
  zero_idx <- which(fit$m == 0 & fit$fitted.values > 0)
  if (length(zero_idx) > 0) {
    expect_true(all(r_d[zero_idx] <= 0))
  }
})

test_that("OLS with all positive m uses log(m) not log(m+1)", {
  d <- positive_data(small_data())
  fit <- quick_fit(d, method = "ols", gamma = NULL)

  Z <- cbind(log(d$N), log(d$n / d$N))
  manual_fit <- lm.fit(Z, log(d$m))
  manual_sigma2 <- sum(manual_fit$residuals^2) / nrow(d)
  manual_ll <- -nrow(d) / 2 * (log(2 * pi * manual_sigma2) + 1)

  expect_equal(as.numeric(logLik(fit)), manual_ll, tolerance = 0.1)
})

test_that("OLS gamma profiling uses log(m) when no zeros", {
  d <- positive_data(small_data())
  fit <- quick_fit(d, method = "ols", gamma = "estimate")
  log_m <- log(d$m)
  resid_manual <- log_m - fit$log_mu
  expect_true(abs(mean(resid_manual)) < 0.5)
})

test_that("cov_gamma parameter works with estimated gamma", {
  d <- small_data()
  d$sex2 <- d$sex
  fit <- quick_fit(d, gamma = "estimate", cov_gamma = ~sex2)
  expect_s3_class(fit, "uncounted")
  expect_true(!is.null(fit$gamma_coefs))
})

# ---- Regression: OLS zero handling scale consistency ----

test_that("OLS with zeros: fitted values are exp(log_mu)", {
  d <- small_data()
  d$m[1] <- 0L
  fit <- quick_fit(d, method = "ols", gamma = 0.005)
  expect_equal(fit$fitted.values, exp(fit$log_mu), tolerance = 1e-10)
  expect_equal(residuals(fit, type = "response"), d$m - fit$fitted.values)
})
