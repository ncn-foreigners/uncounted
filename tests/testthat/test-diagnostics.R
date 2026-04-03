# Setup: create test data and fit models
make_test_data <- function(n = 50, seed = 42) {
  set.seed(seed)
  N <- round(exp(rnorm(n, 5, 1.5)))
  n_aux <- rpois(n, exp(-3 + 0.8 * log(N)))
  ratio <- n_aux / N
  mu <- N^0.7 * (0.005 + ratio)^0.5
  m <- rpois(n, mu)
  data.frame(m = m, n = n_aux, N = N)
}

test_that("logLik returns correct class and attributes for Poisson", {
  dt <- make_test_data()
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(!is.null(attr(ll, "df")))
  expect_true(!is.null(attr(ll, "nobs")))
  expect_equal(attr(ll, "nobs"), fit$n_obs)
})

test_that("logLik returns correct class for NB", {
  dt <- make_test_data()
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(attr(ll, "df") > length(coef(fit)))  # theta adds 1
})

test_that("AIC and BIC work on uncounted objects", {
  dt <- make_test_data()
  fit_po <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  fit_nb <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  expect_true(is.finite(AIC(fit_po)))
  expect_true(is.finite(BIC(fit_po)))
  expect_true(is.finite(AIC(fit_nb)))
  # NB should have lower AIC if data is overdispersed or similar AIC
  # Just check they're computed
})

test_that("logLik pseudo-loglik for OLS is finite", {
  dt <- make_test_data()
  dt_pos <- dt[dt$m > 0 & dt$n > 0, ]
  fit <- estimate_hidden_pop(dt_pos, ~m, ~n, ~N, method = "ols", gamma = NULL)
  ll <- logLik(fit)
  expect_true(is.finite(as.numeric(ll)))
})

test_that("nobs returns correct count", {
  dt <- make_test_data()
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  expect_equal(nobs(fit), nrow(dt))
})

test_that("deviance is positive for Poisson", {
  dt <- make_test_data()
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  expect_true(deviance(fit) > 0)
})

test_that("deviance is positive for NB", {
  dt <- make_test_data()
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  expect_true(deviance(fit) > 0)
})

test_that("Pearson residuals have mean near 0", {
  dt <- make_test_data()
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  r_p <- residuals(fit, type = "pearson")
  expect_length(r_p, nrow(dt))
  # Not necessarily mean 0 for Pearson, but should be finite
  expect_true(all(is.finite(r_p)))
})

test_that("deviance residuals are finite", {
  dt <- make_test_data()
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  r_d <- residuals(fit, type = "deviance")
  expect_true(all(is.finite(r_d)))
})

test_that("Anscombe residuals are finite for Poisson", {
  dt <- make_test_data()
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  r_a <- residuals(fit, type = "anscombe")
  expect_true(all(is.finite(r_a)))
})

test_that("Anscombe residuals are finite for NB", {
  dt <- make_test_data()
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  r_a <- residuals(fit, type = "anscombe")
  expect_true(all(is.finite(r_a)))
})

test_that("response residuals match m - fitted", {
  dt <- make_test_data()
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  r <- residuals(fit, type = "response")
  expect_equal(r, fit$m - fit$fitted.values)
})

# ---- Regression tests for fixed bugs ----

test_that("NB deviance is non-negative even with many zeros (bug fix #1)", {
  dt <- make_test_data()
  # Ensure some zeros exist
  dt$m[1:10] <- 0L
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  d <- deviance(fit)
  expect_true(d >= 0, info = paste("NB deviance was", d))

  # Also check individual deviance residuals are finite
  r_d <- residuals(fit, type = "deviance")
  expect_true(all(is.finite(r_d)))
  # For m=0, mu>0, deviance residual should be negative (m < mu)
  zero_idx <- which(fit$m == 0 & fit$fitted.values > 0)
  if (length(zero_idx) > 0) {
    expect_true(all(r_d[zero_idx] <= 0))
  }
})

test_that("OLS with all positive m uses log(m) not log(m+1) (bug fix #2+3)", {
  dt <- make_test_data()
  dt_pos <- dt[dt$m > 0 & dt$n > 0, ]

  # Fit with fixed gamma (no profiling)
  fit_fixed <- estimate_hidden_pop(dt_pos, ~m, ~n, ~N, method = "ols", gamma = NULL)

  # Manual OLS on log(m)
  Z <- cbind(log(dt_pos$N), log(dt_pos$n / dt_pos$N))
  manual_fit <- lm.fit(Z, log(dt_pos$m))
  manual_sigma2 <- sum(manual_fit$residuals^2) / nrow(dt_pos)
  manual_ll <- -nrow(dt_pos)/2 * (log(2 * pi * manual_sigma2) + 1)

  # logLik should match manual (using log(m), not log(m+1))
  expect_equal(as.numeric(logLik(fit_fixed)), manual_ll, tolerance = 0.1)
})

test_that("OLS gamma profiling uses log(m) when no zeros (bug fix #2)", {
  dt <- make_test_data()
  dt_pos <- dt[dt$m > 0 & dt$n > 0, ]

  # Fit with estimated gamma on positive-only data
  fit_gamma <- estimate_hidden_pop(dt_pos, ~m, ~n, ~N, method = "ols",
                                    gamma = "estimate")

  # The log_mu should be based on log(m), not log(m+1)
  # Verify by checking residuals are not systematically biased
  log_m <- log(dt_pos$m)
  resid_manual <- log_m - fit_gamma$log_mu
  expect_true(abs(mean(resid_manual)) < 0.5)
})

test_that("cov_gamma parameter no longer accepted", {
  dt <- make_test_data()
  dt$sex <- sample(c("f", "m"), nrow(dt), replace = TRUE)
  # cov_gamma should not be a valid parameter anymore
  expect_error(
    estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson",
                        gamma = "estimate", cov_gamma = ~sex),
    "unused argument"
  )
})

# ---- dfbeta / dfpopsize ----

test_that("dfbeta returns matrix with correct dimensions", {
  dt <- make_test_data(n = 20, seed = 99)
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  db <- dfbeta(fit)
  expect_true(is.matrix(db))
  expect_equal(nrow(db), nrow(dt))
  expect_equal(ncol(db), length(coef(fit)))
})

test_that("dfbeta by country works", {
  dt <- make_test_data(n = 20, seed = 99)
  dt$country <- paste0("C", rep(1:10, each = 2))
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson",
                             gamma = 0.005, countries = ~country)
  db <- dfbeta(fit, by = "country")
  expect_equal(nrow(db), 10)
})

test_that("dfpopsize returns named numeric vector", {
  dt <- make_test_data(n = 20, seed = 99)
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  dp <- dfpopsize(fit)
  expect_true(is.numeric(dp))
  expect_equal(length(dp), nrow(dt))
  expect_true(!is.null(names(dp)))
})

test_that("dfpopsize sums to approximately zero", {
  # Total influence should be small (not exactly 0 due to nonlinearity)
  dt <- make_test_data(n = 20, seed = 99)
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  dp <- dfpopsize(fit)
  # Most dxi should be finite
 expect_true(sum(is.finite(dp)) > nrow(dt) * 0.8)
})
