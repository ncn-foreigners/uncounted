# ---- NLS fit statistics use the scale NLS optimises (CRAN) ----
#
# logLik and deviance used the log-scale Gaussian that OLS maximises, so OLS
# always "won" compare_models() even when NLS had the smaller raw-scale RSS.

nls_ols <- function() {
  d <- positive_data()
  list(d = d,
       ols = estimate_hidden_pop(d, ~m, ~n, ~N, method = "ols", gamma = 0.005),
       nls = estimate_hidden_pop(d, ~m, ~n, ~N, method = "nls", gamma = 0.005))
}

test_that("NLS logLik and deviance are raw-scale", {
  x <- nls_ols()
  rss <- sum((x$d$m - fitted(x$nls))^2)
  n <- nrow(x$d)
  expect_equal(deviance(x$nls), rss)
  expect_equal(as.numeric(logLik(x$nls)), -n / 2 * (log(2 * pi * rss / n) + 1))
})

test_that("NLS logLik is the maximum of its own criterion", {
  d <- testdata
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "nls", gamma = 0.005)
  ref <- nls(m ~ N^a * (0.005 + n / N)^b, data = d, start = list(a = 0.7, b = 0.5))
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = 1e-5)
})

test_that("OLS logLik stays on the log scale", {
  x <- nls_ols()
  r <- log(x$d$m) - x$ols$log_mu
  n <- nrow(x$d)
  expect_equal(as.numeric(logLik(x$ols)), -n / 2 * (log(2 * pi * sum(r^2) / n) + 1))
})

test_that("compare_models() warns when mixing pseudo-likelihood families", {
  x <- nls_ols()
  expect_warning(compare_models(x$ols, x$nls), "not comparable")
  iols <- estimate_hidden_pop(x$d, ~m, ~n, ~N, method = "iols", gamma = 0.005)
  expect_warning(compare_models(x$ols, iols), "not comparable")
  ## Same family: no warning about comparability
  o2 <- estimate_hidden_pop(x$d, ~m, ~n, ~N, method = "ols", cov_alpha = ~sex, gamma = 0.005)
  expect_no_warning(compare_models(x$ols, o2))
})

test_that("NLS Pearson and deviance residuals use the raw scale", {
  fit <- nls_ols()$nls
  r <- fit$m - fitted(fit)
  expect_equal(unname(residuals(fit, type = "deviance")), unname(r / sd(r)))
  p <- length(coef(fit))
  expect_equal(unname(residuals(fit, type = "pearson")),
               unname(r / sqrt(sum(r^2) / (fit$n_obs - p))))
})
