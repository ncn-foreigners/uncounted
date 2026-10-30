# ---- Non-integer observed values (CRAN) ----
#
# dpois() returns -Inf for non-integers, so the objective used to be the constant
# 1e20 and the fit returned the OLS start values with convergence = 0.

noninteger_data <- function() {
  d <- testdata
  d$mx <- d$m * 1.0001 + 0.3
  d
}

test_that("Poisson with non-integer m matches the quasi-Poisson GLM", {
  d <- noninteger_data()
  fit <- estimate_hidden_pop(d, ~mx, ~n, ~N, gamma = 0.005)
  ref <- glm(mx ~ 0 + log(N) + log(0.005 + n / N), family = quasipoisson(), data = d)
  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-5)
  expect_true(fit$pseudo_loglik)
  expect_true(is.finite(fit$loglik) && fit$loglik > -1e10)
})

test_that("the Poisson kernel equals dpois() for counts", {
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005)
  expect_false(fit$pseudo_loglik)
  expect_equal(fit$loglik, sum(dpois(testdata$m, fitted(fit), log = TRUE)))
})

test_that("NB with non-integer m gives a clear error, for every estimator", {
  d <- noninteger_data()
  expect_error(estimate_hidden_pop(d, ~mx, ~n, ~N, method = "nb", gamma = 0.005),
               "whole-number")
  expect_error(estimate_hidden_pop(d, ~mx, ~n, ~N, method = "nb", gamma = 0.005,
                                   estimator = "gmm"),
               "whole-number")
})

test_that("pseudo-likelihood fits warn in comparisons and rootogram errors", {
  d <- noninteger_data()
  f0 <- estimate_hidden_pop(d, ~mx, ~n, ~N, gamma = 0.005)
  f1 <- estimate_hidden_pop(d, ~mx, ~n, ~N, cov_alpha = ~sex, gamma = 0.005)
  expect_warning(compare_models(f0, f1), "pseudo-likelihood")
  expect_warning(lrtest(f0, f1), "pseudo-likelihood")
  expect_error(rootogram(f0), "whole-number")
})
