# ---- Constrained estimation on full data (extended, skip on CRAN) ----

skip_on_cran()

test_that("constrained Poisson on full data: alpha in (0,1)", {
  skip_on_cran()
  fit <- quick_fit(testdata, gamma = 0.005, constrained = TRUE)
  expect_true(all(fit$alpha_values > 0 & fit$alpha_values < 1))
  expect_true(all(fit$beta_values > 0))
})

test_that("constrained NB on full data", {
  skip_on_cran()
  fit <- quick_fit(testdata, method = "nb", gamma = 0.005, constrained = TRUE)
  expect_true(all(fit$alpha_values > 0 & fit$alpha_values < 1))
  expect_true(all(fit$beta_values > 0))
  expect_true(fit$theta > 0)
})

test_that("constrained with covariates on full data", {
  skip_on_cran()
  fit <- quick_fit(testdata, gamma = 0.005, constrained = TRUE,
                   cov_alpha = ~0 + sex)
  ps <- popsize(fit)
  expect_equal(nrow(ps), 2)
  expect_true(all(fit$alpha_values > 0 & fit$alpha_values < 1))
})

test_that("constrained with estimated gamma on full data", {
  skip_on_cran()
  fit <- quick_fit(testdata, gamma = "estimate", constrained = TRUE)
  expect_true(fit$gamma > 0)
  expect_true(all(fit$alpha_values > 0 & fit$alpha_values < 1))
})

test_that("constrained popsize on full data: ordered CI", {
  skip_on_cran()
  fit <- quick_fit(testdata, gamma = 0.005, constrained = TRUE)
  ps <- popsize(fit, bias_correction = TRUE)
  expect_true(all(ps$lower < ps$estimate))
  expect_true(all(ps$estimate < ps$upper))
  expect_true(all(ps$estimate_bc > 0))
})
