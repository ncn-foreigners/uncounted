# ---- Plotting tests (CRAN) ----

test_that("plot.uncounted runs without error for Poisson", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  expect_no_error(plot(fit, which = 1:4, ask = FALSE))
})

test_that("plot.uncounted runs without error for NB", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005)
  expect_no_error(plot(fit, which = 1:4, ask = FALSE))
})

test_that("plot.uncounted runs for single panel", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  expect_no_error(plot(fit, which = 1, ask = FALSE))
  expect_no_error(plot(fit, which = 4, ask = FALSE))
})

test_that("plot_explore runs without error", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  expect_no_error(plot_explore(fit))
})

test_that("rootogram works for Poisson", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  result <- rootogram(fit)
  expect_true(sum(result$observed) == nrow(d))
  expect_true(abs(sum(result$expected) - nrow(d)) < 1)
})

test_that("rootogram works for NB", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005)
  result <- rootogram(fit)
  expect_true(sum(result$observed) == nrow(d))
})

test_that("rootogram errors for OLS", {
  d <- positive_data(small_data())
  fit <- quick_fit(d, method = "ols", gamma = NULL)
  expect_error(rootogram(fit), "Poisson and NB")
})

test_that("rootogram styles work", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  expect_no_error(rootogram(fit, style = "standing"))
  expect_no_error(rootogram(fit, style = "suspended"))
})

test_that("rootogram returns capped bins and raw-scale frequencies", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  result <- rootogram(fit, style = "suspended", max_count = 1, sqrt_scale = FALSE)

  expect_identical(result$counts, 0:1)
  expect_false(result$sqrt_scale)
  expect_equal(sum(result$observed), nrow(d))
  expect_true(all(result$expected >= 0))
})

test_that("profile_gamma with plot=TRUE runs without error", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  expect_no_error(profile_gamma(fit, gamma_grid = c(0.001, 0.01, 0.1), plot = TRUE))
})

test_that("profile_gamma with estimated gamma marks estimate on plot", {
  d <- small_data()
  fit <- quick_fit(d, gamma = "estimate")
  expect_no_error(profile_gamma(fit, gamma_grid = c(0.001, 0.01, 0.1), plot = TRUE))
})

test_that("profile helpers return compact data frames without plotting", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)

  gamma_prof <- profile_gamma(fit, gamma_grid = c(0.001, 0.01, 0.1), plot = FALSE)
  alpha_grid <- fit$alpha_coefs[1] + c(-0.05, 0, 0.05)
  beta_grid <- fit$beta_coefs[1] + c(-0.05, 0, 0.05)
  alpha_prof <- profile_alpha(fit, grid = alpha_grid, plot = FALSE)
  beta_prof <- profile_beta(fit, grid = beta_grid, plot = FALSE)
  beta_dispatch <- profile(fit, param = "beta", grid = beta_grid, plot = FALSE)

  expect_identical(names(gamma_prof), c("gamma", "xi", "loglik"))
  expect_identical(names(alpha_prof), c("value", "xi", "loglik"))
  expect_identical(names(beta_prof), c("value", "xi", "loglik"))
  expect_equal(beta_dispatch, beta_prof)
})

test_that("profile_gamma rejects cov_gamma fits", {
  d <- cov_gamma_data()
  fit <- suppressWarnings(
    estimate_hidden_pop(
      data = d,
      observed = ~ m,
      auxiliary = ~ n,
      reference_pop = ~ N,
      method = "poisson",
      gamma = "estimate",
      cov_gamma = ~ 0 + z
    )
  )

  expect_error(profile_gamma(fit, plot = FALSE), "covariate-varying gamma")
})

test_that("plot_explore with constrained model runs without error", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, constrained = TRUE)
  expect_no_error(plot_explore(fit))
})

test_that(".running_mean covers every point", {
  ## w = n %/% 10 used to stop at the 10th of 19 points, and gave NA windows for
  ## n < 10
  x <- 1:19
  rm19 <- .running_mean(x, x)
  expect_length(rm19$x, 10L)
  expect_equal(max(rm19$x), 19)
  expect_equal(weighted.mean(rm19$y, tabulate(ceiling(seq_len(19) * 10 / 19))), mean(x))

  rm5 <- .running_mean(1:5, 5:1)
  expect_equal(rm5$x, 1:5)
  expect_equal(rm5$y, 5:1)
  expect_false(anyNA(rm5$x))

  ## Sorted by x, whatever the input order
  expect_equal(.running_mean(c(3, 1, 2), c(30, 10, 20))$y, c(10, 20, 30))
})
