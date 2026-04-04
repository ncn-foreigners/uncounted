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

test_that("plot_explore with constrained model runs without error", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, constrained = TRUE)
  expect_no_error(plot_explore(fit))
})
