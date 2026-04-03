make_test_data <- function(n = 50, seed = 42) {
  set.seed(seed)
  N <- round(exp(rnorm(n, 5, 1.5)))
  n_aux <- rpois(n, exp(-3 + 0.8 * log(N)))
  ratio <- n_aux / N
  mu <- N^0.7 * (0.005 + ratio)^0.5
  m <- rpois(n, mu)
  data.frame(m = m, n = n_aux, N = N)
}

test_that("plot.uncounted runs without error for Poisson", {
  dt <- make_test_data()
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  expect_no_error(plot(fit, which = 1:4, ask = FALSE))
})

test_that("plot.uncounted runs without error for NB", {
  dt <- make_test_data()
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  expect_no_error(plot(fit, which = 1:4, ask = FALSE))
})

test_that("plot.uncounted runs for single panel", {
  dt <- make_test_data()
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  expect_no_error(plot(fit, which = 1, ask = FALSE))
  expect_no_error(plot(fit, which = 4, ask = FALSE))
})

test_that("plot_explore runs without error", {
  dt <- make_test_data()
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  expect_no_error(plot_explore(fit))
})

test_that("rootogram works for Poisson", {
  dt <- make_test_data()
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  result <- rootogram(fit)
  expect_true(sum(result$observed) == nrow(dt))
  expect_true(abs(sum(result$expected) - nrow(dt)) < 1)  # should sum to ~n
})

test_that("rootogram works for NB", {
  dt <- make_test_data()
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  result <- rootogram(fit)
  expect_true(sum(result$observed) == nrow(dt))
})

test_that("rootogram errors for OLS", {
  dt <- make_test_data()
  dt_pos <- dt[dt$m > 0 & dt$n > 0, ]
  fit <- estimate_hidden_pop(dt_pos, ~m, ~n, ~N, method = "ols", gamma = NULL)
  expect_error(rootogram(fit), "Poisson and NB")
})

test_that("rootogram styles work", {
  dt <- make_test_data()
  fit <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  expect_no_error(rootogram(fit, style = "standing"))
  expect_no_error(rootogram(fit, style = "suspended"))
})
