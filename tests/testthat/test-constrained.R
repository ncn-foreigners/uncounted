# ---- Test data (shared) ----

make_test_data <- function() {
  data.frame(
    m = c(14L, 1L, 56L, 21L, 29L, 57L, 1L, 41L, 258L, 6L,
          19L, 9L, 3L, 124L, 4L, 15L, 30L, 10L, 7L, 686L),
    n = c(2L, 1L, 27L, 9L, 85L, 82L, 1L, 16L, 586L, 2L,
          38L, 1L, 1L, 37L, 4L, 12L, 19L, 6L, 7L, 774L),
    N = c(1023L, 188L, 14655L, 892L, 3672L, 4473L, 44L, 5081L, 131365L, 1357L,
          4391L, 320L, 309L, 6422L, 482L, 1119L, 1784L, 746L, 2732L, 25576L),
    sex = rep(c("M", "F"), each = 10),
    country = paste0("C", 1:20)
  )
}

# ---- Constrained Poisson ----

test_that("constrained Poisson returns alpha in (0,1) and beta > 0", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                             gamma = 0.001, constrained = TRUE)
  expect_s3_class(fit, "uncounted")
  expect_true(isTRUE(fit$constrained))
  expect_true(all(fit$alpha_values > 0 & fit$alpha_values < 1))
  expect_true(all(fit$beta_values > 0))
})

# ---- Constrained NB ----

test_that("constrained NB returns alpha in (0,1) and beta > 0", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "nb",
                             gamma = 0.001, constrained = TRUE)
  expect_true(all(fit$alpha_values > 0 & fit$alpha_values < 1))
  expect_true(all(fit$beta_values > 0))
})

# ---- Constrained popsize: ordered CI and reasonable BC ----

test_that("constrained popsize gives ordered CI and positive BC", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                             gamma = 0.001, constrained = TRUE)
  ps <- popsize(fit, bias_correction = TRUE)
  expect_true(all(ps$lower < ps$estimate))
  expect_true(all(ps$estimate < ps$upper))
  expect_true(all(ps$estimate_bc > 0))
})

# ---- Constrained BC is not hugely negative (regression for delta method fix) ----

test_that("constrained BC does not produce hugely negative values", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                             gamma = 0.001, constrained = TRUE)
  ps <- popsize(fit, bias_correction = TRUE)
  # BC should reduce by less than 50% — anything more signals broken delta method
  expect_true(all(ps$estimate_bc > 0.5 * ps$estimate))
})

# ---- Constrained with covariates ----

test_that("constrained with cov_alpha works", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                             gamma = 0.001, constrained = TRUE,
                             cov_alpha = ~ 0 + sex)
  ps <- popsize(fit)
  expect_equal(nrow(ps), 2)
  expect_true(all(fit$alpha_values > 0 & fit$alpha_values < 1))
})

# ---- Print does not duplicate Constrained line ----

test_that("constrained print has single Constrained line", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                             gamma = 0.001, constrained = TRUE)
  out <- capture.output(print(fit))
  constrained_lines <- grep("^Constrained:", out)
  expect_equal(length(constrained_lines), 1)
})

# ---- Link-scale coefficients match inv_logit transform ----

test_that("constrained alpha_values = inv_logit(X %*% alpha_coefs)", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                             gamma = 0.001, constrained = TRUE)
  eta <- as.numeric(fit$X_alpha %*% fit$alpha_coefs)
  expected_alpha <- 1 / (1 + exp(-eta))
  expect_equal(fit$alpha_values, expected_alpha, tolerance = 1e-10)
})
