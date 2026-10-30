# ---- Starting value and convergence for an estimated gamma (CRAN) ----
#
# Before the profile screen, gamma = "estimate" always started at
# mean(gamma_bounds) = 0.25 and, with covariates, L-BFGS-B stopped at
# gamma ~ 0.29 with a log-likelihood ~14 below the optimum while reporting
# convergence = 0 (xi about 2.5x too small).

fixed_gamma_loglik <- function(..., gamma) {
  estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = gamma, ...)$loglik
}

test_that("Poisson with cov_alpha reaches the optimum, not the gamma ~ 0.29 mode", {
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~year)
  expect_equal(fit$convergence, 0L)
  expect_gt(fit$gamma, 0.002)
  expect_lt(fit$gamma, 0.01)
  ## No fixed gamma does better than the joint optimum
  ll_fixed <- vapply(c(0.003, 0.005, 0.008, 0.29),
                     function(g) fixed_gamma_loglik(cov_alpha = ~year, gamma = g),
                     numeric(1))
  expect_gte(fit$loglik, max(ll_fixed) - 1e-6)
})

test_that("NB with cov_beta and weighted Poisson reach the optimum", {
  nb <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", cov_beta = ~year)
  expect_gt(nb$gamma, 0.002)
  expect_lt(nb$gamma, 0.01)
  expect_gte(nb$loglik, estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb",
                                            cov_beta = ~year, gamma = 0.005)$loglik - 1e-6)

  set.seed(1)
  w <- runif(nrow(testdata), 0.5, 2)
  pw <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~sex, weights = w)
  expect_lt(pw$gamma, 0.01)
  expect_gte(pw$loglik, estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~sex,
                                            weights = w, gamma = 0.005)$loglik - 1e-6)
})

test_that("the old start (0.25) now also reaches the optimum", {
  screened <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~year)
  old <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~year, gamma_start = 0.25)
  expect_equal(old$loglik, screened$loglik, tolerance = 1e-6)
  expect_equal(unname(coef(old)), unname(coef(screened)), tolerance = 1e-4)
})

test_that("optim_info records the start, the screen and the optimiser trace", {
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N)
  oi <- fit$optim_info
  expect_identical(oi$gamma_start_source, "screen")
  expect_length(oi$gamma_screen$grid, 25L)
  expect_equal(oi$gamma_start, oi$gamma_screen$grid[which.min(oi$gamma_screen$deviance)])
  expect_true(all(c("restarts", "max_grad", "message", "gamma_at_bound") %in% names(oi)))
  expect_identical(oi$gamma_at_bound, "none")

  user <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma_start = 0.01)
  expect_identical(user$optim_info$gamma_start_source, "user")
  expect_null(user$optim_info$gamma_screen)
})

test_that("gamma at a bound gives a message and is recorded", {
  ## From gamma_start = 0.4 the optimiser ends in the genuine local optimum
  ## on the upper bound.
  expect_message(
    fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~year, gamma_start = 0.4),
    "upper bound"
  )
  expect_identical(fit$optim_info$gamma_at_bound, "upper")
})

test_that("gamma_start is validated and ignored for fixed gamma", {
  expect_error(estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma_start = 0.9),
               "gamma_start")
  expect_error(estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma_start = c(0.1, 0.2)),
               "gamma_start")
  fixed <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005, gamma_start = 0.3)
  expect_identical(fixed$gamma, 0.005)
})

test_that("refits warm-start from the original gamma", {
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N)
  args <- .refit_args_from_object(fit)
  expect_equal(args$gamma_start, fit$gamma)
  expect_null(.refit_args_from_object(quick_fit(gamma = 0.005))$gamma_start)
})

test_that("profile_dependence() is smooth in delta", {
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N)
  prof <- profile_dependence(fit, seq(-0.15, 0.15, by = 0.05), plot = FALSE)
  ## Before the fix xi jumped between ~67k and ~169k across this grid
  expect_true(all(diff(prof$xi) > 0))
  expect_lt(diff(range(prof$xi)) / prof$xi[4], 0.05)
})

test_that("refit loops do not repeat the gamma-bound message", {
  fit <- suppressMessages(
    estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~year, gamma_start = 0.4)
  )
  expect_identical(fit$optim_info$gamma_at_bound, "upper")
  expect_no_message(bootstrap_popsize(fit, R = 2, seed = 1, verbose = FALSE, ci_type = "cheap"))
  expect_no_message(profile_dependence(fit, c(-0.1, 0.1), plot = FALSE))
  small <- suppressMessages(
    estimate_hidden_pop(small_data(), ~m, ~n, ~N, gamma = "estimate", gamma_start = 0.5)
  )
  expect_no_message(loo(small))
})
