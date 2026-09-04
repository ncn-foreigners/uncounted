# ---- Helper coverage tests (CRAN) ----

test_that("link helpers are internally consistent across supported links", {
  eta <- c(-2, -0.5, 0, 0.5, 2)
  eps <- 1e-6

  expect_identical(uncounted:::.normalize_link_rho("logistic"), "logit")

  for (link_name in c("power", "cloglog", "logit", "probit")) {
    rho <- uncounted:::.rho_from_eta(eta, link_rho = link_name)
    log_rho <- uncounted:::.log_rho_from_eta(eta, link_rho = link_name)
    dlog <- uncounted:::.dlog_rho_deta(eta, link_rho = link_name)
    fd <- (
      uncounted:::.log_rho_from_eta(eta + eps, link_rho = link_name) -
      uncounted:::.log_rho_from_eta(eta - eps, link_rho = link_name)
    ) / (2 * eps)

    expect_equal(log(rho), log_rho, tolerance = 1e-8, info = link_name)
    expect_equal(dlog, fd, tolerance = 1e-5, info = link_name)
  }
})

test_that("rate and log-mu helpers match the manual formulas", {
  ratio <- c(0.1, 0.2, 0.3)
  gamma_values <- c(0.01, 0.02, 0.03)
  rate <- uncounted:::.rate_from_gamma(ratio, gamma_values)
  expect_equal(rate, ratio + gamma_values)
  expect_equal(
    uncounted:::.log_rate(ratio, gamma_values),
    log(ratio + gamma_values)
  )

  alpha_values <- c(0.5, 0.6, 0.7)
  beta_values <- c(0.8, 0.9, 1.0)
  log_N <- log(c(100, 200, 300))
  eta <- beta_values * log(rate)
  expect_equal(
    uncounted:::.compute_log_mu(alpha_values, log_N, beta_values, rate,
                                link_rho = "power"),
    alpha_values * log_N + eta
  )
})

test_that("model-matrix helper treats NULL and empty formulas as intercept-only", {
  d <- small_data()

  X_null <- uncounted:::.build_model_matrix(NULL, d)
  X_zero <- uncounted:::.build_model_matrix(~ 0, d)

  expect_equal(dim(X_null), c(nrow(d), 1))
  expect_equal(dim(X_zero), c(nrow(d), 1))
  expect_true(all(X_null[, 1] == 1))
  expect_true(all(X_zero[, 1] == 1))
  expect_identical(colnames(X_null), "(Intercept)")
  expect_identical(colnames(X_zero), "(Intercept)")
})

test_that("solve_safe falls back to qr.solve and generalized inverse", {
  A_rect <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
  b_rect <- c(1, 1)
  sol_rect <- uncounted:::.solve_safe(A_rect, b_rect)
  expect_length(sol_rect, ncol(A_rect))
  expect_true(all(is.finite(sol_rect)))

  A_sing <- matrix(1, nrow = 2, ncol = 2)
  b_sing <- c(1, 1)
  sol_sing <- uncounted:::.solve_safe(A_sing, b_sing)
  expect_true(all(is.finite(sol_sing)))
  expect_equal(as.numeric(A_sing %*% sol_sing), b_sing, tolerance = 1e-8)
})

test_that("model-based covariance helper matches sigma^2 (X'WX)^-1", {
  X <- cbind("(Intercept)" = 1, x = seq_len(6))
  weights <- seq(1, 6)
  V_model <- uncounted:::.compute_model_vcov(X, weights, sigma2 = 2)
  expect_equal(
    V_model,
    2 * solve(crossprod(X * sqrt(weights))),
    tolerance = 1e-10
  )
})

## Reference fit for the sandwich comparisons below. A weighted lm is used
## because sandwich::vcovHC()/vcovCL() then provide an independent oracle for
## every HC type, and .hat_values_wls() must reproduce hatvalues() exactly.
.sandwich_fixture <- function() {
  set.seed(11)
  n <- 30
  x <- rnorm(n)
  x[1] <- 25    # high leverage, so HC4/HC4m/HC5 exponents actually bind
  x[2] <- 14    # intermediate leverage, separates HC5's max(h) cut-off
  z <- runif(n)
  w <- runif(n, 0.5, 2)
  y <- 1 + 0.5 * x - z + rnorm(n)
  fit <- lm(y ~ x + z, weights = w)
  list(fit = fit, X = model.matrix(fit), e = residuals(fit), w = w, n = n)
}

test_that(".hat_values_wls reproduces hatvalues() for a weighted fit", {
  f <- .sandwich_fixture()
  expect_equal(
    unname(uncounted:::.hat_values_wls(f$X, f$w)),
    unname(stats::hatvalues(f$fit)),
    tolerance = 1e-10
  )
})

test_that(".compute_sandwich_vcov matches sandwich::vcovHC for HC0-HC5", {
  skip_if_not_installed("sandwich")
  f <- .sandwich_fixture()
  h <- uncounted:::.hat_values_wls(f$X, f$w)
  ## Prior weights enter the meat as a_i^2, so the residual carries sqrt(a_i).
  res <- sqrt(f$w) * f$e
  ## Guard the fixture: the high-leverage points must make HC5's cut-off bind,
  ## otherwise HC4/HC5 degenerate to the same exponent and the test is vacuous.
  p <- ncol(f$X)
  cutoff <- max(4, 0.7 * f$n * max(h) / p)
  expect_true(any(f$n * h / p > 4 & f$n * h / p < cutoff))

  for (type in c("HC0", "HC1", "HC2", "HC3", "HC4", "HC4m", "HC5")) {
    V_pkg <- uncounted:::.compute_sandwich_vcov(
      f$X, res, weights = f$w, hat_values = h, vcov_type = type
    )
    V_ref <- sandwich::vcovHC(f$fit, type = type)
    expect_equal(unname(V_pkg), unname(V_ref), tolerance = 1e-8,
                 info = paste("vcov_type =", type))
  }
})

test_that("cluster-robust branch is consistent with the observation-level one", {
  f <- .sandwich_fixture()
  h <- uncounted:::.hat_values_wls(f$X, f$w)
  res <- sqrt(f$w) * f$e

  ## One cluster per observation with CR0 (no small-sample correction) is
  ## exactly HC0; this is what catches a mis-scaled per-observation score.
  V_hc0 <- uncounted:::.compute_sandwich_vcov(
    f$X, res, weights = f$w, hat_values = h, vcov_type = "HC0"
  )
  V_cr0 <- uncounted:::.compute_sandwich_vcov(
    f$X, res, weights = f$w, hat_values = h, vcov_type = "HC0",
    cluster = seq_len(f$n)
  )
  expect_equal(V_cr0, V_hc0, tolerance = 1e-10)
})

test_that("cluster-robust branch matches sandwich::vcovCL", {
  skip_if_not_installed("sandwich")
  f <- .sandwich_fixture()
  h <- uncounted:::.hat_values_wls(f$X, f$w)
  res <- sqrt(f$w) * f$e
  cl <- rep(1:6, each = 5)

  V_cr1 <- uncounted:::.compute_sandwich_vcov(
    f$X, res, weights = f$w, hat_values = h, vcov_type = "HC1", cluster = cl
  )
  expect_equal(
    unname(V_cr1),
    unname(sandwich::vcovCL(f$fit, cluster = cl, type = "HC1", cadjust = TRUE)),
    tolerance = 1e-8
  )

  ## HC2+ is not defined for the clustered meat and must degrade to CR1.
  expect_warning(
    V_cl5 <- uncounted:::.compute_sandwich_vcov(
      f$X, res, weights = f$w, hat_values = h,
      vcov_type = "HC5", cluster = cl
    ),
    "Falling back to CR1"
  )
  expect_equal(V_cl5, V_cr1, tolerance = 1e-10)

  expect_error(
    uncounted:::.compute_sandwich_vcov(
      f$X, res, weights = f$w, hat_values = h,
      vcov_type = "HC1", cluster = rep(1, f$n)
    ),
    "at least 2 clusters"
  )
})

test_that("NB sandwich helper covers HC1, clustering, and single-cluster error", {
  hessian_nll <- diag(c(2, 3, 4))
  score_full <- matrix(
    c(0.5, -0.3, 0.2,
      -0.2, 0.4, -0.1,
      0.3, -0.2, 0.1,
      -0.1, 0.1, -0.2),
    ncol = 3, byrow = TRUE
  )

  V_hc0 <- uncounted:::.compute_nb_sandwich(
    hessian_nll, score_full, n_obs = nrow(score_full), vcov_type = "HC0"
  )
  V_hc1 <- uncounted:::.compute_nb_sandwich(
    hessian_nll, score_full, n_obs = nrow(score_full), vcov_type = "HC1"
  )
  V_cl <- uncounted:::.compute_nb_sandwich(
    hessian_nll, score_full, n_obs = nrow(score_full),
    vcov_type = "HC1", cluster = rep(1:2, each = 2)
  )

  expect_equal(dim(V_hc0), c(3, 3))
  expect_equal(dim(V_hc1), c(3, 3))
  expect_equal(dim(V_cl), c(3, 3))
  expect_true(sum(diag(V_hc1)) > sum(diag(V_hc0)))
  expect_error(
    uncounted:::.compute_nb_sandwich(
      hessian_nll, score_full, n_obs = nrow(score_full),
      vcov_type = "HC1", cluster = rep(1, nrow(score_full))
    ),
    "at least 2 clusters"
  )
})
