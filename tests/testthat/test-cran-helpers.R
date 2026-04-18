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

test_that("model-based and sandwich covariance helpers return finite matrices", {
  X <- cbind("(Intercept)" = 1, x = seq_len(6))
  residuals <- c(-1, 1, -0.5, 0.5, -0.25, 0.25)
  weights <- seq(1, 6)
  hat_values <- uncounted:::.hat_values_wls(X, weights)

  V_model <- uncounted:::.compute_model_vcov(X, weights, sigma2 = 2)
  V_hc4 <- uncounted:::.compute_sandwich_vcov(
    X, residuals, weights = weights, hat_values = hat_values, vcov_type = "HC4"
  )
  expect_equal(
    V_model,
    2 * solve(crossprod(X * sqrt(weights))),
    tolerance = 1e-10
  )
  expect_equal(dim(V_hc4), c(ncol(X), ncol(X)))
  expect_true(all(is.finite(V_hc4)))

  expect_warning(
    V_cl <- uncounted:::.compute_sandwich_vcov(
      X, residuals, weights = weights, hat_values = hat_values,
      vcov_type = "HC5", cluster = rep(1:3, each = 2)
    ),
    "Falling back to CR1"
  )
  expect_equal(dim(V_cl), c(ncol(X), ncol(X)))
  expect_true(all(is.finite(V_cl)))
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
