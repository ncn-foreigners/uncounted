## Tests for covariate-varying gamma (cov_gamma)
## Skip on CRAN — simulation-based and compute-intensive

skip_on_cran()

# ---- Shared simulation data ----
make_cov_gamma_data <- function(n = 300, seed = 42) {
  set.seed(seed)
  sex <- rep(c(0, 1), each = n / 2)
  N <- sample(1000:50000, n, replace = TRUE)
  n_aux <- rpois(n, lambda = N * 0.05)

  alpha_true <- 0.7
  beta_true <- 0.5
  gamma_int <- log(0.005)   # intercept on log scale
  gamma_sex <- 0.8          # sex effect on log scale
  gamma_true <- exp(gamma_int + gamma_sex * sex)

  mu <- N^alpha_true * (gamma_true + n_aux / N)^beta_true
  m <- rpois(n, lambda = mu)

  list(
    data = data.frame(N = N, n = n_aux, m = m, sex = factor(sex)),
    true = list(alpha = alpha_true, beta = beta_true,
                gamma_int = gamma_int, gamma_sex = gamma_sex,
                gamma_true = gamma_true)
  )
}

# ---- Validation errors ----

test_that("cov_gamma errors with fixed gamma", {
  d <- small_data()
  expect_error(
    estimate_hidden_pop(d, ~ m, ~ n, ~ N,
                        method = "poisson", gamma = 0.005, cov_gamma = ~ sex),
    "cov_gamma.*requires.*estimate"
  )
})

test_that("cov_gamma errors with NULL gamma", {
  d <- small_data()
  expect_error(
    estimate_hidden_pop(d, ~ m, ~ n, ~ N,
                        method = "poisson", gamma = NULL, cov_gamma = ~ sex),
    "cov_gamma.*requires.*estimate"
  )
})

test_that("cov_gamma errors with OLS method", {
  d <- small_data()
  expect_error(
    estimate_hidden_pop(d, ~ m, ~ n, ~ N,
                        method = "ols", gamma = "estimate", cov_gamma = ~ sex),
    "only supported for.*poisson.*nb"
  )
})

test_that("cov_gamma errors with NLS method", {
  d <- small_data()
  expect_error(
    estimate_hidden_pop(d, ~ m, ~ n, ~ N,
                        method = "nls", gamma = "estimate", cov_gamma = ~ sex),
    "only supported for.*poisson.*nb"
  )
})

# ---- Poisson smoke test ----

test_that("cov_gamma works with Poisson method", {
  sim <- make_cov_gamma_data()
  fit <- estimate_hidden_pop(sim$data, ~ m, ~ n, ~ N,
                             method = "poisson", gamma = "estimate",
                             cov_gamma = ~ sex)

  expect_s3_class(fit, "uncounted")
  expect_true(!is.null(fit$gamma_coefs))
  expect_equal(length(fit$gamma_coefs), 2)
  expect_true(!is.null(fit$gamma_values))
  expect_equal(length(fit$gamma_values), nrow(sim$data))
  expect_true(all(fit$gamma_values > 0))
  expect_equal(fit$p_gamma, 2L)
  expect_true(fit$gamma_estimated)
})

# ---- NB smoke test ----

test_that("cov_gamma works with NB method", {
  sim <- make_cov_gamma_data()
  fit <- estimate_hidden_pop(sim$data, ~ m, ~ n, ~ N,
                             method = "nb", gamma = "estimate",
                             cov_gamma = ~ sex)

  expect_s3_class(fit, "uncounted")
  expect_true(!is.null(fit$gamma_coefs))
  expect_equal(length(fit$gamma_coefs), 2)
  expect_true(!is.null(fit$theta))
  expect_true(fit$gamma_estimated)
})

# ---- Backward compatibility ----

test_that("cov_gamma = NULL preserves existing Poisson behavior", {
  d <- positive_data()
  fit1 <- estimate_hidden_pop(d, ~ m, ~ n, ~ N,
                              method = "poisson", gamma = "estimate")
  fit2 <- estimate_hidden_pop(d, ~ m, ~ n, ~ N,
                              method = "poisson", gamma = "estimate",
                              cov_gamma = NULL)

  expect_equal(coef(fit1), coef(fit2))
  expect_equal(fit1$gamma, fit2$gamma, tolerance = 1e-8)
  expect_equal(fit1$loglik, fit2$loglik, tolerance = 1e-8)
})

# ---- Coefficient recovery simulation ----

test_that("cov_gamma recovers varying gamma from Poisson simulation", {
  sim <- make_cov_gamma_data(n = 400, seed = 123)

  fit <- estimate_hidden_pop(sim$data, ~ m, ~ n, ~ N,
                             method = "poisson", gamma = "estimate",
                             cov_gamma = ~ sex)

  # Sex coefficient should be positive and roughly correct
  gamma_sex_hat <- fit$gamma_coefs[2]
  expect_true(gamma_sex_hat > 0,
              info = paste("gamma sex coef =", round(gamma_sex_hat, 3)))
  expect_true(abs(gamma_sex_hat - sim$true$gamma_sex) < 2,
              info = paste("gamma sex coef =", round(gamma_sex_hat, 3),
                           "true =", sim$true$gamma_sex))

  # Gamma values should differ between sexes
  gv_female <- unique(round(fit$gamma_values[sim$data$sex == 0], 8))
  gv_male <- unique(round(fit$gamma_values[sim$data$sex == 1], 8))
  expect_equal(length(gv_female), 1)
  expect_equal(length(gv_male), 1)
  expect_false(gv_female == gv_male)
})

test_that("cov_gamma recovers varying gamma from NB simulation", {
  set.seed(99)
  n <- 400
  sex <- rep(c(0, 1), each = n / 2)
  N <- sample(1000:50000, n, replace = TRUE)
  n_aux <- rpois(n, lambda = N * 0.05)
  gamma_true <- exp(log(0.005) + 0.8 * sex)
  mu <- N^0.7 * (gamma_true + n_aux / N)^0.5
  m <- MASS::rnegbin(n, mu = mu, theta = 5)
  d <- data.frame(N = N, n = n_aux, m = m, sex = factor(sex))

  fit <- estimate_hidden_pop(d, ~ m, ~ n, ~ N,
                             method = "nb", gamma = "estimate",
                             cov_gamma = ~ sex)

  expect_true(fit$gamma_coefs[2] > 0)
  expect_true(!is.null(fit$theta))
})

# ---- Vcov dimensions ----

test_that("vcov has correct dimensions with cov_gamma", {
  sim <- make_cov_gamma_data()
  fit <- estimate_hidden_pop(sim$data, ~ m, ~ n, ~ N,
                             method = "poisson", gamma = "estimate",
                             cov_gamma = ~ sex)

  # vcov includes alpha + beta + gamma coefs
  p_total <- fit$p_alpha + fit$p_beta + fit$p_gamma
  expect_equal(nrow(fit$vcov), p_total)
  expect_equal(ncol(fit$vcov), p_total)
  # vcov_full includes same or more
  expect_true(nrow(fit$vcov_full) >= p_total)
})

# ---- print and summary ----

test_that("print works with cov_gamma", {
  sim <- make_cov_gamma_data()
  fit <- estimate_hidden_pop(sim$data, ~ m, ~ n, ~ N,
                             method = "poisson", gamma = "estimate",
                             cov_gamma = ~ sex)
  expect_output(print(fit), "covariate-varying")
})

test_that("summary works with cov_gamma", {
  sim <- make_cov_gamma_data()
  fit <- estimate_hidden_pop(sim$data, ~ m, ~ n, ~ N,
                             method = "poisson", gamma = "estimate",
                             cov_gamma = ~ sex)
  expect_output(summary(fit), "Gamma.*response scale")
})

# ---- tidy output ----

test_that("tidy returns gamma coefficient rows with cov_gamma", {
  skip_if_not_installed("broom")
  sim <- make_cov_gamma_data()
  fit <- estimate_hidden_pop(sim$data, ~ m, ~ n, ~ N,
                             method = "poisson", gamma = "estimate",
                             cov_gamma = ~ sex)
  td <- broom::tidy(fit)
  gamma_rows <- td[grep("gamma:", td$term), ]
  expect_true(nrow(gamma_rows) >= 2)
  expect_true(all(is.finite(gamma_rows$std.error)))
})

# ---- predict with newdata ----

test_that("predict uses per-observation gamma with cov_gamma", {
  sim <- make_cov_gamma_data()
  fit <- estimate_hidden_pop(sim$data, ~ m, ~ n, ~ N,
                             method = "poisson", gamma = "estimate",
                             cov_gamma = ~ sex)
  pred <- predict(fit, newdata = sim$data[1:10, ])
  expect_equal(length(pred), 10)
  expect_true(all(is.finite(pred)))
  expect_true(all(pred > 0))
})

# ---- AIC comparison ----

test_that("AIC favors cov_gamma when DGP has varying gamma", {
  sim <- make_cov_gamma_data(n = 400, seed = 77)

  fit_no <- estimate_hidden_pop(sim$data, ~ m, ~ n, ~ N,
                                method = "poisson", gamma = "estimate")
  fit_cov <- estimate_hidden_pop(sim$data, ~ m, ~ n, ~ N,
                                 method = "poisson", gamma = "estimate",
                                 cov_gamma = ~ sex)
  # Model with cov_gamma should generally have better (lower) AIC
  expect_true(AIC(fit_cov) < AIC(fit_no),
              info = paste("AIC cov_gamma:", round(AIC(fit_cov), 1),
                           "AIC no cov:", round(AIC(fit_no), 1)))
})

# ---- Constrained + cov_gamma ----

test_that("constrained estimation works with cov_gamma", {
  sim <- make_cov_gamma_data()
  fit <- estimate_hidden_pop(sim$data, ~ m, ~ n, ~ N,
                             method = "poisson", gamma = "estimate",
                             cov_gamma = ~ sex, constrained = TRUE)
  expect_s3_class(fit, "uncounted")
  expect_true(all(fit$alpha_values > 0 & fit$alpha_values < 1))
  expect_true(all(fit$beta_values > 0))
  expect_true(all(fit$gamma_values > 0))
})

# ---- P1 regression: one-column cov_gamma (~ 0 + z) ----

test_that("one-column cov_gamma formula works (not misclassified as scalar)", {
  sim <- make_cov_gamma_data()
  # Continuous covariate, no intercept → p_gamma == 1 but has_cov_gamma == TRUE
  sim$data$z <- rnorm(nrow(sim$data))
  fit <- estimate_hidden_pop(sim$data, ~ m, ~ n, ~ N,
                             method = "poisson", gamma = "estimate",
                             cov_gamma = ~ 0 + z)

  expect_true(isTRUE(fit$has_cov_gamma))
  expect_equal(fit$p_gamma, 1L)
  # Gamma coefs should be in all_coefs
  expect_true("gamma:z" %in% names(coef(fit)))
  # vcov should include the gamma coefficient
  expect_true("gamma:z" %in% colnames(vcov(fit)))
  # predict should use per-observation gamma (not scalar)
  pred <- predict(fit, newdata = sim$data[1:5, ])
  expect_equal(length(pred), 5)
  expect_true(all(is.finite(pred)))
})

# ---- P2a regression: tidy() scalar gamma SE uses delta method ----

test_that("tidy() reports correct response-scale SE for scalar gamma", {
  skip_if_not_installed("broom")
  d <- positive_data()
  fit <- estimate_hidden_pop(d, ~ m, ~ n, ~ N,
                             method = "poisson", gamma = "estimate")

  td <- broom::tidy(fit)
  gamma_row <- td[td$term == "gamma", ]
  expect_equal(nrow(gamma_row), 1)

  # SE should be reasonable relative to estimate (not orders of magnitude off)
  gamma_val <- gamma_row$estimate
  gamma_se <- gamma_row$std.error
  expect_true(is.finite(gamma_se))
  # The SE should be on the same order as the estimate (not 100x larger)
  expect_true(gamma_se < 10 * gamma_val,
              info = paste("gamma =", round(gamma_val, 6),
                           "SE =", round(gamma_se, 6)))
})

# ---- P2b regression: lrtest() warns for non-nested gamma submodels ----

test_that("lrtest warns when gamma submodels are not nested", {
  sim <- make_cov_gamma_data()
  sim$data$year <- factor(rep(1:3, length.out = nrow(sim$data)))

  fit1 <- estimate_hidden_pop(sim$data, ~ m, ~ n, ~ N,
                              method = "poisson", gamma = "estimate",
                              cov_gamma = ~ sex)
  fit2 <- estimate_hidden_pop(sim$data, ~ m, ~ n, ~ N,
                              method = "poisson", gamma = "estimate",
                              cov_gamma = ~ year)
  # These have non-nested gamma designs — lrtest should warn
  expect_warning(lrtest(fit1, fit2), "not.*nested")
})

test_that("lrtest warns for fixed-gamma vs cov_gamma = ~ 0 + z (no intercept)", {
  sim <- make_cov_gamma_data()
  sim$data$z <- rnorm(nrow(sim$data))

  fit_fixed <- estimate_hidden_pop(sim$data, ~ m, ~ n, ~ N,
                                   method = "poisson", gamma = 0.005)
  fit_cov <- estimate_hidden_pop(sim$data, ~ m, ~ n, ~ N,
                                 method = "poisson", gamma = "estimate",
                                 cov_gamma = ~ 0 + z)
  # ~ 0 + z has no constant column, so fixed gamma is NOT nested
  expect_warning(lrtest(fit_fixed, fit_cov), "not.*nested")
})

# ---- profile_gamma() on cov_gamma models ----

test_that("profile_gamma errors on cov_gamma models", {
  sim <- make_cov_gamma_data()
  fit <- estimate_hidden_pop(sim$data, ~ m, ~ n, ~ N,
                             method = "poisson", gamma = "estimate",
                             cov_gamma = ~ sex)
  expect_error(profile_gamma(fit, plot = FALSE), "cov_gamma")
})
