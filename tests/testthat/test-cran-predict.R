# ---- predict.uncounted tests ----

test_that("predict returns fitted values with no newdata", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  pred <- predict(fit)
  expect_equal(pred, fitted(fit))
})

test_that("predict with type='link' returns log scale", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  pred_link <- predict(fit, type = "link")
  pred_resp <- predict(fit, type = "response")
  expect_equal(exp(pred_link), pred_resp, tolerance = 1e-10)
})

test_that("predict with newdata reproduces fitted values", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  pred_orig <- predict(fit, newdata = d)
  expect_equal(pred_orig, fitted(fit), tolerance = 1e-10)
})

test_that("predict with newdata subset works", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  pred_sub <- predict(fit, newdata = d[1:5, ])
  expect_length(pred_sub, 5)
  expect_equal(pred_sub, fitted(fit)[1:5], tolerance = 1e-10)
})

test_that("predict for all methods", {
  d <- small_data()
  for (method in c("ols", "poisson", "nb", "iols")) {
    suppressWarnings({
      fit <- quick_fit(d, method = method, gamma = 0.005)
    })
    pred <- predict(fit)
    expect_length(pred, nrow(d))
    expect_true(all(is.finite(pred)), info = method)
  }
})

test_that("predict with covariates reproduces fitted", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  pred <- predict(fit, newdata = d)
  expect_equal(pred, fitted(fit), tolerance = 1e-10)
})

# ---- Regression: factor encoding safety ----

test_that("predict with reordered factor levels gives same result", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  # Reverse factor levels
  d2 <- d
  d2$sex <- factor(d2$sex, levels = rev(levels(factor(d$sex))))
  pred_orig <- predict(fit, newdata = d)
  pred_reord <- predict(fit, newdata = d2)
  expect_equal(pred_orig, pred_reord, tolerance = 1e-10)
})

test_that("predict with single-level factor subset works", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  # Subset to only M
  d_m <- d[d$sex == "M", ]
  pred <- predict(fit, newdata = d_m)
  expect_length(pred, nrow(d_m))
  # Should match fitted for those rows
  idx_m <- which(d$sex == "M")
  expect_equal(pred, fitted(fit)[idx_m], tolerance = 1e-10)
})

test_that("predict caps cov_gamma linear predictor like the fitting path", {
  d <- small_data()
  d$z <- rep(c(0, 1), length.out = nrow(d))

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

  fit$gamma_coefs[] <- 100
  newdata <- d[1:5, , drop = FALSE]
  newdata$z <- 1

  pred_link <- predict(fit, newdata = newdata, type = "link")
  X_gamma_new <- uncounted:::.predict_design(~ 0 + z, newdata, fit$X_gamma,
                                             "gamma", fit$data)
  gamma_vals <- exp(pmin(as.numeric(X_gamma_new %*% fit$gamma_coefs), 10))
  rate_vals <- uncounted:::.rate_from_gamma(newdata$n / newdata$N, gamma_vals)
  X_alpha_new <- uncounted:::.predict_design(NULL, newdata, fit$X_alpha,
                                             "alpha", fit$data)
  X_beta_new <- uncounted:::.predict_design(NULL, newdata, fit$X_beta,
                                            "beta", fit$data)
  manual_link <- uncounted:::.compute_log_mu(
    alpha_values = as.numeric(X_alpha_new %*% fit$alpha_coefs),
    log_N = log(newdata$N),
    beta_values = as.numeric(X_beta_new %*% fit$beta_coefs),
    rate_values = rate_vals,
    link_rho = fit$link_rho
  )

  expect_true(all(is.finite(pred_link)))
  expect_equal(pred_link, as.numeric(manual_link), tolerance = 1e-10)
})

test_that("predict() reproduces fitted() under a dependence offset", {
  ## profile_dependence() refits inside .with_dependence_offset(), which sets
  ## process-global state read by .compute_log_mu(). That state is back to 0
  ## by the time predict() runs, so predict(fit, newdata = <training data>)
  ## used to differ from fitted(fit) by exactly exp(-delta).
  d <- small_data()
  delta <- 0.3

  fit_plain <- quick_fit(d, method = "poisson", gamma = 0.005)
  expect_equal(fit_plain$dependence_offset, 0)

  fit_off <- uncounted:::.with_dependence_offset(
    delta, quick_fit(d, method = "poisson", gamma = 0.005)
  )
  ## The offset is recorded on the fit, and global state is restored.
  expect_equal(fit_off$dependence_offset, delta)
  expect_equal(uncounted:::.current_dependence_offset(), 0)

  ## The offset enters the likelihood, so the fit adapts to it: coefficients
  ## and fitted values genuinely differ from the plain fit. That is what makes
  ## the predict() mismatch below a real inconsistency rather than a constant.
  expect_false(isTRUE(all.equal(unname(coef(fit_off)),
                                unname(coef(fit_plain)))))

  expect_equal(unname(predict(fit_off, newdata = d)), unname(fitted(fit_off)),
               tolerance = 1e-10)
  expect_equal(unname(predict(fit_off, newdata = d, type = "link")),
               unname(fit_off$log_mu), tolerance = 1e-10)

  ## A fit object from before this field existed must still predict.
  fit_legacy <- fit_off
  fit_legacy$dependence_offset <- NULL
  expect_silent(p <- predict(fit_legacy, newdata = d))
  expect_length(p, nrow(d))
})

test_that("predict(newdata) evaluates expressions in reference_pop and auxiliary", {
  d <- testdata
  d$N2 <- d$N / 2
  d$n2 <- d$n / 2
  fit <- estimate_hidden_pop(d, ~m, ~I(n2 * 2), ~I(N2 * 2), gamma = 0.005)
  expect_equal(unname(predict(fit, newdata = d)), unname(fitted(fit)))

  fit_log <- estimate_hidden_pop(d, ~m, ~n, ~I(2 * N2), gamma = 0.005)
  expect_equal(unname(predict(fit_log, newdata = d, type = "link")), unname(fit_log$log_mu))
})
