# ---- Constrained estimation tests (CRAN) ----

test_that("constrained Poisson returns alpha in (0,1) and beta > 0", {
  d <- small_data()
  fit <- quick_fit(d, method = "poisson", gamma = 0.005, constrained = TRUE)
  expect_s3_class(fit, "uncounted")
  expect_true(isTRUE(fit$constrained))
  expect_true(all(fit$alpha_values > 0 & fit$alpha_values < 1))
  expect_true(all(fit$beta_values > 0))
})

test_that("constrained NB returns alpha in (0,1) and beta > 0", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005, constrained = TRUE)
  expect_true(all(fit$alpha_values > 0 & fit$alpha_values < 1))
  expect_true(all(fit$beta_values > 0))
})

test_that("constrained popsize gives ordered CI and positive BC", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, constrained = TRUE)
  ps <- popsize(fit, bias_correction = TRUE)
  expect_true(all(ps$lower < ps$estimate))
  expect_true(all(ps$estimate < ps$upper))
  expect_true(all(ps$estimate_bc > 0))
})

test_that("constrained BC does not produce hugely negative values", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, constrained = TRUE)
  ps <- popsize(fit, bias_correction = TRUE)
  expect_true(all(ps$estimate_bc > 0.5 * ps$estimate))
})

test_that("constrained with cov_alpha works", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, constrained = TRUE, cov_alpha = ~ 0 + sex)
  ps <- popsize(fit)
  expect_equal(nrow(ps), 2)
  expect_true(all(fit$alpha_values > 0 & fit$alpha_values < 1))
})

test_that("constrained print has single Constrained line", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, constrained = TRUE)
  out <- capture.output(print(fit))
  constrained_lines <- grep("^Constrained:", out)
  expect_equal(length(constrained_lines), 1)
})

test_that("constrained alpha_values = inv_logit(X %*% alpha_coefs)", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, constrained = TRUE)
  eta <- as.numeric(fit$X_alpha %*% fit$alpha_coefs)
  expected_alpha <- 1 / (1 + exp(-eta))
  expect_equal(fit$alpha_values, expected_alpha, tolerance = 1e-10)
})

# ---- LOO with constrained model ----

test_that("LOO passes constrained argument to refits", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, constrained = TRUE, countries = ~country)
  loo_res <- loo(fit, by = "obs")
  expect_s3_class(loo_res, "uncounted_loo")
  expect_true(sum(loo_res$converged) > 0)
})

test_that("LOO by country with constrained model", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, constrained = TRUE, countries = ~country)
  loo_res <- loo(fit, by = "country")
  expect_equal(loo_res$n_drops, length(unique(d$country)))
})

# ---- Regression: constrained BC includes logit curvature ----

test_that("constrained BC includes logit curvature (second derivative)", {
  d <- small_data()
  fit_c <- quick_fit(d, gamma = 0.005, constrained = TRUE)

  ps_c <- popsize(fit_c, bias_correction = TRUE)

  # BC ratio should be reasonable (not wildly off)
  bc_ratio <- ps_c$estimate_bc / ps_c$estimate
  expect_true(bc_ratio > 0.8 && bc_ratio < 1.0,
              info = paste("Constrained BC ratio:", round(bc_ratio, 4)))

  # Bias correction should be non-trivial
  expect_true(ps_c$estimate > ps_c$estimate_bc,
              info = "Bias correction should reduce the estimate")
})
