# ---- link_rho and count-estimator tests (CRAN) ----

test_that("bounded link_rho choices keep rho in (0, 1) and preserve xi identity", {
  d <- small_data()

  for (link_name in c("cloglog", "logit", "probit")) {
    fit <- quick_fit(d, method = "poisson", gamma = 0.005, link_rho = link_name)

    expect_true(all(fit$rho_values > 0))
    expect_true(all(fit$rho_values < 1))
    expect_equal(
      fit$fitted.values / fit$rho_values,
      d$N^fit$alpha_values,
      tolerance = 1e-8
    )
  }
})

test_that("predict works with non-power detection links", {
  d <- small_data()
  fit <- quick_fit(d, method = "poisson", gamma = 0.005, link_rho = "probit")

  expect_equal(as.numeric(predict(fit)), as.numeric(fit$fitted.values))
  expect_true(all(is.finite(predict(fit, type = "link"))))
})

test_that("NLS supports the probit detection link", {
  d <- small_data()
  fit <- quick_fit(d, method = "nls", gamma = 0.005, link_rho = "probit")

  expect_identical(fit$link_rho, "probit")
  expect_true(all(fit$rho_values > 0))
  expect_true(all(fit$rho_values < 1))
  expect_true(all(is.finite(fit$fitted.values)))
})

test_that("legacy logistic alias is accepted and normalized to logit", {
  d <- small_data()
  fit <- quick_fit(d, method = "poisson", gamma = 0.005, link_rho = "logistic")

  expect_identical(fit$link_rho, "logit")
  expect_true(all(fit$rho_values > 0))
  expect_true(all(fit$rho_values < 1))
})

test_that("Poisson GMM fit returns a usable uncounted object", {
  d <- small_data()
  fit <- quick_fit(d, method = "poisson", gamma = 0.005,
                   estimator = "gmm", link_rho = "probit")

  expect_s3_class(fit, "uncounted")
  expect_identical(fit$estimator, "gmm")
  expect_equal(fit$link_rho, "probit")
  expect_true(all(is.finite(coef(fit))))
  expect_true(all(is.finite(fit$fitted.values)))
  expect_true(is.na(as.numeric(logLik(fit))))
  expect_true(is.na(AIC(fit)))
  expect_true(is.na(BIC(fit)))
})

test_that("NB EL fit stores theta-aware backend covariance", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005,
                   estimator = "el", link_rho = "probit")

  expect_s3_class(fit, "uncounted")
  expect_identical(fit$estimator, "el")
  expect_identical(fit$link_rho, "probit")
  expect_true(is.finite(fit$theta))
  expect_true(is.finite(fit$theta_se))
  expect_equal(dim(vcov_nb(fit)), c(length(coef(fit)) + 1, length(coef(fit)) + 1))
  expect_true(is.na(as.numeric(logLik(fit))))
})

test_that("vcovHC and vcovCL work for GMM fits", {
  d <- small_data()
  fit <- quick_fit(d, method = "poisson", gamma = 0.005, estimator = "gmm")

  V_hc0 <- sandwich::vcovHC(fit, type = "HC0")
  V_cl <- sandwich::vcovCL(fit, cluster = d$country, type = "HC1")

  expect_equal(dim(V_hc0), rep(length(coef(fit)), 2))
  expect_equal(dim(V_cl), rep(length(coef(fit)), 2))
  expect_warning(
    V_hc3 <- sandwich::vcovHC(fit, type = "HC3"),
    "Falling back to HC1"
  )
  expect_equal(V_hc3, sandwich::vcovHC(fit, type = "HC1"))
})

test_that("vcovFWB works with a Poisson GMM fit", {
  skip_if_not_installed("fwb")
  set.seed(123)
  d <- small_data()
  fit <- quick_fit(d, method = "poisson", gamma = 0.005, estimator = "gmm")

  V_fwb <- fwb::vcovFWB(fit, R = 20)
  expect_equal(dim(V_fwb), rep(length(coef(fit)), 2))
  expect_true(all(diag(V_fwb) > 0))
})

test_that("Poisson GMM vcov_model uses Fisher-style count-model weights", {
  d <- small_data()
  fit <- quick_fit(d, method = "poisson", gamma = 0.005, estimator = "gmm")

  V_expected <- uncounted:::.compute_model_vcov(
    fit$model_matrix_full,
    fit$bread_weights
  )

  expect_equal(as.numeric(fit$vcov_model), as.numeric(V_expected), tolerance = 1e-10)
})

test_that("NB EL vcov_model uses Fisher-style mean-model weights", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005, estimator = "el")

  V_expected <- uncounted:::.compute_model_vcov(
    fit$model_matrix_full,
    fit$bread_weights
  )

  expect_equal(as.numeric(fit$vcov_model), as.numeric(V_expected), tolerance = 1e-10)
})

test_that("plug-in xi does not depend on the covariance choice for GMM fits", {
  d <- small_data()
  fit_hc0 <- quick_fit(d, method = "poisson", gamma = 0.005,
                       estimator = "gmm", vcov = "HC0")
  fit_hc1 <- quick_fit(d, method = "poisson", gamma = 0.005,
                       estimator = "gmm", vcov = "HC1")

  expect_equal(popsize(fit_hc0)$estimate, popsize(fit_hc1)$estimate, tolerance = 1e-10)
})

test_that("non-MLE likelihood-only helpers error clearly", {
  d <- small_data()
  fit_gmm <- quick_fit(d, method = "poisson", gamma = 0.005, estimator = "gmm")
  fit_el <- quick_fit(d, method = "nb", gamma = 0.005, estimator = "el")

  expect_error(lrtest(fit_gmm, fit_el), "estimator = 'mle'")
  expect_error(profile_gamma(fit_gmm, plot = FALSE), "estimator = 'mle'")
  expect_error(profile_alpha(fit_gmm, plot = FALSE), "estimator = 'mle'")
})

test_that("compare_models warns for mixed estimators and keeps likelihood columns NA", {
  d <- small_data()
  fit_mle <- quick_fit(d, method = "poisson", gamma = 0.005)
  fit_gmm <- quick_fit(d, method = "poisson", gamma = 0.005, estimator = "gmm")

  expect_warning(comp <- compare_models(MLE = fit_mle, GMM = fit_gmm),
                 "Mixed-estimator")
  gmm_row <- comp$table[comp$table$Estimator == "GMM", ]
  expect_true(is.na(gmm_row$logLik))
  expect_true(is.na(gmm_row$AIC))
  expect_true(is.na(gmm_row$BIC))
})

test_that("update and loo preserve estimator metadata", {
  d <- small_data()
  fit <- quick_fit(d, method = "poisson", gamma = 0.005,
                   estimator = "gmm", link_rho = "probit",
                   countries = ~country)

  fit_upd <- update(fit, vcov = "HC0")
  loo_res <- loo(fit, by = "obs")

  expect_identical(fit_upd$estimator, "gmm")
  expect_identical(fit_upd$link_rho, "probit")
  expect_s3_class(loo_res, "uncounted_loo")
  expect_true(sum(loo_res$converged) > 0)
})

test_that("bootstrap_popsize works for a small Poisson GMM fit", {
  skip_if_not_installed("fwb")
  set.seed(321)
  d <- small_data()
  fit <- quick_fit(d, method = "poisson", gamma = 0.005,
                   estimator = "gmm", link_rho = "probit")

  boot <- bootstrap_popsize(fit, R = 9, seed = 123, verbose = FALSE)
  expect_s3_class(boot, "uncounted_boot")
  expect_true(boot$n_converged > 0)
})
