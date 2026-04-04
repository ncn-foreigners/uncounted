# ---- Model comparison tests (CRAN) ----

test_that("compare_models works with 2 models", {
  d <- small_data()
  fit_po <- quick_fit(d, method = "poisson", gamma = 0.005)
  fit_nb <- quick_fit(d, method = "nb", gamma = 0.005)
  comp <- compare_models(fit_po, fit_nb)
  expect_s3_class(comp, "uncounted_comparison")
  expect_equal(nrow(comp$table), 2)
  expect_true(all(c("AIC", "BIC", "logLik") %in% names(comp$table)))
})

test_that("compare_models includes pseudo R^2 columns", {
  d <- small_data()
  fit_po <- quick_fit(d, method = "poisson", gamma = 0.005)
  fit_nb <- quick_fit(d, method = "nb", gamma = 0.005)
  comp <- compare_models(fit_po, fit_nb)
  expect_true(all(c("R2_cor", "R2_D", "R2_CW") %in% names(comp$table)))
  expect_true(all(comp$table$R2_cor >= 0 & comp$table$R2_cor <= 1))
})

test_that("pseudo R^2_D is NA for OLS models", {
  d <- positive_data(small_data())
  fit_ols <- quick_fit(d, method = "ols", gamma = NULL)
  fit_po  <- quick_fit(d, method = "poisson", gamma = 0.005)
  comp <- suppressWarnings(compare_models(OLS = fit_ols, Poisson = fit_po))
  ols_row <- comp$table[comp$table$Model == "OLS", ]
  expect_true(is.na(ols_row$R2_D))
  expect_true(is.na(ols_row$R2_CW))
})

test_that("compare_models prints without error", {
  d <- small_data()
  fit_po <- quick_fit(d, method = "poisson", gamma = 0.005)
  fit_nb <- quick_fit(d, method = "nb", gamma = 0.005)
  expect_output(print(compare_models(fit_po, fit_nb)), "Model comparison")
})

test_that("compare_models warns with mixed types", {
  d <- positive_data(small_data())
  fit_ols <- quick_fit(d, method = "ols", gamma = NULL)
  fit_po  <- quick_fit(d, method = "poisson", gamma = 0.005)
  expect_warning(compare_models(fit_ols, fit_po), "pseudo-loglik")
})

test_that("lrtest works for Poisson vs NB", {
  d <- small_data()
  fit_po <- quick_fit(d, method = "poisson", gamma = 0.005)
  fit_nb <- quick_fit(d, method = "nb", gamma = 0.005)
  lr <- lrtest(fit_po, fit_nb)
  expect_s3_class(lr, "uncounted_lrtest")
  expect_true(lr$statistic >= 0)
  expect_true(lr$p_value >= 0 && lr$p_value <= 1)
  expect_true(lr$boundary)
})

test_that("lrtest prints without error", {
  d <- small_data()
  fit_po <- quick_fit(d, method = "poisson", gamma = 0.005)
  fit_nb <- quick_fit(d, method = "nb", gamma = 0.005)
  expect_output(print(lrtest(fit_po, fit_nb)), "Likelihood ratio test")
})

test_that("lrtest errors with OLS models", {
  d <- positive_data(small_data())
  fit1 <- quick_fit(d, method = "ols", gamma = NULL)
  fit2 <- quick_fit(d, method = "ols", gamma = NULL)
  expect_error(lrtest(fit1, fit2), "log-likelihood")
})

# ---- Regression: lrtest warns for non-nested models ----

test_that("lrtest warns for non-nested models", {
  d <- small_data()
  # Non-nested: different covariates with different parameter counts
  # so the df != 0 check doesn't trigger, but nesting check does
  fit1 <- quick_fit(d, gamma = 0.005)  # 2 params (alpha, beta)
  fit2 <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex, cov_beta = ~sex)  # 4 params
  # fit1's coef names (alpha, beta) are NOT a subset of fit2's (alpha:sexF, alpha:sexM, ...)
  expect_warning(lrtest(fit1, fit2), "not be nested")
})
