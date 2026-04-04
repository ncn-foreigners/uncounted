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

test_that("lrtest warns for non-nested covariate models", {
  skip_on_cran()
  d <- testdata
  # sex-only alpha vs year-only alpha (neither spans the other)
  fit1 <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  fit2 <- quick_fit(d, method = "nb", gamma = 0.005, cov_alpha = ~year)
  expect_warning(lrtest(fit1, fit2), "not be nested")
})

test_that("lrtest does NOT warn when fixed gamma nested in estimated gamma", {
  d <- small_data()
  # Fixed gamma is nested in estimated gamma (same covariates)
  fit1 <- quick_fit(d, gamma = 0.005)
  fit2 <- quick_fit(d, method = "nb", gamma = "estimate")
  expect_no_warning(lrtest(fit1, fit2))
})

test_that("lrtest does NOT warn for genuinely nested models", {
  d <- small_data()
  # Nested: intercept-only is nested in sex-covariate model
  fit1 <- quick_fit(d, gamma = 0.005)                     # p_alpha=1, p_beta=1
  fit2 <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)   # p_alpha=2, p_beta=1
  expect_no_warning(lrtest(fit1, fit2))
})

test_that("lrtest does NOT warn for Poisson vs NB (same covariates)", {
  d <- small_data()
  fit_po <- quick_fit(d, method = "poisson", gamma = 0.005)
  fit_nb <- quick_fit(d, method = "nb", gamma = 0.005)
  expect_no_warning(lrtest(fit_po, fit_nb))
})

test_that("lrtest warns for different methods (not Poisson vs NB)", {
  d <- small_data()
  # Different methods AND different param counts (so df != 0)
  fit1 <- quick_fit(d, method = "poisson", gamma = 0.005)
  fit2 <- quick_fit(d, method = "iols", gamma = 0.005, cov_alpha = ~sex)
  expect_warning(lrtest(fit1, fit2), "different methods")
})

test_that("lrtest warns for ~sex vs ~year (non-nested covariates)", {
  skip_on_cran()
  d <- testdata
  fit1 <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  fit2 <- quick_fit(d, method = "nb", gamma = 0.005, cov_alpha = ~year)
  expect_warning(lrtest(fit1, fit2), "not be nested")
})

test_that("lrtest no warn for ~year nested in ~year + sex", {
  skip_on_cran()
  d <- testdata
  fit1 <- quick_fit(d, gamma = 0.005, cov_alpha = ~year)
  fit2 <- quick_fit(d, gamma = 0.005, cov_alpha = ~year + sex)
  expect_no_warning(lrtest(fit1, fit2))
})

test_that("lrtest warns for different fixed gamma values", {
  d <- small_data()
  fit1 <- quick_fit(d, gamma = 0.005)
  fit2 <- quick_fit(d, gamma = 0.01, cov_alpha = ~sex)
  expect_warning(lrtest(fit1, fit2), "not be nested")
})
