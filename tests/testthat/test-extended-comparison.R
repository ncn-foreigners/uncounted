# ---- compare_loo tests (extended, skip on CRAN) ----

test_that("compare_loo returns correct class", {
  skip_on_cran()
  d <- small_data()
  d$country2 <- d$country
  fit_po <- quick_fit(d, method = "poisson", gamma = 0.005)
  fit_nb <- quick_fit(d, method = "nb", gamma = 0.005)
  loo_po <- loo(fit_po, by = "obs")
  loo_nb <- loo(fit_nb, by = "obs")
  comp <- compare_loo(loo_po, loo_nb, labels = c("Poisson", "NB"))
  expect_s3_class(comp, "uncounted_loo_compare")
  expect_true("table" %in% names(comp))
  expect_equal(nrow(comp$table), nrow(d))
})

test_that("compare_loo prints without error", {
  skip_on_cran()
  d <- small_data()
  fit_po <- quick_fit(d, method = "poisson", gamma = 0.005)
  fit_nb <- quick_fit(d, method = "nb", gamma = 0.005)
  loo_po <- loo(fit_po, by = "obs")
  loo_nb <- loo(fit_nb, by = "obs")
  comp <- compare_loo(loo_po, loo_nb, labels = c("Poisson", "NB"))
  expect_output(print(comp), "LOO comparison")
})

test_that("compare_loo scatter plot runs without error", {
  skip_on_cran()
  d <- small_data()
  fit_po <- quick_fit(d, method = "poisson", gamma = 0.005)
  fit_nb <- quick_fit(d, method = "nb", gamma = 0.005)
  loo_po <- loo(fit_po, by = "obs")
  loo_nb <- loo(fit_nb, by = "obs")
  comp <- compare_loo(loo_po, loo_nb)
  expect_no_error(plot(comp, type = "scatter"))
})

test_that("compare_loo bar plot runs without error", {
  skip_on_cran()
  d <- small_data()
  fit_po <- quick_fit(d, method = "poisson", gamma = 0.005)
  fit_nb <- quick_fit(d, method = "nb", gamma = 0.005)
  loo_po <- loo(fit_po, by = "obs")
  loo_nb <- loo(fit_nb, by = "obs")
  comp <- compare_loo(loo_po, loo_nb)
  expect_no_error(plot(comp, type = "bar"))
})

test_that("compare_loo errors with mismatched by", {
  skip_on_cran()
  fit_po <- quick_fit(testdata, method = "poisson", gamma = 0.005,
                      countries = ~country)
  fit_nb <- quick_fit(testdata, method = "nb", gamma = 0.005,
                      countries = ~country)
  loo_po <- loo(fit_po, by = "obs")
  loo_nb <- loo(fit_nb, by = "country")
  expect_error(compare_loo(loo_po, loo_nb), "same 'by'")
})

test_that("compare_loo table is sorted by max influence", {
  skip_on_cran()
  d <- small_data()
  fit_po <- quick_fit(d, method = "poisson", gamma = 0.005)
  fit_nb <- quick_fit(d, method = "nb", gamma = 0.005)
  loo_po <- loo(fit_po, by = "obs")
  loo_nb <- loo(fit_nb, by = "obs")
  comp <- compare_loo(loo_po, loo_nb)
  max_abs <- comp$table$max_abs_pct
  expect_true(all(diff(max_abs) <= 0))
})

# ---- LOO summary and plot coverage ----

test_that("summary.uncounted_loo runs without error", {
  skip_on_cran()
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, countries = ~country)
  loo_res <- loo(fit, by = "obs")
  expect_output(summary(loo_res))
})

test_that("plot(loo_res, type = 'coef') runs without error", {
  skip_on_cran()
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, countries = ~country)
  loo_res <- loo(fit, by = "obs")
  expect_no_error(plot(loo_res, type = "coef"))
})

test_that("compare_loo with data and label_vars works", {
  skip_on_cran()
  d <- small_data()
  fit_po <- quick_fit(d, method = "poisson", gamma = 0.005, countries = ~country)
  fit_nb <- quick_fit(d, method = "nb", gamma = 0.005, vcov = "HC0",
                      countries = ~country)
  loo_po <- loo(fit_po, by = "country")
  loo_nb <- loo(fit_nb, by = "country")
  comp <- compare_loo(loo_po, loo_nb, data = d, label_vars = ~country)
  expect_s3_class(comp, "uncounted_loo_compare")
})

test_that("R^2_D uses null model without covariates", {
  skip_on_cran()
  testdata$grp <- factor(rep(1:2, length.out = nrow(testdata)))
  fit_cov  <- quick_fit(testdata, gamma = 0.005, cov_alpha = ~grp)
  fit_null <- quick_fit(testdata, gamma = 0.005)
  comp <- compare_models(WithCov = fit_cov, NoCov = fit_null)
  expect_true(comp$table$R2_D[comp$table$Model == "WithCov"] >= -0.01)
  expect_true(abs(comp$table$R2_D[comp$table$Model == "NoCov"]) < 0.01)
})
