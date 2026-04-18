# ---- LOO tests (extended, skip on CRAN) ----

skip_on_cran()

test_that("LOO by obs on full data (200 refits)", {
  skip_on_cran()
  fit <- quick_fit(testdata, gamma = 0.005, countries = ~country)
  loo_res <- loo(fit, by = "obs")
  expect_s3_class(loo_res, "uncounted_loo")
  expect_equal(loo_res$n_drops, nrow(testdata))
  expect_true(sum(loo_res$converged) > nrow(testdata) * 0.9)
})

test_that("LOO by country on full data (20 refits)", {
  skip_on_cran()
  fit <- quick_fit(testdata, gamma = 0.005, countries = ~country)
  loo_res <- loo(fit, by = "country")
  expect_equal(loo_res$n_drops, length(unique(testdata$country)))
  expect_true(sum(loo_res$converged) > 15)
})

test_that("LOO with constrained model on full data", {
  skip_on_cran()
  fit <- quick_fit(testdata, gamma = 0.005, constrained = TRUE,
                   countries = ~country)
  loo_res <- loo(fit, by = "obs")
  expect_s3_class(loo_res, "uncounted_loo")
  expect_true(sum(loo_res$converged) > 0)
})

test_that("LOO with NB model on full data", {
  skip_on_cran()
  fit <- quick_fit(testdata, method = "nb", gamma = 0.005,
                   countries = ~country)
  loo_res <- loo(fit, by = "country")
  expect_equal(loo_res$n_drops, length(unique(testdata$country)))
})

test_that("LOO by country with constrained model on full data", {
  skip_on_cran()
  fit <- quick_fit(testdata, gamma = 0.005, constrained = TRUE,
                   countries = ~country)
  loo_res <- loo(fit, by = "country")
  expect_equal(loo_res$n_drops, length(unique(testdata$country)))
})
