# ---- Diagnostics on full data (extended, skip on CRAN) ----

test_that("dfbeta returns correct dimensions on full data", {
  skip_on_cran()
  fit <- quick_fit(testdata, gamma = 0.005)
  db <- dfbeta(fit)
  expect_true(is.matrix(db))
  expect_equal(nrow(db), nrow(testdata))
  expect_equal(ncol(db), length(coef(fit)))
})

test_that("dfbeta by country works on full data", {
  skip_on_cran()
  fit <- quick_fit(testdata, gamma = 0.005, countries = ~country)
  db <- dfbeta(fit, by = "country")
  expect_equal(nrow(db), length(unique(testdata$country)))
})

test_that("dfpopsize returns correct length on full data", {
  skip_on_cran()
  fit <- quick_fit(testdata, gamma = 0.005)
  dp <- dfpopsize(fit)
  expect_true(is.numeric(dp))
  expect_equal(length(dp), nrow(testdata))
  expect_true(!is.null(names(dp)))
  expect_true(sum(is.finite(dp)) > nrow(testdata) * 0.8)
})

test_that("dfpopsize by country on full data", {
  skip_on_cran()
  fit <- quick_fit(testdata, gamma = 0.005, countries = ~country)
  dp <- dfpopsize(fit, by = "country")
  expect_equal(length(dp), length(unique(testdata$country)))
})

test_that("dfbeta with NB on full data", {
  skip_on_cran()
  fit <- quick_fit(testdata, method = "nb", gamma = 0.005)
  db <- dfbeta(fit)
  expect_equal(nrow(db), nrow(testdata))
  expect_equal(ncol(db), length(coef(fit)))
})
