# ---- broom/modelsummary integration tests ----

test_that("tidy.uncounted returns correct structure", {
  skip_if_not_installed("generics")
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  td <- generics::tidy(fit)
  expect_true(is.data.frame(td))
  expect_true(all(c("term", "estimate", "std.error", "statistic", "p.value") %in% names(td)))
  expect_equal(nrow(td), length(coef(fit)))
})

test_that("tidy with conf.int adds CI columns", {
  skip_if_not_installed("generics")
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  td <- generics::tidy(fit, conf.int = TRUE)
  expect_true(all(c("conf.low", "conf.high") %in% names(td)))
})

test_that("glance.uncounted returns 1-row data frame", {
  skip_if_not_installed("generics")
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  gl <- generics::glance(fit)
  expect_true(is.data.frame(gl))
  expect_equal(nrow(gl), 1)
  expect_true(all(c("nobs", "logLik", "AIC", "BIC") %in% names(gl)))
})

test_that("tidy works for all methods", {
  skip_if_not_installed("generics")
  d <- small_data()
  for (method in c("ols", "poisson", "nb", "iols")) {
    suppressWarnings({
      fit <- quick_fit(d, method = method, gamma = 0.005)
    })
    td <- generics::tidy(fit)
    expect_equal(nrow(td), length(coef(fit)), info = method)
    expect_true(all(is.finite(td$estimate)), info = method)
    expect_true(all(is.finite(td$std.error)), info = method)
  }
})

test_that("glance works for all methods", {
  skip_if_not_installed("generics")
  d <- small_data()
  for (method in c("ols", "poisson", "nb", "iols")) {
    suppressWarnings({
      fit <- quick_fit(d, method = method, gamma = 0.005)
    })
    gl <- generics::glance(fit)
    expect_equal(nrow(gl), 1, info = method)
    expect_true(is.finite(gl$nobs), info = method)
  }
})

test_that("modelsummary works with tidy/glance", {
  skip_if_not_installed("modelsummary")
  skip_if_not_installed("generics")
  d <- small_data()
  fit_po <- quick_fit(d, gamma = 0.005)
  suppressWarnings({
    fit_io <- quick_fit(d, method = "iols", gamma = 0.005)
  })
  expect_no_error(
    modelsummary::modelsummary(
      list(Poisson = fit_po, iOLS = fit_io),
      output = "data.frame"
    )
  )
})
