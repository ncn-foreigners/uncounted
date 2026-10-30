# ---- Input checks that used to end in cryptic errors or silent results (CRAN) ----

test_that("NA in the cluster variable is an error for every method", {
  ## NB used to treat NA as one more cluster, with only a rowsum() warning
  d <- testdata
  d$country[1:3] <- NA
  for (method in c("poisson", "nb", "ols", "nls")) {
    dd <- if (method == "ols") positive_data(d) else d
    expect_error(
      estimate_hidden_pop(dd, ~m, ~n, ~N, method = method, gamma = 0.005,
                          cluster = ~country, vcov = "HC1"),
      "'cluster' .* missing value"
    )
  }
})

test_that("gamma = 0 with zero auxiliary counts gives a clear error", {
  ## Used to fail with "NA/NaN/Inf in 'x'" from lm.fit
  expect_true(any(testdata$n == 0))
  expect_error(estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0), "'gamma = 0'.*zero")
  ## Still allowed when every auxiliary count is positive
  fit <- estimate_hidden_pop(positive_data(), ~m, ~n, ~N, gamma = 0)
  expect_s3_class(fit, "uncounted")
})
