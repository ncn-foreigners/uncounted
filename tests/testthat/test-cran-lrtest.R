# ---- lrtest(): conditions of the likelihood ratio test (CRAN) ----

test_that("lrtest() requires the same data and weights", {
  ## Only n_obs used to be compared: weights of 2 gave LR = 1376 for a covariate
  ## that does nothing
  f0 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005)
  f1 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~sex, gamma = 0.005)
  fw <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005,
                            weights = rep(2, nrow(testdata)))
  expect_error(lrtest(fw, f1), "differ in: weights")

  d <- testdata
  d$m <- rev(d$m)
  fy <- estimate_hidden_pop(d, ~m, ~n, ~N, gamma = 0.005)
  expect_error(lrtest(fy, f1), "differ in: observed")

  ## Unit weights are the same as no weights; equal weights pass
  fu <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005,
                            weights = rep(1, nrow(testdata)))
  expect_no_error(lrtest(fu, f1))
  fw1 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~sex, gamma = 0.005,
                             weights = rep(2, nrow(testdata)))
  expect_no_error(lrtest(fw, fw1))
})
