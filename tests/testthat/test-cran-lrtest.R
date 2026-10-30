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
  fy <- suppressWarnings(estimate_hidden_pop(d, ~m, ~n, ~N, gamma = 0.005))
  expect_error(lrtest(fy, f1), "differ in: observed")

  ## Unit weights are the same as no weights; equal weights pass
  fu <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005,
                            weights = rep(1, nrow(testdata)))
  expect_no_error(lrtest(fu, f1))
  fw1 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~sex, gamma = 0.005,
                             weights = rep(2, nrow(testdata)))
  expect_no_error(lrtest(fw, fw1))
})

test_that("an NB model is never treated as nested in a Poisson model", {
  ## Parameter counts used to decide, so NB(intercept) vs Poisson(~country) ran as
  ## "NB nested in Poisson" with no warning
  nb0 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  po_c <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~country, gamma = 0.005)
  expect_error(lrtest(nb0, po_c), "not nested")
  expect_error(lrtest(po_c, nb0), "not nested")
  ## Same parameter count, Poisson with the extra covariate
  po_s <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~sex, gamma = 0.005)
  expect_error(lrtest(po_s, nb0), "not nested")
})

test_that("Poisson nested in NB is accepted in either argument order", {
  po <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005)
  nb <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  nb_s <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", cov_alpha = ~sex,
                              gamma = 0.005)
  r1 <- lrtest(po, nb)
  r2 <- lrtest(nb, po)
  expect_equal(r1$statistic, r2$statistic)
  expect_equal(r1$df, 1)
  expect_match(r2$model1, "POISSON")
  expect_equal(lrtest(nb_s, po)$df, 2)
})
