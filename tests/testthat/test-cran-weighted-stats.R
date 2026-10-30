# ---- Fit statistics respect the weights the same way across methods (CRAN) ----
#
# logLik was weighted for Poisson/NB but not OLS/NLS, and deviance and Pearson X2
# were never weighted.

test_that("weights of 2 double logLik and deviance for every method", {
  d <- positive_data()
  w <- rep(2, nrow(d))
  for (method in c("poisson", "nb", "ols", "nls")) {
    f1 <- estimate_hidden_pop(d, ~m, ~n, ~N, method = method, gamma = 0.005)
    f2 <- estimate_hidden_pop(d, ~m, ~n, ~N, method = method, gamma = 0.005, weights = w)
    ## A Gaussian log-likelihood at sigma2-hat scales exactly with sum(w)
    expect_equal(as.numeric(logLik(f2)), 2 * as.numeric(logLik(f1)), tolerance = 1e-4,
                 info = method)
    expect_equal(deviance(f2), 2 * deviance(f1), tolerance = 1e-3, info = method)
  }
})

test_that("unit weights change nothing", {
  d <- positive_data()
  for (method in c("poisson", "ols", "nls")) {
    f0 <- estimate_hidden_pop(d, ~m, ~n, ~N, method = method, gamma = 0.005)
    f1 <- estimate_hidden_pop(d, ~m, ~n, ~N, method = method, gamma = 0.005,
                              weights = rep(1, nrow(d)))
    expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f0)), tolerance = 1e-6)
    expect_equal(deviance(f1), deviance(f0), tolerance = 1e-6)
  }
})

test_that("weighted logLik equals the profile criterion", {
  d <- testdata
  w <- rep(c(1, 3, 0.5, 2), length.out = nrow(d))
  for (method in c("ols", "nls")) {
    fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = method, gamma = 0.005, weights = w)
    info <- .build_nll_from_object(fit)
    expect_equal(-info$nll(info$par0), as.numeric(logLik(fit)))
  }
})

test_that("Pearson X2 in compare_models() is weighted", {
  d <- testdata
  w <- rep(2, nrow(d))
  f1 <- estimate_hidden_pop(d, ~m, ~n, ~N, gamma = 0.005)
  f2 <- estimate_hidden_pop(d, ~m, ~n, ~N, gamma = 0.005, weights = w)
  tab <- compare_models(A = f1, B = f2)$table
  expect_equal(tab$Pearson_X2[tab$Model == "B"], 2 * tab$Pearson_X2[tab$Model == "A"],
               tolerance = 1e-3)
})
