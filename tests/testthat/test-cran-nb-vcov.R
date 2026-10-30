# ---- NB MLE: covariance layout (CRAN) ----
#
# The optimiser's parameter vector is (alpha, beta, log-theta, gamma...). vcov()
# used to take the first p rows, so with cov_gamma the log-theta variance was
# reported for gamma:(Intercept) and the last gamma coefficient was dropped.

test_that("vcov() lines up with coef() for NB with cov_gamma", {
  for (type in c("HC0", "HC1")) {
    fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", cov_gamma = ~sex,
                               vcov = type)
    v <- vcov(fit)
    expect_identical(rownames(v), names(coef(fit)))
    expect_identical(colnames(v), names(coef(fit)))
    full <- fit$vcov_full
    expect_identical(colnames(full),
                     c("alpha", "beta", "theta", "gamma:(Intercept)", "gamma:sexM"))
    g <- c("gamma:(Intercept)", "gamma:sexM")
    expect_equal(v[g, g], full[g, g])
    expect_equal(v[c("alpha", "beta"), g], full[c("alpha", "beta"), g])
  }
})

test_that("gamma SEs match a direct sandwich for NB with cov_gamma", {
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", cov_gamma = ~sex,
                             vcov = "HC0")
  H_inv <- solve(fit$hessian_nll)
  V <- H_inv %*% crossprod(fit$score_full) %*% H_inv
  expect_equal(unname(sqrt(diag(vcov(fit)))[3:4]), unname(sqrt(diag(V))[4:5]))
})

test_that("tidy() reports a gamma SE for NB with scalar gamma", {
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb")
  row <- subset(tidy(fit, conf.int = TRUE), term == "gamma")
  expect_true(is.finite(row$std.error) && row$std.error > 0)
  expect_lt(row$conf.low, row$estimate)
  expect_gt(row$conf.high, row$estimate)
  ## Delta method from the log-gamma entry
  expect_equal(row$std.error, fit$gamma * sqrt(fit$vcov_full["gamma", "gamma"]))
})

test_that("print() reports gamma SEs for NB with cov_gamma", {
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", cov_gamma = ~sex)
  expect_false(any(grepl("SE(gamma) = NA", capture.output(print(fit)), fixed = TRUE)))
})

test_that("vcov_nb() is named", {
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", cov_gamma = ~sex)
  expect_true("theta" %in% colnames(vcov_nb(fit)))
})
