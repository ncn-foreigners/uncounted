# ---- OLS with estimated gamma (CRAN) ----

test_that("tidy() gamma SE for OLS is on the gamma scale", {
  ## The Jacobian was d log(mu) / d gamma while tidy() assumes the log-gamma
  ## scale, so the SE came out about 29x too small on the fixture.
  d <- positive_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "ols")
  row <- tidy(fit, conf.int = TRUE)
  row <- row[row$term == "gamma", ]

  ## Numerical SE of log(gamma) from the profiled Gaussian log-likelihood
  y <- log(d$m)
  profile_ll <- function(lg) {
    Z <- cbind(log(d$N), log(exp(lg) + d$n / d$N))
    rss <- sum(lm.fit(Z, y)$residuals^2)
    -nrow(d) / 2 * log(rss)
  }
  h <- 1e-3
  lg0 <- log(fit$gamma)
  curv <- -(profile_ll(lg0 + h) - 2 * profile_ll(lg0) + profile_ll(lg0 - h)) / h^2
  se_model <- fit$gamma / sqrt(curv)

  ## vcov_full is a sandwich, so same order of magnitude rather than equal
  expect_gt(row$std.error, se_model / 3)
  expect_lt(row$std.error, se_model * 3)
  expect_equal(row$std.error, fit$gamma * sqrt(fit$vcov_full["gamma", "gamma"]))
  expect_gt(row$conf.low, 0)
})

test_that("coefficient vcov for OLS is not affected by the gamma scale", {
  fit <- estimate_hidden_pop(positive_data(), ~m, ~n, ~N, method = "ols")
  expect_equal(vcov(fit), fit$vcov_full[1:2, 1:2], ignore_attr = TRUE)
  expect_true(all(is.finite(vcov(fit))))
})
