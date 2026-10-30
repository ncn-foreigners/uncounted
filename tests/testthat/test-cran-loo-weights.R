# ---- loo(): refits use the original model, including its weights (CRAN) ----

test_that("loo() and dfbeta() keep the fit's weights", {
  w <- ifelse(testdata$N > median(testdata$N), 1, 0.01)
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005, weights = w)
  l <- loo(fit)
  manual <- estimate_hidden_pop(testdata[-1, ], ~m, ~n, ~N, gamma = 0.005, weights = w[-1])
  unweighted <- estimate_hidden_pop(testdata[-1, ], ~m, ~n, ~N, gamma = 0.005)
  expect_equal(unname(l$coefficients[1, ]), unname(coef(manual)), tolerance = 1e-5)
  expect_false(isTRUE(all.equal(unname(l$coefficients[1, ]), unname(coef(unweighted)),
                                tolerance = 1e-4)))
})

test_that("loo() by country drops the weights of the whole country", {
  w <- rep(c(1, 3, 0.5, 2), length.out = nrow(testdata))
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005, weights = w,
                             countries = ~country)
  l <- loo(fit, by = "country")
  keep <- testdata$country != l$dropped[1]
  manual <- estimate_hidden_pop(testdata[keep, ], ~m, ~n, ~N, gamma = 0.005,
                                weights = w[keep], countries = ~country)
  expect_equal(unname(l$coefficients[1, ]), unname(coef(manual)), tolerance = 1e-5)
})

test_that("a refit that does not converge is reported as such", {
  fit <- quick_fit()
  local_mocked_bindings(
    estimate_hidden_pop = function(...) {
      out <- fit
      out$convergence <- 1L
      out
    }
  )
  l <- loo(fit)
  expect_false(any(l$converged))
  expect_true(all(is.na(l$dxi)))
})

# ---- dxi: influence with the target held fixed ----

test_that("dxi applies the refit's coefficients to all rows", {
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005, countries = ~country)
  l <- loo(fit, by = "country")
  keep <- testdata$country != l$dropped[1]
  refit <- estimate_hidden_pop(testdata[keep, ], ~m, ~n, ~N, gamma = 0.005,
                               countries = ~country)
  xi_full <- sum(popsize(fit)$estimate)
  expect_equal(l$dxi[1], sum(testdata$N^coef(refit)[1]) - xi_full, tolerance = 1e-6)
  expect_equal(l$dxi_subset[1], sum(popsize(refit)$estimate) - xi_full, tolerance = 1e-6)
  expect_equal(unname(dfpopsize(fit, by = "country")), l$dxi)

  ## dxi_subset includes each unit's own share, so it is mostly negative;
  ## influence on a fixed target is not
  expect_gt(mean(l$dxi_subset < 0), 0.75)
  expect_lt(mean(l$dxi < 0), 0.75)
})

test_that("dxi handles covariates and the constrained scale", {
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005, cov_alpha = ~sex,
                             constrained = TRUE, countries = ~country)
  l <- loo(fit, by = "country")
  keep <- testdata$country != l$dropped[1]
  refit <- estimate_hidden_pop(testdata[keep, ], ~m, ~n, ~N, gamma = 0.005,
                               cov_alpha = ~sex, constrained = TRUE, countries = ~country)
  a <- plogis(as.numeric(fit$X_alpha %*% refit$alpha_coefs))
  expect_equal(l$dxi[1], sum(testdata$N^a) - sum(popsize(fit)$estimate), tolerance = 1e-6)
})
