# ---- Start values under constrained = TRUE (CRAN) ----
#
# logit()/log() used to be applied to every start coefficient, so a contrast of
# ~0 was clipped to 0.01 and became -4.6 on the link scale; the optimiser then
# stopped far from the optimum with convergence = 0 (xi about 14x too small).
# When the unconstrained fit lies inside the constraint region, the constrained
# fit must reproduce it.

test_that(".constrained_start transforms the intercept only", {
  X <- cbind("alpha:(Intercept)" = 1, "alpha:z" = 0:1)
  expect_equal(.constrained_start(c(0.7, 0.003), X, qlogis), c(qlogis(0.7), 0))
  X0 <- cbind("alpha:a" = 1:0, "alpha:b" = 0:1)
  expect_equal(.constrained_start(c(0.7, 0.6), X0, qlogis), qlogis(c(0.7, 0.6)))
  X1 <- cbind(alpha = rep(1, 2))
  expect_equal(.constrained_start(0.7, X1, qlogis), qlogis(0.7))
})

test_that("constrained Poisson and NB reproduce an interior unconstrained fit", {
  specs <- list(list(cov_beta = ~year), list(cov_alpha = ~sex + year))
  for (method in c("poisson", "nb")) {
    for (s in specs) {
      args <- c(list(testdata, ~m, ~n, ~N, method = method, gamma = 0.005), s)
      unc <- do.call(estimate_hidden_pop, args)
      con <- do.call(estimate_hidden_pop, c(args, list(constrained = TRUE)))
      expect_true(all(unc$alpha_values > 0 & unc$alpha_values < 1))
      expect_true(all(unc$beta_values > 0))
      expect_equal(con$loglik, unc$loglik, tolerance = 1e-5)
      expect_equal(sum(popsize(con)$estimate), sum(popsize(unc)$estimate),
                   tolerance = 1e-3)
    }
  }
})
