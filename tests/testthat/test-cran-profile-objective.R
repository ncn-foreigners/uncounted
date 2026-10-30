# ---- profile_alpha()/profile_beta(): objective and nuisance parameters (CRAN) ----

peak_at <- function(fit, width = 0.03, block = profile_alpha) {
  est <- unname(coef(fit)[if (identical(block, profile_alpha)) 1L else fit$p_alpha + 1L])
  grid <- seq(est - width, est + width, by = 0.0005)
  p <- block(fit, grid = grid, reoptimize = TRUE, plot = FALSE)
  c(est = est, peak = p$value[which.max(p$loglik)])
}

test_that("the profile peaks at the estimate for NLS and weighted OLS", {
  ## Non-count methods used an unweighted log-scale Gaussian, whose argmax is the
  ## unweighted OLS value whatever the fit.
  nl <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nls", gamma = 0.005)
  x <- peak_at(nl)
  expect_equal(x[["peak"]], x[["est"]], tolerance = 1e-3)

  w <- rep(c(1, 3, 0.5, 2), length.out = nrow(testdata))
  ols <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "ols", gamma = 0.005, weights = w)
  x <- peak_at(ols)
  expect_equal(x[["peak"]], x[["est"]], tolerance = 1e-3)

  x <- peak_at(ols, block = profile_beta)
  expect_equal(x[["peak"]], x[["est"]], tolerance = 1e-3)
})

test_that("the objective at the MLE equals the fit's criterion", {
  po <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = "estimate")
  info <- .build_nll_from_object(po)
  expect_equal(-info$nll(info$par0), po$loglik)
  expect_equal(-info$nll_full(info$par_full), po$loglik)

  nb <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", cov_gamma = ~sex)
  info <- .build_nll_from_object(nb)
  expect_equal(-info$nll_full(info$par_full), nb$loglik)

  nl <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nls", gamma = 0.005)
  expect_equal(-.build_nll_from_object(nl)$nll(c(nl$alpha_coefs, nl$beta_coefs)),
               as.numeric(logLik(nl)))
})

test_that("reoptimize = TRUE also frees gamma and theta", {
  ## Holding gamma at the MLE made the curve too narrow
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N)
  a0 <- unname(coef(fit)[1]) + 0.02
  p <- profile_alpha(fit, grid = a0, reoptimize = TRUE, plot = FALSE)
  joint <- optim(c(coef(fit)[2], log(fit$gamma)), function(par) {
    mu <- testdata$N^a0 * (exp(par[2]) + testdata$n / testdata$N)^par[1]
    -sum(dpois(testdata$m, mu, log = TRUE))
  }, method = "BFGS", control = list(reltol = 1e-12))
  expect_equal(p$loglik, -joint$value, tolerance = 1e-4)
  held <- profile_alpha(fit, grid = a0, plot = FALSE)
  expect_gt(p$loglik, held$loglik)

  nb <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb")
  a0 <- unname(coef(nb)[1]) + 0.02
  p_nb <- profile_alpha(nb, grid = a0, reoptimize = TRUE, plot = FALSE)
  expect_gt(p_nb$loglik, profile_alpha(nb, grid = a0, plot = FALSE)$loglik)
  expect_lte(p_nb$loglik, nb$loglik)
})

test_that("profiling iOLS gives a clear error", {
  fit <- quick_fit(positive_data(), method = "iols", gamma = 0.005)
  expect_error(profile_alpha(fit, plot = FALSE), "iols")
})
