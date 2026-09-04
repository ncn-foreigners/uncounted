## Oracle tests: pin each estimator to an independent implementation.
##
## With gamma fixed and link_rho = "power", the mean model is exactly a
## log-link GLM with two offsets-free regressors,
##
##   log mu_i = alpha * log(N_i) + beta * log(gamma + n_i / N_i),
##
## so glm(), MASS::glm.nb(), lm() and nls() are exact references for the
## coefficients, and sandwich:: for the covariance. Most of the estimation
## tests elsewhere only assert that results are finite and correctly shaped,
## which cannot catch a wrong link, a dropped chain-rule factor, or a
## mis-scaled sandwich. These can.

oracle_design <- function(data = testdata, gamma = 0.005) {
  data.frame(
    m = data$m,
    lN = log(data$N),
    lr = log(gamma + data$n / data$N)
  )
}

test_that("Poisson PMLE reproduces glm(family = poisson) coefficients", {
  gam <- 0.005
  fit <- quick_fit(testdata, method = "poisson", gamma = gam)
  ref <- glm(m ~ 0 + lN + lr, family = poisson(),
             data = oracle_design(testdata, gam))

  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
               tolerance = 1e-8)
  expect_equal(unname(fitted(fit)), unname(fitted(ref)), tolerance = 1e-5)
})

test_that("Poisson robust covariance reproduces sandwich::vcovHC on the glm", {
  gam <- 0.005
  ref <- glm(m ~ 0 + lN + lr, family = poisson(),
             data = oracle_design(testdata, gam))

  for (type in c("HC0", "HC1", "HC2", "HC3", "HC4", "HC4m", "HC5")) {
    fit <- quick_fit(testdata, method = "poisson", gamma = gam, vcov = type)
    expect_equal(
      unname(vcov(fit)),
      unname(sandwich::vcovHC(ref, type = type)),
      tolerance = 1e-4, info = paste("vcov_type =", type)
    )
  }

  ## The model-based covariance must be the ordinary Fisher one, since it is
  ## what the bias correction uses.
  fit <- quick_fit(testdata, method = "poisson", gamma = gam)
  expect_equal(unname(fit$vcov_model), unname(vcov(ref)), tolerance = 1e-5)

  ## hatvalues() must match too, or HC2/HC3 would be silently wrong.
  expect_equal(unname(hatvalues(fit)), unname(hatvalues(ref)), tolerance = 1e-6)
})

test_that("Poisson with covariates reproduces the corresponding glm", {
  gam <- 0.005
  d <- oracle_design(testdata, gam)
  ## alpha_i = a0 + a1 * 1{sex == "M"} multiplies log(N), so the equivalent
  ## glm regressors are lN and lN * sexM. Built explicitly to keep the
  ## coefficient order identical to the fit.
  d$lN_M <- d$lN * as.integer(testdata$sex == "M")
  fit <- quick_fit(testdata, method = "poisson", gamma = gam, cov_alpha = ~ sex)
  ref <- glm(m ~ 0 + lN + lN_M + lr, family = poisson(), data = d)

  expect_equal(names(coef(fit)), c("alpha:(Intercept)", "alpha:sexM", "beta"))
  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
               tolerance = 1e-8)
  expect_equal(unname(fitted(fit)), unname(fitted(ref)), tolerance = 1e-5)
})

test_that("weighted Poisson reproduces the weighted glm", {
  gam <- 0.005
  w <- rep(c(1, 1.5, 2, 0.7), length.out = nrow(testdata))
  d <- oracle_design(testdata, gam)
  fit <- quick_fit(testdata, method = "poisson", gamma = gam, weights = w)
  ref <- glm(m ~ 0 + lN + lr, family = poisson(), data = d, weights = w)

  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-5)
  expect_equal(unname(vcov(fit)), unname(sandwich::vcovHC(ref, type = "HC3")),
               tolerance = 1e-3)
})

test_that("NB MLE reproduces MASS::glm.nb", {
  gam <- 0.005
  fit <- quick_fit(testdata, method = "nb", gamma = gam)
  ref <- suppressWarnings(
    MASS::glm.nb(m ~ 0 + lN + lr, data = oracle_design(testdata, gam))
  )

  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-5)
  expect_equal(unname(fit$theta), unname(ref$theta), tolerance = 1e-3)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
               tolerance = 1e-6)
  ## Mean-parameter block of the inverse Fisher information.
  expect_equal(unname(fit$vcov_model), unname(vcov(ref)), tolerance = 1e-4)
})

test_that("OLS reproduces lm() on the log scale", {
  gam <- 0.005
  d <- oracle_design(testdata, gam)
  ## log(m + 1) is the documented handling when the response has zeros.
  expect_true(any(d$m == 0))
  d$y <- log(d$m + 1)

  fit <- quick_fit(testdata, method = "ols", gamma = gam)
  ref <- lm(y ~ 0 + lN + lr, data = d)

  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-10)
  expect_equal(unname(vcov(fit)),
               unname(sandwich::vcovHC(ref, type = "HC3")), tolerance = 1e-10)
  expect_equal(unname(fit$vcov_model), unname(vcov(ref)), tolerance = 1e-10)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
               tolerance = 1e-10)
})

test_that("NLS reproduces nls() on the original scale", {
  gam <- 0.005
  d <- oracle_design(testdata, gam)
  fit <- quick_fit(testdata, method = "nls", gamma = gam)
  ref <- nls(m ~ exp(a * lN + b * lr), data = d,
             start = list(a = 0.7, b = 0.5))

  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-4)
  ## Same residual sum of squares on the raw scale, which is what NLS minimises.
  expect_equal(sum(residuals(fit)^2), sum(residuals(ref)^2), tolerance = 1e-4)
  expect_equal(unname(fit$vcov_model), unname(vcov(ref)), tolerance = 1e-4)
})

test_that("iOLS solves the Gamma PML score and matches a Gamma GLM", {
  gam <- 0.005
  d_pos <- positive_data(testdata)
  d <- oracle_design(d_pos, gam)

  fit <- estimate_hidden_pop(d_pos, ~ m, ~ n, ~ N, method = "iols",
                             gamma = gam)
  ref <- glm(m ~ 0 + lN + lr, family = Gamma(link = "log"), data = d)

  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-4)

  ## The defining property: sum((m/mu - 1) * z_i) = 0 at the solution.
  mu <- fitted(fit)
  Z <- cbind(d$lN, d$lr)
  score <- colSums((d_pos$m / mu - 1) * Z)
  expect_equal(unname(score), c(0, 0), tolerance = 1e-2)
})

test_that("the estimated gamma is a genuine optimum of its criterion", {
  ## Guards the profiling/joint optimisation of gamma: perturbing gamma away
  ## from the fitted value must not improve the criterion.
  fit <- quick_fit(testdata, method = "poisson", gamma = "estimate")
  g_hat <- fit$gamma
  expect_true(is.finite(g_hat) && g_hat > 0)

  ll_at <- function(g) {
    as.numeric(logLik(quick_fit(testdata, method = "poisson", gamma = g)))
  }
  best <- as.numeric(logLik(fit))
  expect_gte(best, ll_at(g_hat * 0.5) - 1e-6)
  expect_gte(best, ll_at(g_hat * 2) - 1e-6)
})
