# ---- All methods on full data (extended, skip on CRAN) ----

test_that("all methods with estimated gamma on full data", {
  skip_on_cran()
  for (method in c("ols", "poisson", "nb")) {
    fit <- quick_fit(testdata, method = method, gamma = "estimate")
    expect_s3_class(fit, "uncounted")
    expect_true(!is.null(fit$gamma))
    expect_true(fit$gamma > 0)
  }
})

test_that("all methods with covariates on full data", {
  skip_on_cran()
  for (method in c("ols", "poisson", "nb")) {
    fit <- quick_fit(testdata, method = method, gamma = 0.005,
                     cov_alpha = ~0 + sex, cov_beta = ~0 + sex)
    expect_equal(length(coef(fit)), 4)
    expect_true(all(is.finite(coef(fit))))
  }
})

test_that("DGP recovery: alpha estimate close to truth", {
  skip_on_cran()
  fit <- quick_fit(testdata, gamma = dgp$gamma)
  alpha_hat <- coef(fit)[["alpha"]]
  expect_true(abs(alpha_hat - dgp$alpha) < 0.2,
              info = paste("alpha_hat =", round(alpha_hat, 3),
                           "vs truth =", dgp$alpha))
})

test_that("DGP recovery: estimated gamma in reasonable range", {
  skip_on_cran()
  fit <- quick_fit(testdata, gamma = "estimate")
  expect_true(fit$gamma > 0.0001 & fit$gamma < 0.1,
              info = paste("gamma_hat =", fit$gamma,
                           "vs truth =", dgp$gamma))
})

test_that("Poisson on full data: base R comparison", {
  skip_on_cran()
  d <- positive_data()
  g <- 0.005

  glm_fit <- glm(m ~ -1 + log(N) + log(g + n / N), family = poisson(), data = d)
  our_fit <- quick_fit(d, gamma = g)

  expect_equal(as.numeric(coef(our_fit)), as.numeric(coef(glm_fit)),
               tolerance = 1e-5)

  V_glm <- sandwich::vcovHC(glm_fit, type = "HC3")
  V_our <- vcov(our_fit)
  expect_equal(as.numeric(V_our), as.numeric(V_glm), tolerance = 1e-4)
})

test_that("OLS on full positive data: base R comparison", {
  skip_on_cran()
  d <- positive_data()
  g <- 0.005

  lm_fit <- lm(log(m) ~ 0 + log(N) + log(g + n / N), data = d)
  our_fit <- quick_fit(d, method = "ols", gamma = g)

  expect_equal(as.numeric(coef(our_fit)), as.numeric(coef(lm_fit)),
               tolerance = 1e-5)
})

test_that("Poisson with covariates on full data: base R comparison", {
  skip_on_cran()
  d <- positive_data()

  glm_fit <- glm(m ~ -1 + sex:log(N) + log(n / N), family = poisson(), data = d)
  our_fit <- quick_fit(d, gamma = NULL, cov_alpha = ~0 + sex)

  expect_equal(our_fit |> coef() |> (\(x) x[["alpha:sexF"]])(),
               coef(glm_fit)[["sexF:log(N)"]], tolerance = 1e-5)
  expect_equal(our_fit |> coef() |> (\(x) x[["alpha:sexM"]])(),
               coef(glm_fit)[["sexM:log(N)"]], tolerance = 1e-5)
})

test_that("popsize on full data: all methods give positive estimates", {
  skip_on_cran()
  for (method in c("poisson", "nb")) {
    fit <- quick_fit(testdata, method = method, gamma = 0.005)
    ps <- popsize(fit, bias_correction = TRUE)
    expect_true(all(ps$estimate > 0), info = paste("Method:", method))
    expect_true(all(ps$estimate_bc > 0), info = paste("Method:", method))
    expect_true(all(ps$lower < ps$estimate), info = paste("Method:", method))
    expect_true(all(ps$estimate < ps$upper), info = paste("Method:", method))
  }
})
