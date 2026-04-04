# ---- iOLS estimator tests ----

# ---- Basic functionality ----

test_that("iOLS runs with fixed gamma", {
  d <- small_data()
  fit <- quick_fit(d, method = "iols", gamma = 0.005)
  expect_s3_class(fit, "uncounted")
  expect_equal(length(coef(fit)), 2)
  expect_true(all(is.finite(coef(fit))))
})

test_that("iOLS runs with estimated gamma", {
  d <- small_data()
  fit <- quick_fit(d, method = "iols", gamma = "estimate")
  expect_s3_class(fit, "uncounted")
  expect_true(!is.null(fit$gamma))
  expect_true(fit$gamma > 0)
})

test_that("iOLS runs with gamma=NULL (no zeros needed)", {
  d <- positive_data(small_data())
  fit <- quick_fit(d, method = "iols", gamma = NULL)
  expect_s3_class(fit, "uncounted")
  expect_true(all(is.finite(coef(fit))))
})

test_that("iOLS handles zeros in m without error", {
  d <- small_data()
  d$m[1] <- 0L
  fit <- quick_fit(d, method = "iols", gamma = 0.005)
  expect_s3_class(fit, "uncounted")
  expect_true(all(is.finite(coef(fit))))
})

# ---- Comparison with Poisson MLE (should be close for large delta) ----

test_that("iOLS with fixed gamma approximates Poisson MLE", {
  d <- small_data()
  g <- 0.005
  fit_iols <- quick_fit(d, method = "iols", gamma = g)
  fit_pois <- quick_fit(d, method = "poisson", gamma = g)

  # iOLS (which converges to GPML) should be close to PPML
  # Not identical (GPML != PPML) but in the same ballpark
  rel_diff <- abs(coef(fit_iols) - coef(fit_pois)) / (abs(coef(fit_pois)) + 0.01)
  expect_true(all(rel_diff < 0.20),
              info = paste("iOLS:", paste(round(coef(fit_iols), 4), collapse = ", "),
                           "Poisson:", paste(round(coef(fit_pois), 4), collapse = ", ")))
})

test_that("iOLS with gamma=NULL approximates Poisson on positive data", {
  d <- positive_data(small_data())
  fit_iols <- quick_fit(d, method = "iols", gamma = NULL)

  glm_fit <- glm(m ~ -1 + log(N) + log(n / N), family = poisson(), data = d)

  # GPML and PPML should be close on positive-only data
  rel_diff <- abs(coef(fit_iols) - as.numeric(coef(glm_fit))) /
    (abs(as.numeric(coef(glm_fit))) + 0.01)
  expect_true(all(rel_diff < 0.20),
              info = paste("iOLS:", paste(round(coef(fit_iols), 4), collapse = ", "),
                           "glm:", paste(round(coef(glm_fit), 4), collapse = ", ")))
})

# ---- Comparison with OLS (small delta limit) ----

test_that("iOLS coefficients differ from naive log(m+1) OLS", {
  d <- small_data()
  d$m[1] <- 0L
  g <- 0.005
  fit_iols <- quick_fit(d, method = "iols", gamma = g)
  fit_ols <- quick_fit(d, method = "ols", gamma = g)

  # They should NOT be identical (iOLS corrects the log(m+1) bias)
  expect_false(all(abs(coef(fit_iols) - coef(fit_ols)) < 1e-3),
               info = "iOLS should differ from naive OLS with zeros")
})

# ---- S3 methods ----

test_that("iOLS S3 methods work", {
  d <- small_data()
  fit <- quick_fit(d, method = "iols", gamma = 0.005)

  expect_true(is.numeric(coef(fit)))
  expect_true(is.matrix(vcov(fit)))
  expect_true(is.numeric(fitted(fit)))
  expect_true(is.numeric(residuals(fit)))
  expect_output(print(fit))
  expect_output(summary(fit))
})

test_that("iOLS popsize works", {
  d <- small_data()
  fit <- quick_fit(d, method = "iols", gamma = 0.005)
  ps <- popsize(fit)
  expect_true(is.data.frame(ps))
  expect_true(all(ps$estimate > 0))
  expect_true(all(ps$lower < ps$estimate))
  expect_true(all(ps$estimate < ps$upper))
})

# ---- Covariates ----

test_that("iOLS with cov_alpha works", {
  d <- small_data()
  fit <- quick_fit(d, method = "iols", gamma = 0.005, cov_alpha = ~0 + sex)
  expect_true("alpha:sexF" %in% names(coef(fit)))
  expect_true("alpha:sexM" %in% names(coef(fit)))
})

# ---- Sandwich vcov ----

test_that("iOLS bread and estfun work", {
  d <- small_data()
  fit <- quick_fit(d, method = "iols", gamma = 0.005)
  B <- sandwich::bread(fit)
  ef <- sandwich::estfun(fit)
  expect_equal(ncol(ef), length(coef(fit)))
  expect_equal(dim(B), rep(length(coef(fit)), 2))
})

test_that("iOLS sandwich::vcovHC produces valid vcov", {
  d <- small_data()
  fit <- quick_fit(d, method = "iols", gamma = 0.005)
  V <- sandwich::vcovHC(fit, type = "HC3")
  expect_equal(dim(V), rep(length(coef(fit)), 2))
  expect_true(all(diag(V) > 0))
})
