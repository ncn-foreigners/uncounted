# ---- Sandwich infrastructure tests (CRAN) ----

test_that("bread method returns correct dimensions", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  B <- sandwich::bread(fit)
  p <- length(coef(fit))
  expect_equal(dim(B), c(p, p))
  expect_true(all(is.finite(B)))
  expect_equal(model.matrix(fit), fit$model_matrix_full)
  expect_equal(B, bread(fit), tolerance = 1e-12)
})

test_that("estfun method returns correct dimensions", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ef <- sandwich::estfun(fit)
  expect_equal(nrow(ef), nrow(d))
  expect_equal(ncol(ef), length(coef(fit)))
  expect_true(all(is.finite(ef)))
})

test_that("nobs method works", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  expect_equal(nobs(fit), nrow(d))
})

test_that("hatvalues method works", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  h <- hatvalues(fit)
  expect_equal(length(h), nrow(d))
  expect_true(all(h >= 0 & h <= 1))
})

test_that("sandwich::vcovHC produces valid vcov for all methods", {
  d <- small_data()
  for (method in c("ols", "poisson", "nb")) {
    fit <- quick_fit(d, method = method, gamma = 0.005)
    V <- sandwich::vcovHC(fit, type = "HC3")
    expect_equal(dim(V), rep(length(coef(fit)), 2))
    expect_true(all(diag(V) > 0))
  }
})

test_that("sandwich::vcovCL works with cluster", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, countries = ~country)
  V_cl <- sandwich::vcovCL(fit, cluster = d$country, type = "HC1")
  expect_equal(dim(V_cl), rep(length(coef(fit)), 2))
  expect_true(all(diag(V_cl) > 0))
})

test_that("sandwich::vcovCL can return the meat matrix only", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, countries = ~country)
  meat <- sandwich::vcovCL(fit, cluster = d$country, type = "HC1",
                           sandwich = FALSE, fix = FALSE)
  expect_equal(dim(meat), rep(length(coef(fit)), 2))
  expect_true(all(is.finite(meat)))
})

test_that("cluster parameter in estimate_hidden_pop works", {
  d <- small_data()
  fit_hc <- quick_fit(d, gamma = 0.005, vcov = "HC3")
  fit_cl <- quick_fit(d, gamma = 0.005, vcov = "HC1", cluster = ~country)
  se_hc <- sqrt(diag(fit_hc$vcov))
  se_cl <- sqrt(diag(fit_cl$vcov))
  expect_false(all(abs(se_hc - se_cl) < 1e-10))
})

test_that("countries does NOT trigger clustering", {
  d <- small_data()
  fit_no_cl <- quick_fit(d, gamma = 0.005, countries = ~country)
  fit_with_cl <- quick_fit(d, gamma = 0.005, countries = ~country, cluster = ~sex)
  se_no <- sqrt(diag(fit_no_cl$vcov))
  se_cl <- sqrt(diag(fit_with_cl$vcov))
  expect_false(all(abs(se_no - se_cl) < 1e-10))
})

test_that("clustered covariance requires at least two clusters", {
  d <- small_data()
  d$one_cluster <- "A"

  expect_error(
    quick_fit(d, gamma = 0.005, vcov = "HC1", cluster = ~one_cluster),
    "at least 2 clusters"
  )
  expect_error(
    quick_fit(d, method = "nb", gamma = 0.005, vcov = "HC1",
              cluster = ~one_cluster),
    "at least 2 clusters"
  )
})

test_that("vcov as function works", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005,
                   vcov = function(x) sandwich::vcovHC(x, type = "HC0"))
  expect_true(all(diag(fit$vcov) > 0))
})

test_that("vcov as function with vcovCL works", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, countries = ~country,
                   vcov = function(x) sandwich::vcovCL(x,
                     cluster = x$data$country, type = "HC1"))
  expect_true(all(diag(fit$vcov) > 0))
})

test_that("GMM covariance wrappers default to HC1", {
  d <- small_data()
  fit <- quick_fit(d, method = "poisson", gamma = 0.005, estimator = "gmm")
  V_cl_default <- sandwich::vcovCL(fit, cluster = d$country)

  expect_equal(
    sandwich::vcovHC(fit),
    sandwich::vcovHC(fit, type = "HC1"),
    tolerance = 1e-10
  )
  expect_equal(dim(V_cl_default), rep(length(coef(fit)), 2))
  expect_true(all(is.finite(V_cl_default)))
})

test_that("HC0 and HC3 produce different SE", {
  d <- small_data()
  fit_hc0 <- quick_fit(d, gamma = 0.005, vcov = "HC0")
  fit_hc3 <- quick_fit(d, gamma = 0.005, vcov = "HC3")
  expect_false(all(abs(sqrt(diag(fit_hc0$vcov)) -
                       sqrt(diag(fit_hc3$vcov))) < 1e-10))
})

test_that("HC0 and HC3 produce different SE with clustering", {
  d <- small_data()
  fit_hc0 <- quick_fit(d, gamma = 0.005, vcov = "HC0", cluster = ~country)
  fit_hc3 <- quick_fit(d, gamma = 0.005, vcov = "HC3", cluster = ~country)
  expect_false(all(abs(sqrt(diag(fit_hc0$vcov)) -
                       sqrt(diag(fit_hc3$vcov))) < 1e-10))
})

test_that("NB bread/estfun methods work", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005)
  B <- sandwich::bread(fit)
  ef <- sandwich::estfun(fit)
  expect_equal(ncol(ef), length(coef(fit)))
  expect_equal(dim(B), rep(length(coef(fit)), 2))
})

test_that("estimated gamma: bread/estfun include gamma column", {
  d <- small_data()
  fit <- quick_fit(d, gamma = "estimate")
  p_ab <- length(coef(fit))
  expect_equal(ncol(fit$model_matrix_full), p_ab + 1)

  B <- sandwich::bread(fit)
  ef <- sandwich::estfun(fit)
  expect_equal(ncol(B), p_ab + 1)
  expect_equal(ncol(ef), p_ab + 1)
  expect_equal(dim(fit$vcov), rep(p_ab, 2))
})

test_that("fwb::vcovFWB works with uncounted objects", {
  skip_if_not_installed("fwb")
  set.seed(123)
  d <- small_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N,
                              method = "poisson", gamma = 0.005)
  V_fwb <- fwb::vcovFWB(fit, R = 30)
  expect_equal(dim(V_fwb), rep(length(coef(fit)), 2))
  expect_true(all(diag(V_fwb) > 0))
})

test_that("vcov_nb warns and returns stored backend covariance for non-MLE NB fits", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005, estimator = "el")

  expect_warning(
    V_full <- vcov_nb(fit, vcov_type = "HC0", cluster = d$country),
    "ignored"
  )
  expect_equal(
    V_full,
    fit$backend_vcov_full,
    tolerance = 1e-10
  )
})
