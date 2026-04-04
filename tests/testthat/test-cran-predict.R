# ---- predict.uncounted tests ----

test_that("predict returns fitted values with no newdata", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  pred <- predict(fit)
  expect_equal(pred, fitted(fit))
})

test_that("predict with type='link' returns log scale", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  pred_link <- predict(fit, type = "link")
  pred_resp <- predict(fit, type = "response")
  expect_equal(exp(pred_link), pred_resp, tolerance = 1e-10)
})

test_that("predict with newdata reproduces fitted values", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  pred_orig <- predict(fit, newdata = d)
  expect_equal(pred_orig, fitted(fit), tolerance = 1e-10)
})

test_that("predict with newdata subset works", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  pred_sub <- predict(fit, newdata = d[1:5, ])
  expect_length(pred_sub, 5)
  expect_equal(pred_sub, fitted(fit)[1:5], tolerance = 1e-10)
})

test_that("predict for all methods", {
  d <- small_data()
  for (method in c("ols", "poisson", "nb", "iols")) {
    suppressWarnings({
      fit <- quick_fit(d, method = method, gamma = 0.005)
    })
    pred <- predict(fit)
    expect_length(pred, nrow(d))
    expect_true(all(is.finite(pred)), info = method)
  }
})

test_that("predict with covariates reproduces fitted", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  pred <- predict(fit, newdata = d)
  expect_equal(pred, fitted(fit), tolerance = 1e-10)
})
