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

# ---- Regression: factor encoding safety ----

test_that("predict with reordered factor levels gives same result", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  # Reverse factor levels
  d2 <- d
  d2$sex <- factor(d2$sex, levels = rev(levels(factor(d$sex))))
  pred_orig <- predict(fit, newdata = d)
  pred_reord <- predict(fit, newdata = d2)
  expect_equal(pred_orig, pred_reord, tolerance = 1e-10)
})

test_that("predict with single-level factor subset works", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  # Subset to only M
  d_m <- d[d$sex == "M", ]
  pred <- predict(fit, newdata = d_m)
  expect_length(pred, nrow(d_m))
  # Should match fitted for those rows
  idx_m <- which(d$sex == "M")
  expect_equal(pred, fitted(fit)[idx_m], tolerance = 1e-10)
})
