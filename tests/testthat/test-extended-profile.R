## Tests for profile_alpha(), profile_beta(), profile.uncounted()
## Skip on CRAN — involves optimization at each grid point

skip_on_cran()

# ---- profile_alpha basic ----

test_that("profile_alpha returns correct data.frame", {
  d <- positive_data()
  fit <- quick_fit(d, gamma = 0.005)
  result <- profile_alpha(fit, plot = FALSE)
  expect_true(is.data.frame(result))
  expect_equal(names(result), c("value", "xi", "loglik"))
  expect_equal(nrow(result), 30)
  expect_true(all(is.finite(result$xi)))
  expect_true(all(is.finite(result$loglik)))
})

test_that("profile_alpha loglik peaks near MLE", {
  d <- positive_data()
  fit <- quick_fit(d, gamma = 0.005)
  result <- profile_alpha(fit, plot = FALSE)
  # The MLE should be near the max loglik on the grid
  best_idx <- which.max(result$loglik)
  mle_val <- coef(fit)["alpha"]
  expect_true(abs(result$value[best_idx] - mle_val) < 0.1,
              info = paste("Grid max at", round(result$value[best_idx], 3),
                           "MLE at", round(mle_val, 3)))
})

test_that("profile_alpha with custom grid works", {
  d <- positive_data()
  fit <- quick_fit(d, gamma = 0.005)
  grid <- seq(0.3, 1.0, by = 0.1)
  result <- profile_alpha(fit, grid = grid, plot = FALSE)
  expect_equal(nrow(result), length(grid))
  expect_equal(result$value, grid)
})

test_that("profile_alpha with coef_index = 2 on cov_alpha model", {
  d <- positive_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~ sex)
  expect_true(fit$p_alpha >= 2)
  result <- profile_alpha(fit, coef_index = 2, plot = FALSE)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 30)
})

test_that("profile_alpha errors for invalid coef_index", {
  d <- positive_data()
  fit <- quick_fit(d, gamma = 0.005)
  expect_error(profile_alpha(fit, coef_index = 5, plot = FALSE), "coef_index")
})

# ---- profile_alpha with reoptimize ----

test_that("profile_alpha reoptimize = TRUE works", {
  d <- positive_data()
  fit <- quick_fit(d, gamma = 0.005)
  result <- profile_alpha(fit, reoptimize = TRUE,
                          grid = seq(0.5, 0.9, length.out = 5),
                          plot = FALSE)
  expect_true(all(is.finite(result$loglik)))
  # Reoptimized loglik should be >= concentrated loglik at each point
  result_conc <- profile_alpha(fit,
                               grid = seq(0.5, 0.9, length.out = 5),
                               plot = FALSE)
  expect_true(all(result$loglik >= result_conc$loglik - 0.1))
})

# ---- profile_alpha with NB ----

test_that("profile_alpha works with NB method", {
  d <- positive_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005)
  result <- profile_alpha(fit, plot = FALSE)
  expect_true(is.data.frame(result))
  expect_true(all(is.finite(result$loglik)))
})

# ---- profile_alpha with constrained ----

test_that("profile_alpha works with constrained model", {
  d <- positive_data()
  fit <- quick_fit(d, gamma = 0.005, constrained = TRUE)
  result <- profile_alpha(fit, plot = FALSE)
  expect_true(is.data.frame(result))
  # Grid is on link (logit) scale
  expect_true(all(is.finite(result$xi)))
})

# ---- profile_beta basic ----

test_that("profile_beta returns correct data.frame", {
  d <- positive_data()
  fit <- quick_fit(d, gamma = 0.005)
  result <- profile_beta(fit, plot = FALSE)
  expect_true(is.data.frame(result))
  expect_equal(names(result), c("value", "xi", "loglik"))
  expect_equal(nrow(result), 30)
})

test_that("profile_beta xi is constant when reoptimize = FALSE", {
  d <- positive_data()
  fit <- quick_fit(d, gamma = 0.005)
  result <- profile_beta(fit, plot = FALSE)
  # xi = sum(N^alpha) doesn't depend on beta
  expect_equal(length(unique(round(result$xi, 6))), 1)
})

test_that("profile_beta loglik peaks near MLE", {
  d <- positive_data()
  fit <- quick_fit(d, gamma = 0.005)
  result <- profile_beta(fit, plot = FALSE)
  best_idx <- which.max(result$loglik)
  mle_val <- coef(fit)["beta"]
  expect_true(abs(result$value[best_idx] - mle_val) < 0.2)
})

# ---- profile.uncounted S3 method ----

test_that("profile.uncounted dispatches to profile_gamma", {
  d <- positive_data()
  fit <- quick_fit(d, gamma = 0.005)
  result <- profile(fit, param = "gamma",
                    gamma_grid = c(0.001, 0.01, 0.1), plot = FALSE)
  expect_true(is.data.frame(result))
  expect_true("gamma" %in% names(result))
})

test_that("profile.uncounted dispatches to profile_alpha", {
  d <- positive_data()
  fit <- quick_fit(d, gamma = 0.005)
  result <- profile(fit, param = "alpha", plot = FALSE)
  expect_true(is.data.frame(result))
  expect_true("value" %in% names(result))
})

test_that("profile.uncounted dispatches to profile_beta", {
  d <- positive_data()
  fit <- quick_fit(d, gamma = 0.005)
  result <- profile(fit, param = "beta", plot = FALSE)
  expect_true(is.data.frame(result))
  expect_true("value" %in% names(result))
})

# ---- plot output ----

test_that("profile_alpha with plot = TRUE runs without error", {
  d <- positive_data()
  fit <- quick_fit(d, gamma = 0.005)
  expect_no_error(profile_alpha(fit, plot = TRUE))
})
