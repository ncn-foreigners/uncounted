# ---- Dependence robustness tests (CRAN) ----

.synthetic_dependence_profile <- function() {
  data.frame(
    delta = c(-0.4, -0.2, 0, 0.1, 0.3),
    kappa = exp(c(-0.4, -0.2, 0, 0.1, 0.3)),
    xi = c(70, 90, 100, 112, 145)
  )
}

test_that("robustness_dependence works from a profile data frame", {
  prof <- .synthetic_dependence_profile()
  rob <- robustness_dependence(prof, threshold = 95, direction = "decrease")

  expect_s3_class(rob, "uncounted_dependence_robustness")
  expect_equal(rob$baseline_xi, 100)
  expect_equal(rob$target_xi, 95)
  expect_true(rob$reached)
  expect_equal(rob$rv_delta, -0.2)
  expect_equal(rob$rv_kappa, exp(-0.2), tolerance = 1e-12)
  expect_equal(rob$rv_Gamma, exp(0.2), tolerance = 1e-12)
})

test_that("robustness_dependence supports q-based targets", {
  prof <- .synthetic_dependence_profile()
  rob <- robustness_dependence(prof, q = 0.10, direction = "increase")

  expect_equal(rob$target_xi, 110)
  expect_equal(rob$rv_delta, 0.1)
})

test_that("robustness_dependence works from a fitted object", {
  fit <- quick_fit(small_data(), method = "poisson", gamma = 0.005,
                   cov_alpha = ~ sex)

  rob <- robustness_dependence(
    fit,
    q = 0.05,
    direction = "decrease",
    delta_grid = c(-0.2, -0.1, 0, 0.1, 0.2)
  )

  expect_s3_class(rob, "uncounted_dependence_robustness")
  expect_true(is.finite(rob$baseline_xi))
  expect_true(is.logical(rob$reached))
})

test_that("robustness_dependence returns NA values when target is not reached", {
  prof <- .synthetic_dependence_profile()
  rob <- robustness_dependence(prof, threshold = 1000, direction = "increase")

  expect_false(rob$reached)
  expect_true(is.na(rob$rv_delta))
  expect_true(is.na(rob$rv_kappa))
  expect_true(is.na(rob$rv_Gamma))
})

test_that("robustness_dependence validates supplied profiles", {
  bad_prof <- data.frame(delta = c(-0.1, 0.1), kappa = exp(c(-0.1, 0.1)),
                         xi = c(90, 110))

  expect_error(robustness_dependence(bad_prof), "delta = 0")
  expect_error(robustness_dependence(list()), "fit or a dependence-profile")
})

test_that("print.uncounted_dependence_robustness shows a stable header", {
  rob <- robustness_dependence(.synthetic_dependence_profile(), threshold = 95,
                               direction = "decrease")

  expect_output(print(rob), "Dependence robustness analysis")
})
