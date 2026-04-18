# ---- Dependence bounds tests (CRAN) ----

.baseline_total <- function(ps, bias_correction = TRUE) {
  tot <- attr(ps, "total")
  col <- if (bias_correction) "estimate_bc" else "estimate"
  if (!is.null(tot)) tot[[col]] else ps[[col]][1]
}

test_that("dependence_bounds returns expected structure", {
  fit <- quick_fit(gamma = 0.005)
  sens <- dependence_bounds(fit)

  expect_s3_class(sens, "uncounted_dependence_bounds")
  expect_true(all(c(
    "Gamma", "estimate", "lower", "upper",
    "pct_change_lower", "pct_change_upper"
  ) %in% names(sens$table)))
  expect_true(all(c("estimate", "estimate_bc", "se", "lower", "upper") %in%
                    names(sens$baseline)))
})

test_that("Gamma = 1 reproduces the baseline total for single-group fits", {
  fit <- quick_fit(gamma = 0.005)
  ps <- popsize(fit, total = TRUE)
  sens <- dependence_bounds(fit, Gamma = 1)

  expect_equal(sens$table$estimate, .baseline_total(ps, bias_correction = TRUE),
               tolerance = 1e-12)
  expect_equal(sens$table$lower, sens$table$estimate, tolerance = 1e-12)
  expect_equal(sens$table$upper, sens$table$estimate, tolerance = 1e-12)
})

test_that("Gamma values are sorted and deduplicated", {
  fit <- quick_fit(gamma = 0.005)
  sens <- dependence_bounds(fit, Gamma = c(1.25, 1.00, 1.10, 1.10))

  expect_equal(sens$table$Gamma, c(1.00, 1.10, 1.25))
})

test_that("bounds are monotone in Gamma and enclose the baseline", {
  fit <- quick_fit(gamma = 0.005)
  sens <- dependence_bounds(fit, Gamma = c(1, 1.05, 1.10, 1.25))

  expect_true(all(diff(sens$table$lower) <= 0))
  expect_true(all(diff(sens$table$upper) >= 0))
  expect_true(all(sens$table$lower <= sens$table$estimate))
  expect_true(all(sens$table$estimate <= sens$table$upper))
})

test_that("dependence_bounds validates Gamma inputs", {
  fit <- quick_fit(gamma = 0.005)

  expect_error(dependence_bounds(fit, Gamma = 0.99), ">= 1")
  expect_error(dependence_bounds(fit, Gamma = c(1, NA_real_)), "finite values >= 1")
  expect_error(dependence_bounds(fit, Gamma = c(1, Inf)), "finite values >= 1")
  expect_error(dependence_bounds(fit, Gamma = "1.1"), "numeric vector")
  expect_error(dependence_bounds(fit, Gamma = numeric(0)), "non-empty numeric vector")
})

test_that("dependence_bounds validates threshold inputs", {
  fit <- quick_fit(gamma = 0.005)

  expect_error(dependence_bounds(fit, threshold = -1), "positive finite number")
  expect_error(dependence_bounds(fit, threshold = c(1, 2)), "positive finite number")
  expect_error(dependence_bounds(fit, threshold = NA_real_), "positive finite number")
})

test_that("dependence_bounds works for Poisson and NB fits", {
  d <- small_data()

  for (method in c("poisson", "nb")) {
    fit <- quick_fit(d, method = method, gamma = 0.005)
    sens <- dependence_bounds(fit, Gamma = c(1, 1.1))

    expect_equal(nrow(sens$table), 2, info = paste("Method:", method))
    expect_true(all(sens$table$estimate > 0), info = paste("Method:", method))
  }
})

test_that("dependence_bounds uses total xi for multi-group fits", {
  fit <- quick_fit(gamma = 0.005, cov_alpha = ~ sex)
  ps <- popsize(fit, total = TRUE)
  sens <- dependence_bounds(fit, Gamma = 1)

  expect_equal(sens$table$estimate, attr(ps, "total")$estimate_bc, tolerance = 1e-12)
  expect_equal(sens$baseline$estimate, attr(ps, "total")$estimate, tolerance = 1e-12)
})

test_that("bias_correction = FALSE uses the plug-in total", {
  fit <- quick_fit(gamma = 0.005, cov_alpha = ~ sex)
  ps <- popsize(fit, bias_correction = FALSE, total = TRUE)
  sens <- dependence_bounds(fit, Gamma = 1, bias_correction = FALSE)

  expect_equal(sens$table$estimate, attr(ps, "total")$estimate, tolerance = 1e-12)
  expect_true(is.na(sens$baseline$estimate_bc))
})

test_that("threshold handling returns finite gamma_star when crossed", {
  fit <- quick_fit(gamma = 0.005)
  baseline <- dependence_bounds(fit, Gamma = 1)$table$estimate

  sens_lower <- dependence_bounds(
    fit,
    Gamma = c(1, 1.1, 1.25),
    threshold = baseline * 0.95,
    threshold_side = "lower"
  )
  sens_upper <- dependence_bounds(
    fit,
    Gamma = c(1, 1.1, 1.25),
    threshold = baseline * 1.05,
    threshold_side = "upper"
  )

  expect_equal(sens_lower$gamma_star, 1.1, tolerance = 1e-12)
  expect_equal(sens_upper$gamma_star, 1.1, tolerance = 1e-12)
})

test_that("threshold handling returns NA when grid does not cross threshold", {
  fit <- quick_fit(gamma = 0.005)
  baseline <- dependence_bounds(fit, Gamma = 1)$table$estimate
  sens <- dependence_bounds(
    fit,
    Gamma = c(1, 1.1, 1.25),
    threshold = baseline * 2,
    threshold_side = "upper"
  )

  expect_true(is.na(sens$gamma_star))
})

test_that("sensitivity_dependence is a deprecated wrapper", {
  fit <- quick_fit(gamma = 0.005)
  bounds <- dependence_bounds(fit, Gamma = c(1, 1.1))
  expect_warning(
    sens <- sensitivity_dependence(fit, Gamma = c(1, 1.1)),
    "dependence_bounds"
  )

  expect_equal(sens$table, bounds$table, tolerance = 1e-12)
  expect_s3_class(sens, "uncounted_dependence_bounds")
})

test_that("print.uncounted_dependence_bounds shows a stable header", {
  fit <- quick_fit(gamma = 0.005)
  sens <- dependence_bounds(fit, Gamma = c(1, 1.1))

  expect_output(print(sens), "Dependence bounds analysis")
})
