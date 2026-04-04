# ---- Population size estimation tests (CRAN) ----

test_that("popsize returns correct structure", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ps <- popsize(fit)
  expect_true(is.data.frame(ps))
  expect_true(all(c("group", "estimate", "estimate_bc", "lower", "upper") %in% names(ps)))
  expect_true(nrow(ps) >= 1)
})

test_that("bias correction reduces estimate (Jensen's inequality)", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ps <- popsize(fit, bias_correction = TRUE)
  expect_true(all(ps$estimate_bc <= ps$estimate))
})

test_that("BC estimate is positive and <= plugin for Poisson and NB", {
  d <- small_data()
  for (method in c("poisson", "nb")) {
    fit <- quick_fit(d, method = method, gamma = 0.005)
    ps <- popsize(fit, bias_correction = TRUE)
    expect_true(all(ps$estimate_bc > 0), info = paste("Method:", method))
    expect_true(all(ps$estimate_bc <= ps$estimate), info = paste("Method:", method))
  }
})

test_that("no BC when bias_correction=FALSE: estimate == estimate_bc", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ps <- popsize(fit, bias_correction = FALSE)
  expect_equal(ps$estimate, ps$estimate_bc)
})

test_that("popsize with covariates gives multiple rows", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~0 + sex)
  ps <- popsize(fit)
  expect_equal(nrow(ps), 2)
})

test_that("popsize CI is ordered", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ps <- popsize(fit)
  expect_true(all(ps$lower < ps$estimate))
  expect_true(all(ps$estimate < ps$upper))
})

test_that("share_pct column sums to 100", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  ps <- popsize(fit)
  expect_true("share_pct" %in% names(ps))
  expect_equal(sum(ps$share_pct), 100, tolerance = 0.01)
})

test_that("delta-method total CI for multi-group models", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  ps <- popsize(fit, total = TRUE)
  total <- attr(ps, "total")
  expect_false(is.null(total))
  expect_true(total$estimate > 0)
  expect_true(total$se > 0)
  expect_true(total$lower < total$estimate)
  expect_true(total$upper > total$estimate)
  expect_equal(total$estimate, sum(ps$estimate), tolerance = 0.01)
})

test_that("popsize default has no total", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  ps <- popsize(fit)
  expect_null(attr(ps, "total"))
})

test_that("popsize by aggregates correctly", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ps_by <- popsize(fit, by = ~ sex, total = TRUE)
  ps_all <- popsize(fit, total = TRUE)

  total_by <- attr(ps_by, "total")$estimate
  total_attr <- attr(ps_all, "total")
  total_all <- if (!is.null(total_attr)) total_attr$estimate else sum(ps_all$estimate)
  expect_equal(total_by, total_all, tolerance = 1e-6)
  expect_equal(nrow(ps_by), 2)
})

test_that("popsize by: sum of group estimates = total", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ps <- popsize(fit, by = ~ sex, total = TRUE)
  expect_equal(sum(ps$estimate), attr(ps, "total")$estimate, tolerance = 1e-6)
})

test_that("popsize by: observed column sums m per group", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ps <- popsize(fit, by = ~ sex)
  expect_equal(ps$observed[ps$group == "M"], sum(d$m[d$sex == "M"]))
  expect_equal(ps$observed[ps$group == "F"], sum(d$m[d$sex == "F"]))
})

test_that("popsize by errors on missing variable", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  expect_error(popsize(fit, by = ~ nonexistent), "not found")
})

test_that("popsize by has share_pct summing to 100", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ps <- popsize(fit, by = ~ sex)
  expect_equal(sum(ps$share_pct), 100, tolerance = 0.01)
})

test_that("popsize by works with NB model", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005)
  ps <- popsize(fit, by = ~ sex)
  expect_equal(nrow(ps), 2)
  expect_true(all(ps$estimate > 0))
})

test_that("profile_gamma returns correct structure", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  result <- profile_gamma(fit, gamma_grid = c(0.001, 0.01, 0.1), plot = FALSE)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("gamma", "xi", "loglik") %in% names(result)))
  expect_equal(nrow(result), 3)
  expect_true(all(result$xi > 0))
  expect_true(all(is.finite(result$loglik)))
})

# ---- Regression: printed total CI uses delta-method, not summed bounds ----

test_that("total CI uses delta-method, not summed subgroup bounds", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  ps <- popsize(fit, total = TRUE)
  tot <- attr(ps, "total")

  # Delta-method bounds should differ from naive sums (correlated groups)
  summed_lower <- sum(ps$lower)
  summed_upper <- sum(ps$upper)
  expect_false(isTRUE(all.equal(tot$lower, summed_lower)),
               info = "Total lower should NOT equal sum of subgroup lowers")
  expect_false(isTRUE(all.equal(tot$upper, summed_upper)),
               info = "Total upper should NOT equal sum of subgroup uppers")

  # Delta-method CI should be wider than summed CI (positive correlation)
  expect_true(tot$upper - tot$lower > summed_upper - summed_lower,
              info = "Delta-method total CI should be wider than summed CI")
})
