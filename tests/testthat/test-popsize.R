# ---- Test data (shared) ----

make_test_data <- function() {
  data.frame(
    m = c(14L, 1L, 56L, 21L, 29L, 57L, 1L, 41L, 258L, 6L,
          19L, 9L, 3L, 124L, 4L, 15L, 30L, 10L, 7L, 686L),
    n = c(2L, 1L, 27L, 9L, 85L, 82L, 1L, 16L, 586L, 2L,
          38L, 1L, 1L, 37L, 4L, 12L, 19L, 6L, 7L, 774L),
    N = c(1023L, 188L, 14655L, 892L, 3672L, 4473L, 44L, 5081L, 131365L, 1357L,
          4391L, 320L, 309L, 6422L, 482L, 1119L, 1784L, 746L, 2732L, 25576L),
    sex = rep(c("M", "F"), each = 10),
    country = paste0("C", 1:20)
  )
}

# ---- BC positive and less than plugin ----

test_that("BC estimate is positive and <= plugin for Poisson and NB", {
  d <- make_test_data()
  for (method in c("poisson", "nb")) {
    fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = method, gamma = 0.001)
    ps <- popsize(fit, bias_correction = TRUE)
    expect_true(all(ps$estimate_bc > 0),
                info = paste("Method:", method))
    expect_true(all(ps$estimate_bc <= ps$estimate),
                info = paste("Method:", method))
  }
})

# ---- No BC when bias_correction=FALSE ----

test_that("no BC when bias_correction=FALSE: estimate == estimate_bc", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)
  ps <- popsize(fit, bias_correction = FALSE)
  expect_equal(ps$estimate, ps$estimate_bc)
})

# ---- Group share adds up to 100% ----

test_that("share_pct column sums to 100", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                              gamma = 0.001, cov_alpha = ~sex)
  ps <- popsize(fit)
  expect_true("share_pct" %in% names(ps))
  expect_equal(sum(ps$share_pct), 100, tolerance = 0.01)
})

# ---- Delta-method total CI ----

test_that("delta-method total CI is attached for multi-group models", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                              gamma = 0.001, cov_alpha = ~sex)
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
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                              gamma = 0.001, cov_alpha = ~sex)
  ps <- popsize(fit)
  expect_null(attr(ps, "total"))
})

test_that("delta-method total CI not attached for single-group models", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)
  ps <- popsize(fit)
  expect_null(attr(ps, "total"))
})

# ---- Stratified popsize (by parameter) ----

test_that("popsize with by aggregates correctly", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.005)

  ps_by <- popsize(fit, by = ~ sex, total = TRUE)
  ps_all <- popsize(fit, total = TRUE)

  # Total should match
  total_by <- attr(ps_by, "total")$estimate
  total_attr <- attr(ps_all, "total")
  total_all <- if (!is.null(total_attr)) total_attr$estimate else sum(ps_all$estimate)
  expect_equal(total_by, total_all, tolerance = 1e-6)

  # Two groups
  expect_equal(nrow(ps_by), 2)
  expect_true(all(c("F", "M") %in% ps_by$group))
})

test_that("popsize by: sum of group estimates = total", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.005)

  ps <- popsize(fit, by = ~ sex, total = TRUE)
  expect_equal(sum(ps$estimate), attr(ps, "total")$estimate, tolerance = 1e-6)
})

test_that("popsize by: observed column sums m per group", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.005)

  ps <- popsize(fit, by = ~ sex)
  expect_equal(ps$observed[ps$group == "M"], sum(d$m[d$sex == "M"]))
  expect_equal(ps$observed[ps$group == "F"], sum(d$m[d$sex == "F"]))
})

test_that("popsize by with cov_alpha gives different grouping than default", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                              gamma = 0.005, cov_alpha = ~ sex)

  ps_default <- popsize(fit, total = TRUE)             # grouped by sex (from cov_alpha)
  ps_country <- popsize(fit, by = ~ country, total = TRUE)  # grouped by country

  expect_equal(nrow(ps_default), 2)   # M, F
  expect_equal(nrow(ps_country), 20)  # 20 countries
  # Totals must match
  expect_equal(attr(ps_default, "total")$estimate,
               attr(ps_country, "total")$estimate, tolerance = 1e-6)
})

test_that("popsize by with single group equals total", {
  d <- make_test_data()
  d$one <- "all"
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.005)

  ps <- popsize(fit, by = ~ one)
  expect_equal(nrow(ps), 1)
  total_all <- sum(popsize(fit)$estimate)
  expect_equal(ps$estimate, total_all, tolerance = 1e-6)
})

test_that("popsize by errors on missing variable", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  expect_error(popsize(fit, by = ~ nonexistent), "not found")
})

test_that("popsize by has share_pct summing to 100", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.005)

  ps <- popsize(fit, by = ~ sex)
  expect_equal(sum(ps$share_pct), 100, tolerance = 0.01)
})

test_that("popsize by works with NB model", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "nb", gamma = 0.005)

  ps <- popsize(fit, by = ~ sex)
  expect_equal(nrow(ps), 2)
  expect_true(all(ps$estimate > 0))
})

# ---- profile_gamma ----

test_that("profile_gamma returns a data frame with gamma, xi, loglik", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson")
  result <- profile_gamma(fit, gamma_grid = c(0.001, 0.01, 0.1), plot = FALSE)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("gamma", "xi", "loglik") %in% names(result)))
  expect_equal(nrow(result), 3)
  expect_true(all(result$xi > 0))
  expect_true(all(is.finite(result$loglik)))
})
