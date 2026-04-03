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

# ---- Bootstrap structure ----

test_that("bootstrap_popsize returns correct structure", {
  skip_if_not_installed("fwb")
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)
  boot <- bootstrap_popsize(fit, R = 19, seed = 42, verbose = FALSE)

  expect_s3_class(boot, "uncounted_boot")
  expect_true(is.matrix(boot$t))
  expect_equal(nrow(boot$t), 19)
  expect_true(boot$n_converged > 0)
  expect_true(is.data.frame(boot$popsize))
  expect_true(all(c("group", "estimate", "lower", "upper") %in% names(boot$popsize)))
})

# ---- Bootstrap with constrained model ----

test_that("bootstrap works with constrained model", {
  skip_if_not_installed("fwb")
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                             gamma = 0.001, constrained = TRUE)
  boot <- bootstrap_popsize(fit, R = 19, seed = 42, verbose = FALSE)
  expect_true(boot$n_converged >= 5)
  expect_true(all(is.finite(boot$popsize$estimate)))
})

# ---- Bootstrap CI brackets point estimate ----

test_that("bootstrap CI brackets point estimate", {
  skip_if_not_installed("fwb")
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)
  boot <- bootstrap_popsize(fit, R = 49, seed = 42, verbose = FALSE)
  expect_true(all(boot$popsize$lower < boot$popsize$estimate))
  expect_true(all(boot$popsize$estimate < boot$popsize$upper))
})

# ---- Bootstrap print and summary ----

test_that("bootstrap print and summary run without error", {
  skip_if_not_installed("fwb")
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)
  boot <- bootstrap_popsize(fit, R = 19, seed = 42, verbose = FALSE)
  expect_output(print(boot))
  expect_output(summary(boot))
})

# ---- Bootstrap median differs from plugin ----

test_that("bootstrap median differs from plugin", {
  skip_if_not_installed("fwb")
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)
  boot <- bootstrap_popsize(fit, R = 49, seed = 42, verbose = FALSE)
  # Median and plugin should not be exactly identical
  expect_false(identical(boot$t0, apply(boot$t, 2, median, na.rm = TRUE)))
})

# ---- Cluster bootstrap ----

test_that("cluster bootstrap runs", {
  skip_if_not_installed("fwb")
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                             gamma = 0.001, countries = ~country)
  boot <- bootstrap_popsize(fit, R = 19, cluster = ~country,
                            seed = 42, verbose = FALSE)
  expect_true(boot$cluster)
  expect_true(boot$n_converged > 0)
})

# ---- BC CI type ----

test_that("bias-corrected bootstrap CI runs", {
  skip_if_not_installed("fwb")
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)
  boot <- bootstrap_popsize(fit, R = 49, ci_type = "bc",
                            seed = 42, verbose = FALSE)
  expect_equal(boot$ci_type, "bc")
  expect_true(all(is.finite(boot$popsize$lower)))
})
