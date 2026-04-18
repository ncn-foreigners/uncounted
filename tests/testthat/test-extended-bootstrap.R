# ---- Bootstrap tests (extended, skip on CRAN) ----

skip_on_cran()

test_that("bootstrap R=199 on full data", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  fit <- quick_fit(testdata, gamma = 0.005)
  boot <- bootstrap_popsize(fit, R = 199, seed = 42, verbose = FALSE)

  expect_s3_class(boot, "uncounted_boot")
  expect_true(is.matrix(boot$t))
  expect_equal(nrow(boot$t), 199)
  expect_true(boot$n_converged > 199 * 0.75)
  expect_true(is.data.frame(boot$popsize))
  expect_true(all(c("group", "estimate", "lower", "upper") %in% names(boot$popsize)))
})

test_that("bootstrap CI brackets point estimate", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  fit <- quick_fit(testdata, gamma = 0.005)
  boot <- bootstrap_popsize(fit, R = 199, seed = 42, verbose = FALSE)
  expect_true(all(boot$popsize$lower < boot$popsize$estimate))
  expect_true(all(boot$popsize$estimate < boot$popsize$upper))
})

test_that("bootstrap with constrained model", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  fit <- quick_fit(testdata, gamma = 0.005, constrained = TRUE)
  boot <- bootstrap_popsize(fit, R = 99, seed = 42, verbose = FALSE)
  expect_true(boot$n_converged >= 50)
  expect_true(all(is.finite(boot$popsize$estimate)))
})

test_that("bootstrap with NB model", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  fit <- quick_fit(testdata, method = "nb", gamma = 0.005)
  boot <- bootstrap_popsize(fit, R = 99, seed = 42, verbose = FALSE)
  expect_true(boot$n_converged >= 50)
})

test_that("cluster bootstrap R=199", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  fit <- quick_fit(testdata, gamma = 0.005, countries = ~country)
  boot <- bootstrap_popsize(fit, R = 199, cluster = ~country,
                            seed = 42, verbose = FALSE)
  expect_true(boot$cluster)
  expect_true(boot$n_converged > 0)
})

test_that("bias-corrected bootstrap CI", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  fit <- quick_fit(testdata, gamma = 0.005)
  boot <- bootstrap_popsize(fit, R = 199, ci_type = "bc",
                            seed = 42, verbose = FALSE)
  expect_equal(boot$ci_type, "bc")
  expect_true(all(is.finite(boot$popsize$lower)))
})

test_that("bootstrap print and summary run without error", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  fit <- quick_fit(testdata, gamma = 0.005)
  boot <- bootstrap_popsize(fit, R = 49, seed = 42, verbose = FALSE)
  expect_output(print(boot))
  expect_output(summary(boot))
})

test_that("bootstrap median differs from plugin", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  fit <- quick_fit(testdata, gamma = 0.005)
  boot <- bootstrap_popsize(fit, R = 99, seed = 42, verbose = FALSE)
  expect_false(identical(boot$t0, apply(boot$t, 2, median, na.rm = TRUE)))
})

# ---- Bootstrap coverage: point_estimate variants ----

test_that("bootstrap point_estimate='plugin' works", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  fit <- quick_fit(testdata, gamma = 0.005)
  boot <- bootstrap_popsize(fit, R = 49, seed = 42, verbose = FALSE,
                            point_estimate = "plugin")
  expect_true(all(is.finite(boot$popsize$estimate)))
})

test_that("bootstrap point_estimate='mean' works", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  fit <- quick_fit(testdata, gamma = 0.005)
  boot <- bootstrap_popsize(fit, R = 49, seed = 42, verbose = FALSE,
                            point_estimate = "mean")
  expect_true(all(is.finite(boot$popsize$estimate)))
})

# ---- Bootstrap coverage: total ----

test_that("bootstrap with total=TRUE includes total", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  fit <- quick_fit(testdata, gamma = 0.005, cov_alpha = ~sex)
  boot <- bootstrap_popsize(fit, R = 49, seed = 42, verbose = FALSE,
                            total = TRUE)
  expect_true(!is.null(boot$total))
  expect_true(boot$total$plugin > 0)
  expect_true(boot$total$lower < boot$total$upper)
  # Total row should appear in popsize and popsize_full data frames
  expect_true("Total" %in% boot$popsize$group)
  expect_true("Total" %in% boot$popsize_full$group)
})

test_that("bootstrap total=TRUE uses correct full-gradient BC", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  fit <- quick_fit(testdata, gamma = 0.005, cov_alpha = ~sex)
  boot <- bootstrap_popsize(fit, R = 49, seed = 42, verbose = FALSE,
                            total = TRUE)
  # Compare to popsize(total=TRUE) for correct BC
  ps0 <- popsize(fit, bias_correction = TRUE, total = TRUE)
  expect_equal(boot$total$plugin_bc, attr(ps0, "total")$estimate_bc)
})

# ---- Bootstrap: boot_params ----

test_that("bootstrap exposes boot_params matrix", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  fit <- quick_fit(testdata, gamma = 0.005)
  boot <- bootstrap_popsize(fit, R = 29, seed = 42, verbose = FALSE)
  expect_true(!is.null(boot$boot_params))
  expect_true(is.matrix(boot$boot_params))
  expect_equal(nrow(boot$boot_params), 29)
  expect_true("alpha" %in% colnames(boot$boot_params))
  expect_true("beta" %in% colnames(boot$boot_params))
})

test_that("bootstrap boot_params includes gamma for estimated gamma", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  d <- positive_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                             gamma = "estimate")
  boot <- bootstrap_popsize(fit, R = 29, seed = 42, verbose = FALSE)
  expect_true("gamma" %in% colnames(boot$boot_params))
})

test_that("bootstrap boot_params includes theta for NB", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  fit <- quick_fit(testdata, method = "nb", gamma = 0.005)
  boot <- bootstrap_popsize(fit, R = 29, seed = 42, verbose = FALSE)
  expect_true("theta" %in% colnames(boot$boot_params))
})

# ---- Bootstrap coverage: by parameter ----

test_that("bootstrap with by parameter works", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  fit <- quick_fit(testdata, gamma = 0.005)
  boot <- bootstrap_popsize(fit, R = 29, seed = 42, verbose = FALSE,
                            by = ~ sex)
  expect_true(nrow(boot$popsize) >= 2)
  expect_true(any(c("M", "F") %in% boot$popsize$group))
})
