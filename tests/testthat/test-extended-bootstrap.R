# ---- Bootstrap tests (extended, skip on CRAN) ----

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
