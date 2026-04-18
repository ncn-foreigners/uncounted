# ---- Lightweight bootstrap tests (CRAN) ----

test_that("bootstrap_popsize works on a tiny grouped example", {
  skip_if_not_installed("fwb")
  fit <- quick_fit(small_data(), gamma = 0.005, cov_alpha = ~sex)

  boot <- bootstrap_popsize(
    fit,
    R = 11,
    seed = 123,
    verbose = FALSE,
    by = ~ sex,
    total = TRUE,
    point_estimate = "plugin"
  )

  expect_s3_class(boot, "uncounted_boot")
  expect_equal(nrow(boot$t), 11)
  expect_true(all(c("group", "estimate", "lower", "upper") %in% names(boot$popsize)))
  expect_true("Total" %in% boot$popsize$group)
  expect_false(is.null(boot$total))
  expect_output(print(boot), "Bootstrap population size estimation")
  expect_output(summary(boot), "Bootstrap distribution summary")
})
