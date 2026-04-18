# ---- LOO method tests (CRAN) ----

test_that("loo print, summary, and plot methods work on a tiny fit", {
  d <- small_data()[1:8, ]
  fit <- quick_fit(d, gamma = 0.005, countries = ~country)
  loo_res <- loo(fit, by = "obs")

  expect_output(print(loo_res), "Leave-one-out sensitivity analysis")
  expect_output(summary(loo_res), "Coefficient stability")
  expect_no_error(plot(loo_res, type = "xi"))
  expect_no_error(plot(loo_res, type = "coef"))
})

test_that("compare_loo supports labels, printing, plotting, and validation", {
  d <- small_data()[1:8, ]
  fit_1 <- quick_fit(d, gamma = 0.005, countries = ~country)
  fit_2 <- quick_fit(d, gamma = 0.005, countries = ~country, cov_alpha = ~sex)
  loo_1 <- loo(fit_1, by = "obs")
  loo_2 <- loo(fit_2, by = "obs")
  comp <- compare_loo(
    loo_1, loo_2,
    labels = c("Base", "Sex"),
    data = d,
    label_vars = c("country", "sex")
  )

  expect_s3_class(comp, "uncounted_loo_compare")
  expect_true(all(c("label", "max_abs_pct") %in% names(comp$table)))
  expect_output(print(comp), "Top")
  expect_no_error(plot(comp, type = "scatter", label_top = 2))
  expect_no_error(plot(comp, type = "bar", n = 5))

  expect_error(
    loo(quick_fit(d, gamma = 0.005), by = "country"),
    "requires 'countries'"
  )
  expect_error(
    compare_loo(loo_1, loo(fit_1, by = "country")),
    "same 'by'"
  )
})
