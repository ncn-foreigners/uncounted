# ---- popsize plot and compare_popsize tests ----

test_that("popsize returns uncounted_popsize class", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  ps <- popsize(fit)
  expect_s3_class(ps, "uncounted_popsize")
  expect_s3_class(ps, "data.frame")
})

test_that("plot.uncounted_popsize runs without error", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  ps <- popsize(fit)
  expect_no_error(plot(ps))
})

test_that("plot.uncounted_popsize type='compare' runs", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  ps <- popsize(fit, bias_correction = TRUE)
  expect_no_error(plot(ps, type = "compare"))
})

test_that("plot popsize with by variable", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ps <- popsize(fit, by = ~sex)
  expect_no_error(plot(ps))
})

# ---- compare_popsize ----

test_that("compare_popsize works with two models", {
  d <- small_data()
  fit_po <- quick_fit(d, gamma = 0.005)
  fit_io <- quick_fit(d, method = "iols", gamma = 0.005)
  comp <- compare_popsize(fit_po, fit_io, labels = c("Poisson", "iOLS"))
  expect_s3_class(comp, "uncounted_popsize_compare")
  expect_true("model" %in% names(comp$table))
})

test_that("compare_popsize with by variable", {
  d <- small_data()
  fit_po <- quick_fit(d, gamma = 0.005)
  fit_io <- quick_fit(d, method = "iols", gamma = 0.005)
  comp <- compare_popsize(fit_po, fit_io, by = ~sex,
                          labels = c("Poisson", "iOLS"))
  expect_true(nrow(comp$table) >= 4)
})

test_that("compare_popsize plot runs", {
  d <- small_data()
  fit_po <- quick_fit(d, gamma = 0.005)
  fit_io <- quick_fit(d, method = "iols", gamma = 0.005)
  comp <- compare_popsize(fit_po, fit_io)
  expect_no_error(plot(comp))
})

test_that("compare_popsize print works", {
  d <- small_data()
  fit_po <- quick_fit(d, gamma = 0.005)
  fit_io <- quick_fit(d, method = "iols", gamma = 0.005)
  comp <- compare_popsize(fit_po, fit_io)
  expect_output(print(comp))
})

test_that("compare_popsize with three models", {
  d <- small_data()
  fit_po <- quick_fit(d, gamma = 0.005)
  suppressWarnings({
    fit_nb <- quick_fit(d, method = "nb", gamma = 0.005)
  })
  fit_io <- quick_fit(d, method = "iols", gamma = 0.005)
  comp <- compare_popsize(fit_po, fit_nb, fit_io,
                          labels = c("Poisson", "NB", "iOLS"))
  expect_equal(comp$n_models, 3)
  expect_equal(length(unique(comp$table$model)), 3)
})
