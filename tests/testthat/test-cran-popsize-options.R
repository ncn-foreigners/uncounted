# ---- popsize() options coverage ----

# ---- bias_correction = FALSE ----

test_that("bias_correction=FALSE: estimate_bc is NA", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ps <- popsize(fit, bias_correction = FALSE)
  expect_true(all(is.na(ps$estimate_bc)))
})

test_that("bias_correction=TRUE: estimate_bc differs from estimate", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ps <- popsize(fit, bias_correction = TRUE)
  expect_false(all(ps$estimate == ps$estimate_bc))
})

# ---- total = TRUE ----

test_that("total=TRUE adds total attribute", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  ps <- popsize(fit, total = TRUE)
  tot <- attr(ps, "total")
  expect_false(is.null(tot))
  expect_true(tot$estimate > 0)
  expect_true(tot$lower < tot$estimate)
  expect_true(tot$estimate < tot$upper)
})

test_that("total=FALSE has no total attribute", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  ps <- popsize(fit, total = FALSE)
  expect_null(attr(ps, "total"))
})

test_that("total=TRUE prints total in output", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  ps <- popsize(fit, total = TRUE)
  out <- capture.output(print(ps))
  expect_true(any(grepl("Total", out)))
})

test_that("total=TRUE with bias_correction=FALSE", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  ps <- popsize(fit, total = TRUE, bias_correction = FALSE)
  tot <- attr(ps, "total")
  expect_false(is.null(tot))
  expect_true(is.na(tot$estimate_bc))
})

# ---- by parameter ----

test_that("by=~sex gives 2 rows", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ps <- popsize(fit, by = ~sex)
  expect_equal(nrow(ps), 2)
  expect_true(all(c("F", "M") %in% ps$group))
})

test_that("by=~country gives one row per country", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ps <- popsize(fit, by = ~country)
  expect_equal(nrow(ps), length(unique(d$country)))
})

test_that("by + total + bias_correction all work together", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ps <- popsize(fit, by = ~sex, total = TRUE, bias_correction = TRUE)
  expect_equal(nrow(ps), 2)
  expect_false(is.null(attr(ps, "total")))
  expect_true(all(ps$estimate_bc > 0))
  expect_true(all(ps$estimate_bc <= ps$estimate))
})

# ---- level parameter ----

test_that("level=0.90 gives narrower CI than level=0.95", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ps_90 <- popsize(fit, level = 0.90)
  ps_95 <- popsize(fit, level = 0.95)
  width_90 <- ps_90$upper - ps_90$lower
  width_95 <- ps_95$upper - ps_95$lower
  expect_true(all(width_90 < width_95))
})

# ---- all methods ----

test_that("popsize works for all methods with all options", {
  d <- small_data()
  for (method in c("ols", "poisson", "nb", "iols")) {
    suppressWarnings({
      fit <- quick_fit(d, method = method, gamma = 0.005, cov_alpha = ~sex)
    })
    ps <- popsize(fit, by = ~sex, total = TRUE, bias_correction = TRUE,
                  level = 0.90)
    expect_equal(nrow(ps), 2, info = method)
    expect_false(is.null(attr(ps, "total")), info = method)
    expect_true(all(ps$estimate > 0), info = method)
    expect_true(all(ps$estimate_bc > 0), info = method)
  }
})

# ---- print method ----

test_that("print.uncounted_popsize works without total", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  ps <- popsize(fit)
  expect_output(print(ps))
})

test_that("print.uncounted_popsize shows total when present", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~sex)
  ps <- popsize(fit, total = TRUE)
  out <- capture.output(print(ps))
  expect_true(any(grepl("Total", out)))
})
