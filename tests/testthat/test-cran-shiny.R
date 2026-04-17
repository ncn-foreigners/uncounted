# ---- Shiny app tests (CRAN) ----

test_that("run_app is exported and callable", {
  expect_true("run_app" %in% getNamespaceExports("uncounted"))
  expect_true(is.function(uncounted::run_app))
})

test_that("run_app checks for required packages", {
  expect_true("..." %in% names(formals(run_app)))
})

test_that("app UI function returns valid shiny tag", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")

  ui <- uncounted:::.app_ui()
  expect_true(inherits(ui, "shiny.tag") || inherits(ui, "shiny.tag.list"))
})

test_that("app server function has correct signature", {
  expect_true(is.function(uncounted:::.app_server))
  args <- names(formals(uncounted:::.app_server))
  expect_true(all(c("input", "output", "session") %in% args))
})

test_that(".auto_factor preserves mapped measure columns", {
  d <- data.frame(
    m = c(0, 1, 2),
    n = c(0, 1, 2),
    N = c(10, 10, 20),
    year = c(2019, 2020, 2021)
  )

  out <- uncounted:::.auto_factor(d, exclude = c("m", "n", "N"))

  expect_false(is.factor(out$m))
  expect_false(is.factor(out$n))
  expect_false(is.factor(out$N))
  expect_true(is.factor(out$year))
})

test_that("comparison plot helper falls back when AIC is unavailable", {
  tab <- data.frame(
    Model = c("MLE", "GMM"),
    AIC = c(12, NA_real_),
    Deviance = c(5, 8),
    RMSE = c(1.2, 1.8)
  )

  spec <- uncounted:::.comparison_plot_spec(tab)

  expect_identical(spec$metric, "Deviance")
  expect_equal(spec$relative, c(0, 3))
})
