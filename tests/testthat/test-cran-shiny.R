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

test_that("shiny helper functions return the supported choices", {
  expect_identical(
    uncounted:::.shiny_link_choices("poisson"),
    c("power", "cloglog", "logit", "probit")
  )
  expect_identical(uncounted:::.shiny_link_choices("ols"), "power")
  expect_identical(
    uncounted:::.shiny_estimator_choices("nb"),
    c("mle", "gmm", "el")
  )
  expect_identical(uncounted:::.shiny_estimator_choices("nls"), "mle")
  expect_identical(
    uncounted:::.shiny_select_default("logit", c("power", "logit"), "power"),
    "logit"
  )
  expect_identical(
    uncounted:::.shiny_select_default("bad", c("power", "logit"), "power"),
    "power"
  )
  expect_identical(
    uncounted:::.formula_missing_vars(~ year + sex + missing_var,
                                      c("year", "sex", "m", "n", "N")),
    "missing_var"
  )
})

test_that(".auto_factor leaves degenerate numeric columns numeric", {
  d <- data.frame(
    value = c(1, 1, 1),
    binary = c(0, 1, 0)
  )

  out <- uncounted:::.auto_factor(d)
  expect_false(is.factor(out$value))
  expect_true(is.factor(out$binary))
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

test_that("comparison plot helper errors when no metric is finite", {
  tab <- data.frame(
    Model = c("A", "B"),
    AIC = c(NA_real_, NA_real_),
    Deviance = c(NA_real_, NA_real_),
    RMSE = c(NA_real_, NA_real_)
  )

  expect_error(
    uncounted:::.comparison_plot_spec(tab),
    "No finite comparison metric"
  )
})

test_that("Shiny server loads uploaded data and renders a summary", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")

  csv_path <- tempfile(fileext = ".csv")
  utils::write.csv(small_data(), csv_path, row.names = FALSE)

  shiny::testServer(uncounted:::.app_server, {
    session$setInputs(
      data_source = "upload",
      file_upload = data.frame(
        name = basename(csv_path),
        size = file.info(csv_path)$size,
        type = "text/csv",
        datapath = csv_path
      ),
      col_m = "m",
      col_n = "n",
      col_N = "N",
      col_country = "country"
    )

    expect_equal(nrow(rv$data), nrow(small_data()))
    expect_match(output$data_summary, "Rows:")
    expect_match(output$data_summary, "Zeros in m:")
  })
})

test_that("uploaded m, n and N stay numeric with few distinct values", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")

  ## A country panel where N is constant per country: every column has at most
  ## 20 distinct values, so auto-factoring at load made the data unfittable
  d <- expand.grid(country = paste0("C", 1:6), year = 2019:2021)
  d$N <- rep(c(500, 1000, 2000, 4000, 8000, 16000), times = 3)
  d$n <- rep(c(3, 5, 8, 12, 20, 30), times = 3)
  d$m <- rep(c(1, 2, 4, 6, 9, 14), times = 3)
  csv_path <- tempfile(fileext = ".csv")
  utils::write.csv(d, csv_path, row.names = FALSE)

  shiny::testServer(uncounted:::.app_server, {
    session$setInputs(
      data_source = "upload",
      file_upload = data.frame(
        name = basename(csv_path),
        size = file.info(csv_path)$size,
        type = "text/csv",
        datapath = csv_path
      ),
      col_m = "m",
      col_n = "n",
      col_N = "N",
      col_country = "country"
    )
    dat <- processed_data()
    expect_true(is.numeric(dat$m) && is.numeric(dat$n) && is.numeric(dat$N))
    ## Other small-range numeric columns are still treated as factors
    expect_true(is.factor(dat$year))
    expect_s3_class(estimate_hidden_pop(dat, ~m, ~n, ~N, gamma = 0.005), "uncounted")
  })
})
