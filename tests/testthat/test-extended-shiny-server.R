# ---- Shiny server flows (extended, skip on CRAN) ----

skip_on_cran()

test_that("Shiny server fits models, compares them, and computes population size", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")

  shiny::testServer(uncounted:::.app_server, {
    rv$data <- small_data()
    session$setInputs(data_source = "upload")

    session$setInputs(
      col_m = "m",
      col_n = "n",
      col_N = "N",
      col_country = "country",
      col_country_fit = "country",
      model_name = "Poisson",
      method = "poisson",
      link_rho = "probit",
      estimator = "mle",
      cov_alpha = "~ sex",
      cov_beta = "~ 1",
      gamma_opt = "fixed",
      gamma_val = 0.005,
      vcov_type = "HC3",
      constrained = FALSE,
      use_cluster = FALSE,
      output_detail = "summary",
      fit_btn = 1
    )

    expect_equal(length(rv$models), 1)
    expect_true("Poisson" %in% names(rv$models))
    session$setInputs(selected_model = "Poisson")
    expect_true(nchar(output$model_output) > 0)

    session$setInputs(
      model_name = "NB",
      method = "nb",
      link_rho = "power",
      estimator = "mle",
      fit_btn = 2
    )

    expect_equal(length(rv$models), 2)
    expect_true(all(c("Poisson", "NB") %in% names(rv$models)))
    expect_match(output$compare_table, "Model")
    expect_no_error(output$compare_plot)

    session$setInputs(lr_m1 = "Poisson", lr_m2 = "NB", lr_test_btn = 1)
    expect_match(output$lr_test_result, "Likelihood ratio|Error:")

    session$setInputs(
      ps_model = "Poisson",
      ps_by = "~ sex",
      ps_bias_corr = TRUE,
      ps_total = TRUE,
      ps_compute = 1
    )

    expect_true(is.data.frame(rv$ps_result))
    expect_false(is.null(attr(rv$ps_result, "total")))
    expect_no_error(output$ps_plot)
  })
})

test_that("Shiny diagnostics branches run for fitted models", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")

  fit <- quick_fit(small_data(), gamma = 0.005, countries = ~country)
  fit_cov_gamma <- suppressWarnings(
    estimate_hidden_pop(
      data = cov_gamma_data(),
      observed = ~ m,
      auxiliary = ~ n,
      reference_pop = ~ N,
      method = "poisson",
      gamma = "estimate",
      cov_gamma = ~ 0 + z,
      countries = ~ country
    )
  )

  shiny::testServer(uncounted:::.app_server, {
    rv$models <- list(Base = fit, CovGamma = fit_cov_gamma)
    session$setInputs(data_source = "upload")

    session$setInputs(diag_model = "Base", diag_type = "residuals")
    expect_no_error(output$diag_plot)

    session$setInputs(diag_type = "rootogram", rootogram_style = "standing")
    expect_no_error(output$diag_plot)

    session$setInputs(diag_type = "gamma")
    expect_no_error(output$diag_plot)

    session$setInputs(diag_type = "explore")
    expect_no_error(output$diag_plot)

    session$setInputs(diag_type = "loo_country")
    expect_no_error(output$diag_plot)
    expect_match(output$loo_table, "Leave-one-out sensitivity analysis")

    session$setInputs(diag_model = "CovGamma", diag_type = "gamma")
    expect_no_error(output$diag_plot)
  })
})

test_that("Shiny popsize result follows its model and plots a single group", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")

  shiny::testServer(uncounted:::.app_server, {
    rv$data <- small_data()
    session$setInputs(data_source = "upload")
    session$setInputs(
      col_m = "m", col_n = "n", col_N = "N",
      col_country = "country", col_country_fit = "country",
      model_name = "A", method = "poisson", link_rho = "power", estimator = "mle",
      cov_alpha = "~ 1", cov_beta = "~ 1", gamma_opt = "fixed", gamma_val = 0.005,
      vcov_type = "HC3", constrained = FALSE, use_cluster = FALSE,
      output_detail = "summary", fit_btn = 1
    )
    expect_true("A" %in% names(rv$models))

    ## No by and no cov_alpha: one group, which used to draw nothing
    session$setInputs(ps_model = "A", ps_by = "", ps_bias_corr = TRUE,
                      ps_total = FALSE, ps_compute = 1)
    expect_equal(nrow(rv$ps_result), 1L)
    expect_no_error(output$ps_plot)

    ## Removing another model keeps the result; removing its own clears it
    session$setInputs(model_name = "B", fit_btn = 2)
    session$setInputs(selected_model = "B", remove_btn = 1)
    expect_false(is.null(rv$ps_result))
    session$setInputs(selected_model = "A", remove_btn = 2)
    expect_null(rv$ps_result)
  })
})

test_that(".ps_plot_heights uses the bias-corrected estimate when available", {
  expect_equal(.ps_plot_heights(data.frame(estimate = 2, estimate_bc = 1)), 1)
  expect_equal(.ps_plot_heights(data.frame(estimate = 2, estimate_bc = NA_real_)), 2)
  expect_equal(.ps_plot_heights(data.frame(estimate = 2)), 2)
})
