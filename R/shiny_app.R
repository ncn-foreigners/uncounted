utils::globalVariables("irregular_migration")

#' Launch the uncounted Shiny Application
#'
#' An interactive Shiny application for estimating unauthorised population
#' sizes using the \code{uncounted} package. The app provides a multi-tab
#' interface for data upload, model fitting, population size estimation,
#' model comparison, and diagnostics.
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#'
#' @details
#' The application has six tabs:
#' \describe{
#'   \item{Data}{Upload CSV data or use the built-in \code{irregular_migration}
#'     dataset. Select which columns map to observed counts (m), auxiliary
#'     counts (n), reference population (N), and country identifier.}
#'   \item{Model Fitting}{Specify and fit multiple models with different
#'     methods (Poisson, NB, OLS, NLS) and covariate structures. View
#'     coefficient summaries.}
#'   \item{Population Size}{Compute population size estimates stratified by
#'     user-specified grouping variables. Download results as CSV.}
#'   \item{Model Comparison}{Compare fitted models using AIC, BIC, and
#'     pseudo-R\eqn{^2} measures. Run likelihood ratio tests.}
#'   \item{Diagnostics}{Residual plots, rootograms, gamma profiles, and
#'     leave-one-out sensitivity analysis.}
#' }
#'
#' @examples
#' if (interactive()) {
#'   run_app()
#' }
#'
#' @export
run_app <- function(...) {
  for (pkg in c("shiny", "DT", "bslib")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required. Install with install.packages('",
           pkg, "')", call. = FALSE)
    }
  }

  app <- shiny::shinyApp(ui = .app_ui(), server = .app_server)
  shiny::runApp(app, ...)
}


# ── Helpers ──────────────────────────────────────────────────────────────────

## Auto-convert integer/numeric columns with <= 20 unique values to factor.
## This prevents year (2019-2024) from being treated as continuous.
.auto_factor <- function(d, exclude = NULL) {
  exclude <- unique(stats::na.omit(exclude))
  exclude <- exclude[nzchar(exclude)]
  for (col in names(d)) {
    if (col %in% exclude) next
    if (is.integer(d[[col]]) || is.numeric(d[[col]])) {
      n_unique <- length(unique(d[[col]]))
      if (n_unique <= 20 && n_unique > 1) {
        d[[col]] <- as.factor(d[[col]])
      }
    }
  }
  d
}

# Choose a comparison metric that is finite for every displayed model.
.comparison_plot_spec <- function(tab) {
  candidates <- list(
    list(metric = "AIC",
         values = tab$AIC,
         title = "AIC relative to best model",
         xlab = expression(Delta * AIC)),
    list(metric = "Deviance",
         values = tab$Deviance,
         title = "Deviance relative to best model",
         xlab = expression(Delta * Deviance)),
    list(metric = "RMSE",
         values = tab$RMSE,
         title = "RMSE relative to best model",
         xlab = "Delta RMSE")
  )

  for (spec in candidates) {
    vals <- spec$values
    if (!is.null(vals) && all(is.finite(vals))) {
      spec$relative <- vals - min(vals)
      return(spec)
    }
  }

  stop("No finite comparison metric available for plotting.", call. = FALSE)
}

.shiny_link_choices <- function(method) {
  if (is.null(method)) method <- "poisson"
  if (method %in% c("poisson", "nb", "nls")) {
    c("power", "cloglog", "logit", "probit")
  } else {
    c("power")
  }
}

.shiny_estimator_choices <- function(method) {
  if (is.null(method)) method <- "poisson"
  if (method %in% c("poisson", "nb")) {
    c("mle", "gmm", "el")
  } else {
    c("mle")
  }
}

.shiny_select_default <- function(selected, choices, default) {
  if (!is.null(selected) && selected %in% choices) {
    selected
  } else {
    default
  }
}

.formula_missing_vars <- function(formula, data_names) {
  if (is.null(formula)) return(character(0))
  setdiff(all.vars(formula), data_names)
}

# ── UI ───────────────────────────────────────────────────────────────────────

.app_ui <- function() {
  pkg_version <- tryCatch(
    as.character(utils::packageVersion("uncounted")),
    error = function(e) "1.1.0"
  )

  bslib::page_navbar(
    title = "uncounted",
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
    header = shiny::tags$style(".shiny-output-error { color: red; font-weight: bold; }"),

    # ── Tab 1: Data ──────────────────────────────────────────────────────────
    bslib::nav_panel("Data",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = "Data Source",
          shiny::radioButtons("data_source", "Source:",
            choices = c("Built-in (irregular_migration)" = "builtin",
                        "Upload CSV" = "upload")),
          shiny::conditionalPanel("input.data_source == 'upload'",
            shiny::fileInput("file_upload", "Choose CSV:", accept = ".csv")),
          shiny::hr(),
          shiny::h5("Variable Mapping"),
          shiny::uiOutput("var_observed"),
          shiny::uiOutput("var_auxiliary"),
          shiny::uiOutput("var_refpop"),
          shiny::uiOutput("var_countries")
        ),
        shiny::h4("Data Preview"),
        shiny::verbatimTextOutput("data_summary"),
        DT::dataTableOutput("data_table")
      )
    ),

    # ── Tab 2: Model Fitting ─────────────────────────────────────────────────
    bslib::nav_panel("Model Fitting",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = "Specification",
          shiny::textInput("model_name", "Model name:", value = "Model 1"),
          shiny::selectInput("method", "Method:",
            choices = c("poisson", "nb", "ols", "nls")),
          shiny::uiOutput("link_rho_ui"),
          shiny::uiOutput("estimator_ui"),
          shiny::textInput("cov_alpha", "cov_alpha formula:", value = "~ 1",
            placeholder = "~ year + sex"),
          shiny::textInput("cov_beta", "cov_beta formula:", value = "~ 1",
            placeholder = "~ year"),
          shiny::radioButtons("gamma_opt", "Gamma:",
            choices = c("Estimate" = "estimate", "Fixed" = "fixed", "None" = "none")),
          shiny::uiOutput("gamma_fixed_ui"),
          shiny::checkboxInput(
            "constrained",
            "Constrained (alpha via inverse-logit, beta via exp)",
            value = FALSE
          ),
          shiny::selectInput("vcov_type", "Robust SE type:",
            choices = c("HC3", "HC0", "HC1", "HC2", "HC4", "HC5"),
            selected = "HC3"),
          shiny::checkboxInput("use_cluster", "Cluster-robust SE", value = FALSE),
          shiny::uiOutput("var_cluster_fit"),
          shiny::uiOutput("var_countries_fit"),
          shiny::actionButton("fit_btn", "Fit Model",
            class = "btn-primary btn-lg w-100"),
          shiny::hr(),
          shiny::h5("Fitted Models"),
          shiny::uiOutput("model_list_ui"),
          shiny::actionButton("remove_btn", "Remove Selected",
            class = "btn-danger btn-sm")
        ),
        shiny::radioButtons("output_detail", NULL, inline = TRUE,
          choices = c("Summary" = "summary", "Print" = "print"),
          selected = "summary"),
        shiny::verbatimTextOutput("model_output")
      )
    ),

    # ── Tab 3: Population Size ───────────────────────────────────────────────
    bslib::nav_panel("Population Size",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = "Options",
          shiny::uiOutput("ps_model_select"),
          shiny::textInput("ps_by", "Stratify by (formula):", value = "",
            placeholder = "~ year or ~ year + sex"),
          shiny::checkboxInput("ps_bias_corr", "Bias correction", value = TRUE),
          shiny::checkboxInput("ps_total", "Include total", value = FALSE),
          shiny::actionButton("ps_compute", "Compute",
            class = "btn-primary w-100"),
          shiny::hr(),
          shiny::downloadButton("ps_download", "Download CSV")
        ),
        shiny::h4("Population Size Estimates"),
        DT::dataTableOutput("ps_table"),
        shiny::hr(),
        shiny::plotOutput("ps_plot", height = "400px")
      )
    ),

    # ── Tab 4: Model Comparison ──────────────────────────────────────────────
    bslib::nav_panel("Model Comparison",
      shiny::h4("Model Comparison Table"),
      shiny::verbatimTextOutput("compare_table"),
      shiny::hr(),
      shiny::plotOutput("compare_plot", height = "350px"),
      shiny::hr(),
      shiny::h4("Likelihood Ratio Test"),
      bslib::layout_columns(
        shiny::uiOutput("lr_model1_select"),
        shiny::uiOutput("lr_model2_select"),
        shiny::actionButton("lr_test_btn", "Run LR Test",
          class = "btn-primary")
      ),
      shiny::verbatimTextOutput("lr_test_result")
    ),

    # ── Tab 5: Diagnostics ───────────────────────────────────────────────────
    bslib::nav_panel("Diagnostics",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = "Options",
          shiny::uiOutput("diag_model_select"),
          shiny::selectInput("diag_type", "Diagnostic:",
            choices = c("Residual plots" = "residuals",
                        "Rootogram" = "rootogram",
                        "Gamma profile" = "gamma",
                        "LOO (by country)" = "loo_country",
                        "Exploratory (log-log)" = "explore")),
          shiny::conditionalPanel("input.diag_type == 'rootogram'",
            shiny::selectInput("rootogram_style", "Style:",
              choices = c("hanging", "suspended", "standing")))
        ),
        shiny::plotOutput("diag_plot", height = "550px"),
        shiny::hr(),
        shiny::verbatimTextOutput("loo_table")
      )
    ),

    # ── Tab 6: About ───────────────────────────────────────────────────────
    bslib::nav_panel("About",
      shiny::div(class = "container", style = "max-width: 800px; margin-top: 20px;",
        shiny::h2("uncounted: Estimating Unauthorised Population Size"),
        shiny::p("This application implements the community-anchored modelling framework
          for estimating unauthorised population sizes from aggregated administrative data."),
        shiny::h4("Model"),
        shiny::p(shiny::HTML("The core model assumes:<br>
          <code>E(m<sub>i</sub>) = N<sub>i</sub><sup>&alpha;</sup> (&gamma; + n<sub>i</sub>/N<sub>i</sub>)<sup>&beta;</sup></code><br>
          where <b>m</b> = observed unauthorised counts, <b>N</b> = reference (authorised) population,
          <b>n</b> = auxiliary detection covariate, <b>&alpha;</b> = community anchor elasticity,
          <b>&beta;</b> = exposure elasticity, <b>&gamma;</b> = baseline rate offset.")),
        shiny::h4("Estimation Methods"),
        shiny::tags$ul(
          shiny::tags$li(shiny::strong("Poisson PML"), " (recommended): consistent under correct
            mean specification, HC-robust standard errors"),
          shiny::tags$li(shiny::strong("Negative Binomial"), ": accommodates overdispersion via
            dispersion parameter theta"),
          shiny::tags$li(shiny::strong("OLS"), ": log-linearised (biased with many zeros)"),
          shiny::tags$li(shiny::strong("NLS"), ": nonlinear least squares")
        ),
        shiny::h4("Population Size"),
        shiny::p(shiny::HTML("The estimated population size is: <code>&xi; = &Sigma; N<sub>i</sub><sup>&alpha;</sup></code>,
          with delta-method CIs and optional bias correction.")),
        shiny::h4("Workflow"),
        shiny::tags$ol(
          shiny::tags$li("Upload data or use the built-in dataset"),
          shiny::tags$li("Map variables (observed, auxiliary, reference population)"),
          shiny::tags$li("Fit one or more models with different specifications"),
          shiny::tags$li("Compare models (AIC, BIC, LR test)"),
          shiny::tags$li("Estimate population size with stratification and CIs"),
          shiny::tags$li("Assess quality with diagnostics (residuals, LOO, rootograms)")
        ),
        shiny::hr(),
        shiny::p(shiny::HTML(sprintf(
          paste0("Package: <code>uncounted</code> v%s | ",
                 "Author: Maciej Ber&#281;sewicz | ",
                 "<a href='https://github.com/ncn-foreigners/paper-unauthorized-pop' target='_blank'>GitHub</a>"),
          pkg_version
        )))
      )
    )
  )
}


# ── Server ───────────────────────────────────────────────────────────────────

.app_server <- function(input, output, session) {

  rv <- shiny::reactiveValues(
    data = NULL,
    models = list(),
    selected_model = NULL,
    ps_result = NULL
  )

  # ── Data loading ─────────────────────────────────────────────────────────
  shiny::observe({
    if (input$data_source == "builtin") {
      rv$data <- .auto_factor(irregular_migration)
    } else if (!is.null(input$file_upload)) {
      d <- tryCatch(
        utils::read.csv(input$file_upload$datapath, stringsAsFactors = TRUE),
        error = function(e) { shiny::showNotification(e$message, type = "error"); NULL }
      )
      rv$data <- if (!is.null(d)) .auto_factor(d) else NULL
    }
  })

  cols <- shiny::reactive({
    if (is.null(rv$data)) return(character(0))
    names(rv$data)
  })

  processed_data <- shiny::reactive({
    shiny::req(rv$data)
    .auto_factor(rv$data, exclude = c(input$col_m, input$col_n, input$col_N))
  })

  output$var_observed <- shiny::renderUI({
    shiny::selectInput("col_m", "Observed (m):", choices = cols(),
                       selected = if ("m" %in% cols()) "m" else cols()[1])
  })
  output$var_auxiliary <- shiny::renderUI({
    shiny::selectInput("col_n", "Auxiliary (n):", choices = cols(),
                       selected = if ("n" %in% cols()) "n" else cols()[1])
  })
  output$var_refpop <- shiny::renderUI({
    shiny::selectInput("col_N", "Reference pop (N):", choices = cols(),
                       selected = if ("N" %in% cols()) "N" else cols()[1])
  })
  output$var_countries <- shiny::renderUI({
    ch <- c("(none)" = "", cols())
    sel <- if ("country_code" %in% cols()) "country_code" else ""
    shiny::selectInput("col_country", "Countries (cluster):", choices = ch,
                       selected = sel)
  })

  output$gamma_fixed_ui <- shiny::renderUI({
    shiny::req(input$gamma_opt == "fixed")
    shiny::numericInput("gamma_val", "Gamma value:", value = 0.01,
                        min = 0, step = 0.001)
  })

  output$link_rho_ui <- shiny::renderUI({
    choices <- .shiny_link_choices(input$method)
    input_link <- if (identical(input$link_rho, "logistic")) "logit" else input$link_rho
    selected <- .shiny_select_default(input_link, choices, "power")
    shiny::selectInput("link_rho", "Detection link:", choices = choices,
                       selected = selected)
  })

  output$estimator_ui <- shiny::renderUI({
    choices <- .shiny_estimator_choices(input$method)
    selected <- .shiny_select_default(input$estimator, choices, "mle")
    shiny::selectInput("estimator", "Estimator:", choices = choices,
                       selected = selected)
  })

  output$var_countries_fit <- shiny::renderUI({
    ## Countries variable for grouping (popsize, LOO)
    exclude <- c(input$col_m, input$col_n, input$col_N)
    ch <- c("(none)" = "", setdiff(cols(), exclude))
    sel <- if (!is.null(input$col_country) && input$col_country %in% ch) input$col_country else ""
    shiny::selectInput("col_country_fit", "Countries (grouping):", choices = ch,
                       selected = sel)
  })

  output$var_cluster_fit <- shiny::renderUI({
    shiny::req(isTRUE(input$use_cluster))
    ## Cluster variable for cluster-robust SE
    exclude <- c(input$col_m, input$col_n, input$col_N)
    ch <- setdiff(cols(), exclude)
    sel <- if (!is.null(input$col_country) && input$col_country %in% ch) input$col_country else ch[1]
    shiny::selectInput("col_cluster", "Cluster variable:", choices = ch,
                       selected = sel)
  })

  output$data_summary <- shiny::renderPrint({
    shiny::req(processed_data())
    d <- processed_data()
    cat("Rows:", nrow(d), "| Columns:", ncol(d), "\n")
    cat("Column names:", paste(names(d), collapse = ", "), "\n")
    if (!is.null(input$col_m) && input$col_m %in% names(d)) {
      cat("Zeros in m:", sum(d[[input$col_m]] == 0), "/", nrow(d), "\n")
    }
    if (!is.null(input$col_n) && input$col_n %in% names(d)) {
      cat("Zeros in n:", sum(d[[input$col_n]] == 0), "/", nrow(d), "\n")
    }
  })

  output$data_table <- DT::renderDataTable({
    shiny::req(processed_data())
    DT::datatable(utils::head(processed_data(), 200), options = list(pageLength = 10,
      scrollX = TRUE))
  })

  # ── Model fitting ───────────────────────────────────────────────────────
  shiny::observeEvent(input$fit_btn, {
    shiny::req(processed_data(), input$col_m, input$col_n, input$col_N)

    gamma_arg <- switch(input$gamma_opt,
      estimate = "estimate",
      fixed = if (!is.null(input$gamma_val)) input$gamma_val else 0.01,
      none = NULL
    )
    countries_arg <- if (!is.null(input$col_country_fit) && input$col_country_fit != "") {
      stats::as.formula(paste("~", input$col_country_fit))
    } else NULL

    cov_a <- tryCatch(stats::as.formula(input$cov_alpha), error = function(e) NULL)
    cov_b <- tryCatch(stats::as.formula(input$cov_beta), error = function(e) NULL)

    if (is.null(cov_a)) {
      shiny::showNotification("Invalid cov_alpha formula", type = "error"); return()
    }
    if (is.null(cov_b)) {
      shiny::showNotification("Invalid cov_beta formula", type = "error"); return()
    }

    ## Validate that formula variables exist in data
    check_vars <- function(f, label) {
      if (is.null(f)) return(TRUE)
      missing <- .formula_missing_vars(f, names(processed_data()))
      if (length(missing) > 0) {
        shiny::showNotification(
          paste0(label, ": variable(s) not in data: ", paste(missing, collapse = ", ")),
          type = "error")
        return(FALSE)
      }
      TRUE
    }
    if (!check_vars(cov_a, "cov_alpha")) return()
    if (!check_vars(cov_b, "cov_beta")) return()

    ## Resolve all args to plain R objects (no Shiny input$ references)
    obs_f <- stats::as.formula(paste("~", input$col_m))
    aux_f <- stats::as.formula(paste("~", input$col_n))
    ref_f <- stats::as.formula(paste("~", input$col_N))
    alpha_f <- if (input$cov_alpha == "~ 1") NULL else cov_a
    beta_f  <- if (input$cov_beta == "~ 1") NULL else cov_b
    method_val <- input$method
    link_rho_val <- input$link_rho
    estimator_val <- input$estimator
    vcov_val <- input$vcov_type
    constr_val <- isTRUE(input$constrained)
    dat <- processed_data()

    ## Cluster-robust SE
    cluster_arg <- NULL
    if (isTRUE(input$use_cluster) && !is.null(input$col_cluster) &&
        input$col_cluster != "") {
      cluster_arg <- stats::as.formula(paste("~", input$col_cluster))
    }

    shiny::withProgress(message = "Fitting model...", {
      fit <- tryCatch(
        estimate_hidden_pop(
          data = dat, observed = obs_f, auxiliary = aux_f,
          reference_pop = ref_f, method = method_val,
          cov_alpha = alpha_f, cov_beta = beta_f,
          gamma = gamma_arg, link_rho = link_rho_val,
          estimator = estimator_val, vcov = vcov_val,
          constrained = constr_val, countries = countries_arg,
          cluster = cluster_arg
        ),
        warning = function(w) {
          shiny::showNotification(paste("Warning:", w$message), type = "warning")
          suppressWarnings(estimate_hidden_pop(
            data = dat, observed = obs_f, auxiliary = aux_f,
            reference_pop = ref_f, method = method_val,
            cov_alpha = alpha_f, cov_beta = beta_f,
            gamma = gamma_arg, link_rho = link_rho_val,
            estimator = estimator_val, vcov = vcov_val,
            constrained = constr_val, countries = countries_arg,
            cluster = cluster_arg
          ))
        },
        error = function(e) {
          shiny::showNotification(paste("Error:", e$message), type = "error")
          NULL
        }
      )
    })

    if (!is.null(fit)) {
      nm <- input$model_name
      if (nm %in% names(rv$models)) nm <- paste0(nm, " (", length(rv$models) + 1, ")")
      rv$models[[nm]] <- fit
      rv$selected_model <- nm
      shiny::showNotification(paste("Model", nm, "fitted successfully"), type = "message")
    }
  })

  model_names <- shiny::reactive(names(rv$models))

  output$model_list_ui <- shiny::renderUI({
    shiny::req(length(rv$models) > 0)
    shiny::selectInput("selected_model", "Select model:", choices = model_names(),
                       selected = rv$selected_model)
  })

  shiny::observeEvent(input$remove_btn, {
    shiny::req(input$selected_model)
    rv$models[[input$selected_model]] <- NULL
    if (length(rv$models) > 0) rv$selected_model <- names(rv$models)[1]
  })

  current_fit <- shiny::reactive({
    shiny::req(input$selected_model, rv$models[[input$selected_model]])
    rv$models[[input$selected_model]]
  })

  output$model_output <- shiny::renderPrint({
    shiny::validate(shiny::need(length(rv$models) > 0, "Fit a model first"))
    shiny::req(current_fit())
    if (input$output_detail == "summary") summary(current_fit()) else print(current_fit())
  })

  # ── Population size ────────────────────────────────────────────────────
  output$ps_model_select <- shiny::renderUI({
    shiny::req(length(rv$models) > 0)
    shiny::selectInput("ps_model", "Model:", choices = model_names())
  })

  shiny::observeEvent(input$ps_compute, {
    shiny::req(input$ps_model, rv$models[[input$ps_model]])

    by_arg <- if (nchar(trimws(input$ps_by)) > 0) {
      tryCatch(stats::as.formula(input$ps_by), error = function(e) {
        shiny::showNotification("Invalid 'by' formula", type = "error"); NULL
      })
    } else NULL

    rv$ps_result <- tryCatch(
      popsize(rv$models[[input$ps_model]], by = by_arg,
              bias_correction = input$ps_bias_corr, total = input$ps_total),
      error = function(e) {
        shiny::showNotification(paste("Error:", e$message), type = "error"); NULL
      }
    )
  })

  output$ps_table <- DT::renderDataTable({
    shiny::req(rv$ps_result)
    DT::datatable(rv$ps_result, options = list(pageLength = 20))
  })

  output$ps_plot <- shiny::renderPlot({
    shiny::req(rv$ps_result)
    ps <- rv$ps_result
    if (nrow(ps) > 1) {
      n_groups <- nrow(ps)
      par(mar = c(7, 4, 3, 1))
      bp <- barplot(ps$estimate, names.arg = ps$group, las = 2,
                    ylab = expression(hat(xi)), main = "Population Size Estimates",
                    col = "steelblue", border = NA)
      if (!is.null(ps$lower) && !is.null(ps$upper)) {
        arrows(bp, ps$lower, bp, ps$upper, angle = 90, code = 3, length = 0.05)
      }
    }
  })

  output$ps_download <- shiny::downloadHandler(
    filename = function() paste0("popsize_", Sys.Date(), ".csv"),
    content = function(file) {
      shiny::req(rv$ps_result)
      utils::write.csv(rv$ps_result, file, row.names = FALSE)
    }
  )

  # ── Model comparison ───────────────────────────────────────────────────
  output$compare_table <- shiny::renderPrint({
    shiny::validate(shiny::need(length(rv$models) >= 2, "Fit at least 2 models to compare"))
    do.call(compare_models, rv$models)
  })

  output$compare_plot <- shiny::renderPlot({
    shiny::validate(shiny::need(length(rv$models) >= 2, "Fit at least 2 models to compare"))
    comp <- do.call(compare_models, rv$models)
    tab <- comp$table
    spec <- tryCatch(
      .comparison_plot_spec(tab),
      error = function(e) {
        plot.new()
        text(0.5, 0.5, e$message, cex = 1.1)
        NULL
      }
    )
    if (is.null(spec)) return(invisible())
    ord <- order(spec$relative)
    tab <- tab[ord, , drop = FALSE]
    rel_vals <- spec$relative[ord]
    par(mar = c(5, 10, 3, 2))
    barplot(rel_vals, names.arg = tab$Model, horiz = TRUE,
            las = 1, xlab = spec$xlab,
            main = spec$title,
            col = ifelse(rel_vals == 0, "steelblue", "grey70"),
            border = NA)
  })

  output$lr_model1_select <- shiny::renderUI({
    shiny::req(length(rv$models) >= 2)
    shiny::selectInput("lr_m1", "Model 1:", choices = model_names())
  })
  output$lr_model2_select <- shiny::renderUI({
    shiny::req(length(rv$models) >= 2)
    shiny::selectInput("lr_m2", "Model 2:", choices = model_names())
  })

  output$lr_test_result <- shiny::renderPrint({
    shiny::req(input$lr_test_btn)
    shiny::isolate({
      shiny::req(input$lr_m1, input$lr_m2, input$lr_m1 != input$lr_m2)
      tryCatch(
        lrtest(rv$models[[input$lr_m1]], rv$models[[input$lr_m2]]),
        error = function(e) cat("Error:", e$message, "\n")
      )
    })
  })

  # ── Diagnostics ────────────────────────────────────────────────────────
  output$diag_model_select <- shiny::renderUI({
    shiny::req(length(rv$models) > 0)
    shiny::selectInput("diag_model", "Model:", choices = model_names())
  })

  output$diag_plot <- shiny::renderPlot({
    shiny::validate(shiny::need(length(rv$models) > 0, "Fit a model first (Model Fitting tab)"))
    shiny::req(input$diag_model)
    shiny::validate(shiny::need(input$diag_model %in% names(rv$models), "Select a fitted model"))
    fit <- rv$models[[input$diag_model]]

    switch(input$diag_type,
      residuals = {
        par(mfrow = c(2, 2))
        plot(fit, ask = FALSE)
      },
      rootogram = {
        rootogram(fit, style = input$rootogram_style, max_count = 50)
      },
      gamma = {
        if (isTRUE(fit$has_cov_gamma)) {
          plot.new()
          text(0.5, 0.5,
               "Gamma profile not available\nfor covariate-varying gamma models",
               cex = 1.2)
        } else if (!is.null(fit$gamma)) {
          profile_gamma(fit, gamma_grid = seq(1e-4, 0.3, length.out = 30))
        } else {
          plot.new()
          text(0.5, 0.5, "Gamma not estimated for this model", cex = 1.5)
        }
      },
      loo_country = {
        shiny::withProgress(message = "Running LOO...", {
          loo_res <- loo(fit, by = "country", verbose = FALSE)
        })
        plot(loo_res, type = "xi", n = 20)
      },
      explore = {
        par(mfrow = c(1, 2))
        plot_explore(fit)
      }
    )
  })

  output$loo_table <- shiny::renderPrint({
    shiny::req(input$diag_type == "loo_country")
    shiny::validate(shiny::need(length(rv$models) > 0, ""))
    shiny::req(input$diag_model, input$diag_model %in% names(rv$models))
    fit <- rv$models[[input$diag_model]]
    loo_res <- loo(fit, by = "country", verbose = FALSE)
    print(loo_res, n = 15)
  })
}
