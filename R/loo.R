#' Leave-One-Out Sensitivity Analysis
#'
#' Refits the model dropping one observation (or one country/group) at a time.
#' Returns the impact on coefficient estimates, fitted values, and
#' \eqn{\xi} (total population size estimate).
#'
#' @details
#' **Purpose.** LOO sensitivity analysis assesses how much each observation
#' (or country) influences the estimated population size. An observation
#' with a large \eqn{|\Delta\xi|} relative to the full-model \eqn{\hat{\xi}}
#' is "influential" and may warrant closer inspection.
#'
#' **\code{by = "obs"} vs \code{by = "country"}.** With \code{by = "obs"},
#' each row is dropped individually (n refits). This identifies individual
#' data points that drive the estimate, useful for detecting outliers or
#' data errors. With \code{by = "country"}, all rows belonging to one
#' country are dropped simultaneously (one refit per country). This measures
#' each country's overall contribution to \eqn{\hat{\xi}} and is more
#' relevant for assessing structural sensitivity: would the conclusion
#' change if a country were excluded?
#'
#' **Interpretation of \code{dxi}.** \eqn{\Delta\xi_i = \hat{\xi}_{(-i)} -
#' \hat{\xi}}: the change in total estimated population when observation
#' (or country) \eqn{i} is removed.
#' \itemize{
#'   \item \eqn{\Delta\xi_i < 0}: dropping \eqn{i} decreases the estimate
#'     (observation was pulling the estimate up).
#'   \item \eqn{\Delta\xi_i > 0}: dropping \eqn{i} increases the estimate
#'     (observation was pulling the estimate down).
#' }
#'
#' @param object An `"uncounted"` object (fitted model).
#' @param by What to leave out: `"obs"` drops one row at a time,
#'   `"country"` drops all rows for one country at a time (requires
#'   `countries` to have been specified in the original fit).
#' @param verbose Logical; print progress?
#' @param ... Additional arguments (ignored).
#'
#' @return An object of class `"uncounted_loo"` with components:
#' \describe{
#'   \item{coefficients}{Matrix (\code{n_drops} x \code{p}): coefficients from
#'     each refit.}
#'   \item{xi}{Data frame with \eqn{\xi} estimates from each refit.}
#'   \item{dropped}{Character vector: label of dropped observation or country.}
#'   \item{dfbeta}{Matrix (\code{n_drops} x \code{p}): change in coefficients
#'     \eqn{\hat{\beta}_{(-i)} - \hat{\beta}}.}
#'   \item{dxi}{Numeric vector: change in population size
#'     \eqn{\hat{\xi}_{(-i)} - \hat{\xi}}.}
#'   \item{full_coefs}{Full-model coefficients.}
#'   \item{full_ps}{Full-model population size data frame.}
#'   \item{full_ps_total}{Full-model total \eqn{\hat{\xi}}.}
#'   \item{by}{The \code{by} argument used.}
#'   \item{converged}{Logical vector: whether each refit converged.}
#'   \item{n_drops}{Number of leave-one-out iterations.}
#' }
#'
#' @examples
#' # Simulate synthetic data
#' set.seed(42)
#' n_obs <- 15
#' sim_data <- data.frame(
#'   country = rep(paste0("C", 1:5), each = 3),
#'   year    = rep(2018:2020, 5),
#'   N       = rpois(n_obs, lambda = 500000),
#'   n       = rpois(n_obs, lambda = 1000),
#'   m       = rpois(n_obs, lambda = 50)
#' )
#'
#' fit <- estimate_hidden_pop(
#'   data = sim_data, observed = ~m, auxiliary = ~n,
#'   reference_pop = ~N, method = "poisson",
#'   countries = ~country
#' )
#'
#' # LOO by observation (drops one row at a time)
#' loo_obs <- loo(fit, by = "obs")
#' print(loo_obs)       # top 10 most influential observations
#' summary(loo_obs)     # coefficient and xi stability
#' plot(loo_obs)        # bar plot of dxi
#' plot(loo_obs, type = "coef")  # DFBETA plots per coefficient
#'
#' # LOO by country (drops all rows for one country)
#' loo_ctry <- loo(fit, by = "country")
#' print(loo_ctry)
#'
#' @export
loo <- function(object, ...) {
  UseMethod("loo")
}

#' @rdname loo
#' @export
loo.uncounted <- function(object, by = c("obs", "country"),
                          verbose = FALSE, ...) {
  by <- match.arg(by)
  data <- object$data
  call <- object$call
  n_obs <- nrow(data)

  # Determine groups to drop
  if (by == "country") {
    if (is.null(object$countries_var)) {
      stop("'by = \"country\"' requires 'countries' in the original fit.")
    }
    groups <- as.character(object$countries_var)
    unique_groups <- unique(groups)
    n_drops <- length(unique_groups)
    drop_labels <- unique_groups
  } else {
    n_drops <- n_obs
    drop_labels <- as.character(seq_len(n_obs))
  }

  # Full model results
  full_coefs <- object$coefficients
  full_ps <- popsize(object)
  p <- length(full_coefs)

  # Storage
  coef_mat <- matrix(NA_real_, nrow = n_drops, ncol = p,
                     dimnames = list(drop_labels, names(full_coefs)))
  ps_list <- vector("list", n_drops)
  converged <- logical(n_drops)

  # Extract call arguments for refitting
  # Formulas from call (safe to eval); non-formula args from fit object
  cl <- object$call
  method <- object$method
  vcov_type <- object$vcov_type
  cluster_formula <- if (!is.null(cl$cluster)) eval(cl$cluster) else NULL
  observed <- eval(cl$observed)
  auxiliary <- eval(cl$auxiliary)
  reference_pop <- eval(cl$reference_pop)
  cov_alpha <- if (!is.null(cl$cov_alpha)) eval(cl$cov_alpha) else NULL
  cov_beta <- if (!is.null(cl$cov_beta)) eval(cl$cov_beta) else NULL
  cov_gamma <- if (!is.null(cl$cov_gamma)) eval(cl$cov_gamma) else NULL
  gamma_arg <- if (object$gamma_estimated) "estimate"
               else if (!is.null(object$gamma)) object$gamma
               else NULL
  gamma_bounds <- tryCatch(
    if (!is.null(cl$gamma_bounds)) eval(cl$gamma_bounds) else c(1e-10, 0.5),
    error = function(e) c(1e-10, 0.5)
  )
  theta_start <- if (!is.null(object$theta)) object$theta else 1
  countries <- if (!is.null(cl$countries)) eval(cl$countries) else NULL
  constrained_arg <- isTRUE(object$constrained)

  warned_rank <- FALSE
  for (i in seq_len(n_drops)) {
    if (verbose && i %% 10 == 0) {
      message(sprintf("LOO %d/%d", i, n_drops))
    }

    # Determine rows to drop
    if (by == "country") {
      drop_idx <- which(groups == unique_groups[i])
    } else {
      drop_idx <- i
    }

    data_i <- data[-drop_idx, , drop = FALSE]

    # Check for rank deficiency before refitting
    rank_ok <- TRUE
    if (!is.null(cov_alpha)) {
      X_i <- model.matrix(cov_alpha, data = data_i)
      if (qr(X_i)$rank < ncol(X_i)) rank_ok <- FALSE
    }
    if (rank_ok && !is.null(cov_beta)) {
      X_i <- model.matrix(cov_beta, data = data_i)
      if (qr(X_i)$rank < ncol(X_i)) rank_ok <- FALSE
    }
    if (rank_ok && !is.null(cov_gamma)) {
      X_i <- model.matrix(cov_gamma, data = data_i)
      if (qr(X_i)$rank < ncol(X_i)) rank_ok <- FALSE
    }
    if (!rank_ok) {
      if (!warned_rank) {
        warning("LOO: dropping '", drop_labels[i],
                "' creates a rank-deficient design matrix. ",
                "Skipping refits for units that remove covariate levels.",
                call. = FALSE)
        warned_rank <- TRUE
      }
      converged[i] <- FALSE
      next
    }

    # Refit
    fit_i <- tryCatch(
      estimate_hidden_pop(
        data = data_i,
        observed = observed,
        auxiliary = auxiliary,
        reference_pop = reference_pop,
        method = method,
        cov_alpha = cov_alpha,
        cov_beta = cov_beta,
        gamma = gamma_arg,
        cov_gamma = cov_gamma,
        gamma_bounds = gamma_bounds,
        theta_start = theta_start,
        vcov = vcov_type,
        constrained = constrained_arg,
        countries = countries,
        cluster = cluster_formula
      ),
      error = function(e) NULL
    )

    if (!is.null(fit_i)) {
      coef_mat[i, ] <- fit_i$coefficients
      ps_list[[i]] <- tryCatch(popsize(fit_i), error = function(e) NULL)
      converged[i] <- TRUE
    } else {
      converged[i] <- FALSE
    }
  }

  # Combine popsize results
  ps_df <- do.call(rbind, lapply(seq_along(ps_list), function(i) {
    if (!is.null(ps_list[[i]])) {
      cbind(dropped = drop_labels[i], ps_list[[i]])
    }
  }))

  # Compute differences
  dfbeta <- sweep(coef_mat, 2, full_coefs, "-")

  # Population size differences (summed across groups if multiple)
  full_ps_total <- sum(full_ps$estimate)
  dxi <- sapply(ps_list, function(x) {
    if (!is.null(x)) sum(x$estimate) - full_ps_total else NA_real_
  })

  out <- list(
    coefficients = coef_mat,
    xi = ps_df,
    dropped = drop_labels,
    dfbeta = dfbeta,
    dxi = dxi,
    full_coefs = full_coefs,
    full_ps = full_ps,
    full_ps_total = full_ps_total,
    by = by,
    converged = converged,
    n_drops = n_drops
  )

  class(out) <- "uncounted_loo"
  out
}


# ---- S3 methods for LOO ----

#' Print LOO Sensitivity Results
#'
#' Displays a summary header (drop type, convergence count, full-model
#' \eqn{\hat{\xi}}, and the LOO range) followed by a table of the top 10
#' most influential observations/countries ranked by \eqn{|\Delta\xi|}.
#' The table shows:
#' \describe{
#'   \item{dropped}{The observation index or country label that was removed.}
#'   \item{dxi}{\eqn{\Delta\xi}: change in total population size estimate.}
#'   \item{pct_change}{Percentage change relative to full-model \eqn{\hat{\xi}}.}
#' }
#'
#' @param x An `"uncounted_loo"` object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.uncounted_loo <- function(x, ...) {
  cat("Leave-one-out sensitivity analysis\n")
  cat("Dropped by:", x$by, "\n")
  cat("N iterations:", x$n_drops, "\n")
  cat("Converged:", sum(x$converged), "/", x$n_drops, "\n\n")

  cat("Full model xi:", round(x$full_ps_total, 1), "\n")
  cat("LOO xi range:",
      round(x$full_ps_total + min(x$dxi, na.rm = TRUE), 1), "to",
      round(x$full_ps_total + max(x$dxi, na.rm = TRUE), 1), "\n\n")

  # Most influential observations
  ord <- order(abs(x$dxi), decreasing = TRUE, na.last = TRUE)
  top_n <- min(10, sum(!is.na(x$dxi)))
  top <- ord[seq_len(top_n)]

  cat("Most influential (by |delta xi|):\n")
  df <- data.frame(
    dropped = x$dropped[top],
    dxi = round(x$dxi[top], 1),
    pct_change = round(100 * x$dxi[top] / x$full_ps_total, 2)
  )
  print(df, row.names = FALSE)

  invisible(x)
}

#' Summary of LOO Sensitivity Results
#'
#' Prints coefficient stability (full-model value, LOO mean, SD, min, max
#' for each coefficient) and xi stability (full-model value, LOO mean,
#' LOO range, maximum absolute percentage change). Useful for assessing
#' whether the model is driven by a few observations or is broadly stable.
#'
#' @param object An `"uncounted_loo"` object.
#' @param ... Additional arguments (ignored).
#'
#' @export
summary.uncounted_loo <- function(object, ...) {
  cat("Leave-one-out sensitivity analysis\n")
  cat("Dropped by:", object$by, "\n")
  cat("N iterations:", object$n_drops,
      "(converged:", sum(object$converged), ")\n\n")

  cat("=== Coefficient stability ===\n")
  coef_summary <- data.frame(
    full = object$full_coefs,
    loo_mean = colMeans(object$coefficients, na.rm = TRUE),
    loo_sd = apply(object$coefficients, 2, sd, na.rm = TRUE),
    loo_min = apply(object$coefficients, 2, min, na.rm = TRUE),
    loo_max = apply(object$coefficients, 2, max, na.rm = TRUE)
  )
  print(round(coef_summary, 6))

  cat("\n=== Xi stability ===\n")
  dxi <- object$dxi[!is.na(object$dxi)]
  xi_full <- object$full_ps_total
  cat("Full xi:", round(xi_full, 1), "\n")
  cat("LOO mean:", round(xi_full + mean(dxi), 1), "\n")
  cat("LOO range:", round(xi_full + min(dxi), 1), "to",
      round(xi_full + max(dxi), 1), "\n")
  cat("Max |%change|:", round(100 * max(abs(dxi)) / xi_full, 2), "%\n")

  invisible(object)
}

#' Plot LOO Sensitivity Results
#'
#' Produces diagnostic plots for the leave-one-out analysis.
#' \describe{
#'   \item{\code{type = "xi"}}{Horizontal bar plot of \eqn{\Delta\xi} for each
#'     dropped observation/country, sorted by value. Blue bars indicate positive
#'     change (dropping the obs increases the estimate); red bars indicate
#'     negative change (dropping the obs decreases the estimate).}
#'   \item{\code{type = "coef"}}{DFBETA plots: one panel per coefficient
#'     showing \eqn{\Delta\beta_j} for each observation. Points beyond
#'     \eqn{\pm 2 \, \mathrm{SD}} are highlighted in red.}
#' }
#'
#' @param x An `"uncounted_loo"` object.
#' @param type `"xi"` (default) or `"coef"`.
#' @param ... Additional arguments passed to \code{plot}/\code{barplot}.
#'
#' @export
plot.uncounted_loo <- function(x, type = c("xi", "coef"), ...) {
  type <- match.arg(type)

  if (type == "xi") {
    dxi <- x$dxi
    ok <- !is.na(dxi)
    ord <- order(dxi[ok])
    n <- sum(ok)

    oldpar <- par(mar = c(4, 8, 3, 1))
    on.exit(par(oldpar))

    barplot(
      dxi[ok][ord],
      names.arg = x$dropped[ok][ord],
      horiz = TRUE, las = 1,
      main = "LOO: Change in xi",
      xlab = expression(Delta * xi),
      col = ifelse(dxi[ok][ord] > 0, "steelblue", "tomato"),
      border = NA,
      cex.names = 0.7
    )
    abline(v = 0, lty = 2)
  } else {
    # Coefficient DFBETA plot
    dfb <- x$dfbeta
    p <- ncol(dfb)
    oldpar <- par(mfrow = c(1, p), mar = c(4, 4, 3, 1))
    on.exit(par(oldpar))

    for (j in seq_len(p)) {
      plot(dfb[, j], pch = 16, cex = 0.6,
           main = colnames(dfb)[j],
           ylab = expression(Delta * beta),
           xlab = "Observation",
           col = ifelse(abs(dfb[, j]) > 2 * sd(dfb[, j], na.rm = TRUE),
                        "red", "grey50"))
      abline(h = 0, lty = 2)
    }
  }
}


#' Compare LOO Results from Two Models
#'
#' Takes two \code{"uncounted_loo"} objects and produces a comparison table
#' and diagnostic plots showing how each observation/country influences
#' the two models differently.
#'
#' @details
#' **Purpose.** When comparing two model specifications (e.g., Poisson vs NB,
#' or different covariate sets), it is useful to check whether the same
#' observations are influential under both models, or whether influence
#' patterns diverge.
#'
#' **Scatter plot interpretation (4 quadrants).** The scatter plot
#' (\code{type = "scatter"}) shows \eqn{\%\Delta\xi} for Model 1 (x-axis)
#' vs Model 2 (y-axis). The four quadrants reveal:
#' \itemize{
#'   \item **Top-right (+, +):** dropping this obs increases \eqn{\hat{\xi}}
#'     under both models (obs pulls estimate down in both).
#'   \item **Bottom-left (-, -):** dropping this obs decreases \eqn{\hat{\xi}}
#'     under both models (obs pulls estimate up in both).
#'   \item **Top-left (-, +):** divergent influence -- obs pulls the estimate
#'     in opposite directions across models.
#'   \item **Bottom-right (+, -):** same as top-left but reversed.
#' }
#' Points near the diagonal have similar influence under both models.
#' Points far from the diagonal have divergent influence and warrant
#' investigation.
#'
#' **Bar plot.** The bar plot (\code{type = "bar"}) shows the top N most
#' influential observations side by side for both models, ranked by
#' \code{max(|pct_1|, |pct_2|)}.
#'
#' @param loo1 An \code{"uncounted_loo"} object (first model).
#' @param loo2 An \code{"uncounted_loo"} object (second model).
#' @param labels Character vector of length 2 with model names.
#'   Defaults to \code{c("Model 1", "Model 2")}.
#' @param data Optional data.table/data.frame with columns matching
#'   \code{loo1$dropped} identifiers. If provided and LOO is by \code{"obs"},
#'   columns named \code{label_vars} are used to create readable labels.
#' @param label_vars Character vector of column names in \code{data} to use
#'   for observation labels (e.g., \code{c("country_code", "year", "sex")}).
#'
#' @return An object of class \code{"uncounted_loo_compare"} with components:
#' \describe{
#'   \item{table}{Data frame sorted by influence, with columns: \code{dropped},
#'     \code{label}, \code{dxi_<Model1>}, \code{pct_<Model1>},
#'     \code{dxi_<Model2>}, \code{pct_<Model2>}, \code{max_abs_pct}.}
#'   \item{labels}{Model labels.}
#'   \item{by}{The LOO type (\code{"obs"} or \code{"country"}).}
#'   \item{dxi_1, dxi_2}{Raw \eqn{\Delta\xi} vectors for each model.}
#'   \item{pct_1, pct_2}{Percentage \eqn{\Delta\xi} vectors.}
#'   \item{full_xi_1, full_xi_2}{Full-model \eqn{\hat{\xi}} for each model.}
#' }
#'
#' @examples
#' \donttest{
#' # Simulate synthetic data
#' set.seed(42)
#' n_obs <- 15
#' sim_data <- data.frame(
#'   country = rep(paste0("C", 1:5), each = 3),
#'   year    = rep(2018:2020, 5),
#'   N       = rpois(n_obs, lambda = 500000),
#'   n       = rpois(n_obs, lambda = 1000),
#'   m       = rpois(n_obs, lambda = 50)
#' )
#'
#' fit_po <- estimate_hidden_pop(
#'   data = sim_data, observed = ~m, auxiliary = ~n,
#'   reference_pop = ~N, method = "poisson",
#'   countries = ~country
#' )
#' fit_nb <- estimate_hidden_pop(
#'   data = sim_data, observed = ~m, auxiliary = ~n,
#'   reference_pop = ~N, method = "nb",
#'   countries = ~country
#' )
#'
#' loo_po <- loo(fit_po, by = "country")
#' loo_nb <- loo(fit_nb, by = "country")
#' comp <- compare_loo(loo_po, loo_nb, labels = c("Poisson", "NB"))
#' print(comp)
#' plot(comp, type = "scatter")
#' plot(comp, type = "bar")
#' }
#'
#' @export
compare_loo <- function(loo1, loo2,
                        labels = c("Model 1", "Model 2"),
                        data = NULL, label_vars = NULL) {
  if (!inherits(loo1, "uncounted_loo") || !inherits(loo2, "uncounted_loo"))
    stop("Both arguments must be 'uncounted_loo' objects.")

  if (loo1$by != loo2$by)
    stop("Both LOO objects must use the same 'by' (obs or country).")

  if (loo1$n_drops != loo2$n_drops)
    stop("Both LOO objects must have the same number of drops.")

  n <- loo1$n_drops
  dropped <- loo1$dropped

  # Build labels for observations
  if (!is.null(data) && !is.null(label_vars) && loo1$by == "obs") {
    obs_labels <- apply(data[, label_vars, drop = FALSE], 1, function(row) {
      paste0(row, collapse = ", ")
    })
    obs_labels <- paste0(obs_labels)
  } else {
    obs_labels <- dropped
  }

  dxi_1 <- loo1$dxi
  dxi_2 <- loo2$dxi
  pct_1 <- 100 * dxi_1 / loo1$full_ps_total
  pct_2 <- 100 * dxi_2 / loo2$full_ps_total
  max_abs_pct <- pmax(abs(pct_1), abs(pct_2), na.rm = TRUE)

  tab <- data.frame(
    dropped = dropped,
    label = obs_labels,
    dxi_1 = dxi_1,
    pct_1 = round(pct_1, 2),
    dxi_2 = dxi_2,
    pct_2 = round(pct_2, 2),
    max_abs_pct = round(max_abs_pct, 2),
    stringsAsFactors = FALSE
  )
  names(tab)[3:4] <- paste0(c("dxi_", "pct_"), labels[1])
  names(tab)[5:6] <- paste0(c("dxi_", "pct_"), labels[2])

  tab <- tab[order(-max_abs_pct), ]
  rownames(tab) <- NULL

  out <- list(
    table = tab,
    labels = labels,
    by = loo1$by,
    dxi_1 = dxi_1,
    dxi_2 = dxi_2,
    pct_1 = pct_1,
    pct_2 = pct_2,
    obs_labels = obs_labels,
    full_xi_1 = loo1$full_ps_total,
    full_xi_2 = loo2$full_ps_total
  )
  class(out) <- "uncounted_loo_compare"
  out
}


#' Print LOO Comparison Results
#'
#' Displays the full-model \eqn{\hat{\xi}} for each model and a table of
#' the top \code{n} most influential observations/countries ranked by the
#' maximum absolute percentage change across both models.
#'
#' @param x An \code{"uncounted_loo_compare"} object.
#' @param n Number of top observations to display (default 15).
#' @param ... Additional arguments (ignored).
#'
#' @export
print.uncounted_loo_compare <- function(x, n = 15, ...) {
  cat("LOO comparison:", x$labels[1], "vs", x$labels[2], "\n")
  cat("Dropped by:", x$by, "\n")
  cat("Full xi --", x$labels[1], ":", round(x$full_xi_1),
      "|", x$labels[2], ":", round(x$full_xi_2), "\n\n")

  cat("Top", n, "most influential (by max |%change|):\n")
  tab <- head(x$table, n)
  cols <- c("label", names(tab)[3:7])
  ptab <- tab[, cols]
  # Round numeric columns
  num_cols <- sapply(ptab, is.numeric)
  ptab[num_cols] <- lapply(ptab[num_cols], function(col) round(col, 2))
  print(ptab, row.names = FALSE, right = FALSE)
  invisible(x)
}


#' Plot LOO Comparison: Scatter and Barplot
#'
#' Produces diagnostic plots comparing observation influence across two models.
#'
#' @details
#' **Scatter plot (\code{type = "scatter"}).** Each point is one
#' observation/country. The x-axis is \eqn{\%\Delta\xi} for Model 1,
#' the y-axis for Model 2. The diagonal line marks equal influence.
#' Points near the origin are non-influential under both models.
#' Points far from the diagonal have divergent influence -- they affect
#' the two models differently. The top \code{label_top} most influential
#' points are labeled. See \code{\link{compare_loo}} for quadrant
#' interpretation.
#'
#' **Bar plot (\code{type = "bar"}).** Horizontal side-by-side bar chart of
#' the top \code{n} most influential observations (by maximum absolute
#' percentage change across both models). Red bars are Model 1, blue
#' bars are Model 2. Useful for quickly seeing which observations matter
#' most and whether the models agree on direction.
#'
#' @param x An \code{"uncounted_loo_compare"} object.
#' @param type \code{"scatter"} (default) or \code{"bar"}.
#' @param n Number of top observations to show in barplot (default 20).
#' @param label_top Number of points to label in scatter (default 10).
#' @param ... Additional arguments passed to \code{plot}.
#'
#' @export
plot.uncounted_loo_compare <- function(x, type = c("scatter", "bar"),
                                        n = 20, label_top = 10, ...) {
  type <- match.arg(type)

  if (type == "scatter") {
    pct1 <- x$pct_1
    pct2 <- x$pct_2
    ok <- is.finite(pct1) & is.finite(pct2)
    pct1 <- pct1[ok]; pct2 <- pct2[ok]
    labs <- x$obs_labels[ok]

    lim <- max(abs(c(pct1, pct2)), na.rm = TRUE) * 1.1
    plot(pct1, pct2,
         xlim = c(-lim, lim), ylim = c(-lim, lim),
         xlab = paste0("d.xi ", x$labels[1], " (%)"),
         ylab = paste0("d.xi ", x$labels[2], " (%)"),
         main = paste("LOO influence:", x$labels[1], "vs", x$labels[2]),
         pch = 1, col = adjustcolor("gray30", 0.5), ...)
    abline(0, 1, lty = 2, col = "gray60")
    abline(h = 0, col = "gray80")
    abline(v = 0, col = "gray80")

    # Label top influential
    max_abs <- pmax(abs(pct1), abs(pct2))
    top_idx <- order(max_abs, decreasing = TRUE)[seq_len(min(label_top, length(max_abs)))]
    text(pct1[top_idx], pct2[top_idx], labels = labs[top_idx],
         cex = 0.6, pos = sample(1:4, length(top_idx), replace = TRUE),
         col = "black")

  } else {
    # Barplot
    tab <- head(x$table, n)
    n_show <- nrow(tab)

    dxi1 <- tab[[3]]  # dxi for model 1
    dxi2 <- tab[[5]]  # dxi for model 2

    mat <- rbind(dxi1, dxi2)

    oldpar <- par(mar = c(4, 10, 3, 1))
    on.exit(par(oldpar))

    bp <- barplot(mat, beside = TRUE, horiz = TRUE,
                  names.arg = tab$label,
                  las = 1, cex.names = 0.65,
                  col = c(adjustcolor("tomato", 0.7),
                          adjustcolor("steelblue", 0.7)),
                  border = NA,
                  main = paste("LOO:", x$labels[1], "vs", x$labels[2]),
                  xlab = expression(Delta * xi))
    abline(v = 0, lty = 2)
    legend("bottomright", legend = x$labels,
           fill = c(adjustcolor("tomato", 0.7),
                    adjustcolor("steelblue", 0.7)),
           border = NA, bty = "n", cex = 0.8)
  }
}
