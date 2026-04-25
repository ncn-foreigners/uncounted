#' Bootstrap Confidence Intervals for Population Size
#'
#' Computes bootstrap confidence intervals for the estimated population size
#' using the fractional weighted bootstrap (FWB) of Xu et al. (2020).
#' Supports cluster bootstrap to account for within-country correlation
#' across years and sex.
#'
#' @details
#' **Fractional weighted bootstrap (FWB).** Instead of resampling rows, FWB
#' generates random non-negative weights from a Dirichlet distribution and
#' refits the model with weighted likelihood (Xu et al., 2020). This avoids
#' duplicate-row issues with discrete data and is better suited to GLMs.
#' The \pkg{fwb} package provides the implementation.
#'
#' **Cluster bootstrap.** When \code{cluster} is specified, all observations
#' within the same cluster (e.g., country) receive the same FWB weight.
#' This preserves within-cluster correlation structure (e.g., the same
#' country observed across multiple years and sex categories). Cluster
#' bootstrap is recommended whenever observations are not independent.
#'
#' **Why bootstrap median as point estimate.** The transformation
#' \eqn{\xi = \sum N_i^{\alpha}} is convex in \eqn{\alpha} for \eqn{N > 1},
#' so by Jensen's inequality \eqn{E[\hat{\xi}] \geq \xi}. The bootstrap
#' distribution of \eqn{\hat{\xi}} is therefore typically right-skewed,
#' and the bootstrap mean overestimates \eqn{\xi}. The bootstrap median
#' is a more robust central-tendency summary that is less sensitive to
#' this skewness, and is the default point estimate.
#'
#' **Percentile CI (\code{ci_type = "perc"}).** The interval is
#' \eqn{[\hat{\xi}^{*}_{\alpha/2}, \; \hat{\xi}^{*}_{1-\alpha/2}]},
#' i.e., the \eqn{\alpha/2} and \eqn{1-\alpha/2} quantiles of the
#' bootstrap distribution.
#'
#' **Bias-corrected percentile CI (\code{ci_type = "bc"}).** Adjusts the
#' quantile probabilities for median bias. Let
#' \eqn{z_0 = \Phi^{-1}(\mathrm{Pr}(\hat{\xi}^* < \hat{\xi}_0))}.
#' The adjusted quantile probabilities are
#' \eqn{p_L = \Phi(2z_0 + z_{\alpha/2})} and
#' \eqn{p_U = \Phi(2z_0 - z_{\alpha/2})}, yielding
#' \eqn{[\hat{\xi}^{*}_{p_L}, \; \hat{\xi}^{*}_{p_U}]}.
#'
#' **Total across groups.** When multiple groups exist, the total is computed
#' by summing per-replicate group estimates before taking quantiles. This
#' correctly accounts for between-group correlation within each bootstrap
#' replicate, unlike summing marginal quantiles.
#'
#' @references
#' Xu, L., Gotwalt, C., Hong, Y., King, C. B., & Meeker, W. Q. (2020).
#' Applications of the fractional-random-weight bootstrap.
#' \emph{The American Statistician}, 74(4), 345--358.
#'
#' @param object An `"uncounted"` object (fitted model).
#' @param R Number of bootstrap replications (default 199).
#' @param cluster One-sided formula for cluster variable (e.g., `~ country_code`).
#'   When provided, FWB weights are generated at the cluster level
#'   (all observations in a cluster get the same weight).
#'   When `NULL`, observation-level FWB is used.
#' @param level Confidence level (default 0.95).
#' @param ci_type Type of bootstrap CI: `"perc"` (percentile) or `"bc"`
#'   (bias-corrected percentile). Default `"perc"`.
#' @param point_estimate Which point estimate to report:
#'   `"median"` (default, recommended) -- bootstrap median, robust to skewness;
#'   `"plugin"` -- plug-in \eqn{\hat{\xi} = \sum N^{\hat{\alpha}}};
#'   `"mean"` -- bootstrap mean.
#' @param seed Random seed for reproducibility.
#' @param verbose Print progress bar? Default TRUE.
#' @param by Optional formula for stratified population size estimation
#'   (e.g., \code{~ year}, \code{~ country}). When provided, bootstrap
#'   population sizes are computed per \code{by}-group instead of per
#'   \code{cov_alpha} group. See \code{\link{popsize}} for details.
#' @param total Logical; if \code{TRUE} and multiple groups exist, include
#'   a total row computed by summing per-replicate group estimates before
#'   taking quantiles. Default \code{FALSE}.
#'
#' @return An object of class `"uncounted_boot"` with components:
#' \describe{
#'   \item{t}{Matrix (\code{R} x \code{n_groups}): bootstrap population size
#'     estimates per group. Each row is one replicate.}
#'   \item{t0}{Numeric vector: plug-in point estimates from original fit.}
#'   \item{t0_bc}{Numeric vector: analytical bias-corrected point estimates
#'     from \code{popsize()} (multiplicative lognormal for unconstrained fits,
#'     Taylor approximation for constrained fits).}
#'   \item{popsize}{Data frame with columns: \code{group}, \code{estimate},
#'     \code{lower}, \code{upper}, where \code{estimate} uses the chosen
#'     \code{point_estimate} type.}
#'   \item{popsize_full}{Data frame with all point estimate types:
#'     \code{plugin}, \code{plugin_bc}, \code{boot_median}, \code{boot_mean},
#'     \code{lower}, \code{upper}.}
#'   \item{total}{List with total across groups (\code{plugin}, \code{plugin_bc},
#'     \code{median}, \code{mean}, \code{lower}, \code{upper}), or \code{NULL}
#'     if only one group.}
#'   \item{R}{Number of replications requested.}
#'   \item{ci_type}{CI type used (\code{"perc"} or \code{"bc"}).}
#'   \item{point_estimate}{Point estimate type used.}
#'   \item{level}{Confidence level.}
#'   \item{n_converged}{Number of bootstrap replications where the model converged.}
#'   \item{cluster}{Logical; whether cluster bootstrap was used.}
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
#' \donttest{
#' # Cluster bootstrap (recommended when data has repeated country obs)
#' boot <- bootstrap_popsize(fit, R = 99, cluster = ~country, seed = 123)
#' boot                # prints table with median as point estimate
#' boot$popsize        # data.frame for further use
#' boot$popsize_full   # all estimate types side by side
#'
#' # Bias-corrected percentile CI
#' boot_bc <- bootstrap_popsize(fit, R = 99, ci_type = "bc", seed = 123)
#' }
#'
#' @export
bootstrap_popsize <- function(object, R = 199, cluster = NULL,
                              level = 0.95, ci_type = c("perc", "bc"),
                              point_estimate = c("median", "plugin", "mean"),
                              seed = NULL, verbose = TRUE, by = NULL,
                              total = FALSE) {

  if (!inherits(object, "uncounted")) {
    stop("object must be of class 'uncounted'")
  }

  if (!requireNamespace("fwb", quietly = TRUE)) {
    stop("Package 'fwb' is required for bootstrap. Install with: install.packages('fwb')")
  }

  ci_type <- match.arg(ci_type)
  point_estimate <- match.arg(point_estimate)

  # Extract call arguments for refitting
  # Use fit object fields where available; eval formulas from call (formulas are safe)
  # Avoid eval on non-formula args (symbols like use_gamma fail in remote envs)
  cl <- object$call
  data <- object$data
  observed <- eval(cl$observed)
  auxiliary <- eval(cl$auxiliary)
  reference_pop <- eval(cl$reference_pop)
  method <- object$method
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
  constrained_arg <- isTRUE(object$constrained)
  countries <- if (!is.null(cl$countries)) eval(cl$countries) else NULL
  estimator_arg <- object$estimator
  link_rho_arg <- object$link_rho

  # Point estimates from original fit (using by if provided)
  ps0 <- popsize(object, by = by, bias_correction = TRUE, total = total)
  t0 <- ps0$estimate
  t0_bc <- ps0$estimate_bc
  group_meta <- attr(ps0, "groups")

  # If by is specified, build group factor for aggregating bootstrap replicates
  by_group_factor <- NULL
  if (!is.null(by)) {
    by_vars <- all.vars(by)
    by_data <- data[, by_vars, drop = FALSE]
    by_group_factor <- interaction(by_data, drop = TRUE, sep = ", ")
  }

  # Cluster variable
  cluster_var <- NULL
  if (!is.null(cluster)) {
    cluster_var <- eval(cluster[[2]], data)
  }

  # Pre-allocate storage for per-replicate parameters
  param_names <- names(coef(object))
  if (object$gamma_estimated && !isTRUE(object$has_cov_gamma)) {
    param_names <- c(param_names, "gamma")
  }
  if (!is.null(object$theta)) {
    param_names <- c(param_names, "theta")
  }
  bp_env <- new.env(parent = emptyenv())
  bp_env$mat <- matrix(NA_real_, nrow = R, ncol = length(param_names),
                       dimnames = list(NULL, param_names))
  bp_env$idx <- 0L

  # Statistic function for fwb
  stat_fn <- function(data, w) {
    fit_w <- tryCatch(
      uncounted::estimate_hidden_pop(
        data = data,
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
        link_rho = link_rho_arg,
        estimator = estimator_arg,
        vcov = "HC0",
        weights = w,
        constrained = constrained_arg,
        countries = countries
      ),
      error = function(e) NULL
    )

    bp_env$idx <- bp_env$idx + 1L
    if (is.null(fit_w) || !isTRUE(fit_w$convergence == 0)) {
      return(rep(NA_real_, length(t0)))
    }

    # Capture per-replicate parameters (guard against extra fwb calls)
    if (bp_env$idx <= R) {
      row_i <- coef(fit_w)
      if (fit_w$gamma_estimated && !isTRUE(fit_w$has_cov_gamma)) {
        row_i <- c(row_i, gamma = fit_w$gamma)
      }
      if (!is.null(fit_w$theta)) {
        row_i <- c(row_i, theta = fit_w$theta)
      }
      bp_env$mat[bp_env$idx, ] <- row_i
    }

    ps_w <- uncounted::popsize(fit_w, by = by, bias_correction = FALSE)
    ps_w$estimate
  }

  # Run FWB
  if (!is.null(seed)) set.seed(seed)

  boot_result <- fwb::fwb(
    data = as.data.frame(data),
    statistic = stat_fn,
    R = R,
    cluster = if (!is.null(cluster_var)) cluster_var else NULL,
    verbose = verbose
  )

  # Extract bootstrap matrix
  boot_t <- boot_result$t  # R x n_groups
  if (ncol(boot_t) == length(ps0$group)) {
    colnames(boot_t) <- ps0$group
  }

  # Compute per-group statistics
  alpha_ci <- 1 - level
  n_groups <- length(t0)

  ci_lower <- numeric(n_groups)
  ci_upper <- numeric(n_groups)
  boot_median <- numeric(n_groups)
  boot_mean <- numeric(n_groups)

  for (j in seq_len(n_groups)) {
    boot_j <- boot_t[, j]
    boot_j <- boot_j[is.finite(boot_j)]

    if (length(boot_j) < 10) {
      ci_lower[j] <- ci_upper[j] <- boot_median[j] <- boot_mean[j] <- NA_real_
      next
    }

    boot_median[j] <- median(boot_j)
    boot_mean[j] <- mean(boot_j)

    if (ci_type == "perc") {
      ci_lower[j] <- quantile(boot_j, alpha_ci / 2)
      ci_upper[j] <- quantile(boot_j, 1 - alpha_ci / 2)
    } else {
      # Bias-corrected percentile (BC)
      z0 <- qnorm(mean(boot_j < t0[j]))
      za <- qnorm(alpha_ci / 2)
      p_lower <- pnorm(2 * z0 + za)
      p_upper <- pnorm(2 * z0 - za)
      ci_lower[j] <- quantile(boot_j, p_lower)
      ci_upper[j] <- quantile(boot_j, p_upper)
    }
  }

  # Select point estimate
  est_selected <- switch(point_estimate,
    median = boot_median,
    plugin = t0,
    mean = boot_mean
  )

  # Full results (all estimate types)
  ps_full <- data.frame(
    group = ps0$group,
    plugin = t0,
    plugin_bc = t0_bc,
    boot_median = boot_median,
    boot_mean = boot_mean,
    lower = ci_lower,
    upper = ci_upper,
    stringsAsFactors = FALSE
  )

  # Main results (selected point estimate)
  ps_boot <- data.frame(
    group = ps0$group,
    estimate = est_selected,
    lower = ci_lower,
    upper = ci_upper,
    stringsAsFactors = FALSE
  )

  # Compute proper total from per-replicate sums (medians/quantiles don't add)
  total_boot <- NULL
  if (total && n_groups > 1) {
    boot_total <- rowSums(boot_t)
    boot_total <- boot_total[is.finite(boot_total)]
    if (length(boot_total) >= 10) {
      total_median <- median(boot_total)
      total_mean <- mean(boot_total)
      if (ci_type == "perc") {
        total_lower <- quantile(boot_total, alpha_ci / 2)
        total_upper <- quantile(boot_total, 1 - alpha_ci / 2)
      } else {
        z0 <- qnorm(mean(boot_total < sum(t0)))
        za <- qnorm(alpha_ci / 2)
        total_lower <- quantile(boot_total, pnorm(2 * z0 + za))
        total_upper <- quantile(boot_total, pnorm(2 * z0 - za))
      }
      # Use correct full-gradient BC total from popsize(total=TRUE)
      total_bc <- if (!is.null(attr(ps0, "total"))) {
        attr(ps0, "total")$estimate_bc
      } else {
        sum(t0_bc)
      }
      total_boot <- list(
        plugin = sum(t0), plugin_bc = total_bc,
        median = total_median, mean = total_mean,
        lower = as.numeric(total_lower),
        upper = as.numeric(total_upper)
      )
    }
  }

  # Append Total row to popsize and popsize_full when total=TRUE
  if (!is.null(total_boot)) {
    total_est_selected <- switch(point_estimate,
      median = total_boot$median,
      plugin = total_boot$plugin,
      mean = total_boot$mean
    )
    ps_boot <- rbind(ps_boot, data.frame(
      group = "Total", estimate = total_est_selected,
      lower = total_boot$lower, upper = total_boot$upper,
      stringsAsFactors = FALSE
    ))
    ps_full <- rbind(ps_full, data.frame(
      group = "Total", plugin = total_boot$plugin,
      plugin_bc = total_boot$plugin_bc,
      boot_median = total_boot$median, boot_mean = total_boot$mean,
      lower = total_boot$lower, upper = total_boot$upper,
      stringsAsFactors = FALSE
    ))
  }

  attr(ps_boot, "groups") <- .bootstrap_popsize_groups_with_total(group_meta, ps_boot$group)
  attr(ps_full, "groups") <- .bootstrap_popsize_groups_with_total(group_meta, ps_full$group)

  out <- list(
    t = boot_t,
    t0 = t0,
    t0_bc = t0_bc,
    popsize = ps_boot,
    popsize_full = ps_full,
    total = total_boot,
    R = R,
    ci_type = ci_type,
    point_estimate = point_estimate,
    level = level,
    boot_params = bp_env$mat,
    n_converged = sum(apply(boot_t, 1, function(x) all(is.finite(x)))),
    cluster = !is.null(cluster),
    groups = group_meta
  )

  class(out) <- "uncounted_boot"
  out
}


#' Bootstrap Exceedance Probability for Population Size
#'
#' Computes the empirical bootstrap tail area for a threshold question such as
#' \dQuote{What fraction of bootstrap replications imply \eqn{\xi} above a
#' given number?}
#'
#' @param object An \code{"uncounted_boot"} object returned by
#'   \code{\link{bootstrap_popsize}}.
#' @param threshold A single finite numeric threshold.
#' @param group Optional group label. When \code{NULL} (default), the helper
#'   uses the total bootstrap distribution if available; otherwise it uses the
#'   only group when the bootstrap object contains a single group.
#' @param direction Which tail area to compute: \code{"above"} for
#'   \code{P^*(xi > c)} or \code{"below"} for \code{P^*(xi < c)}.
#' @param x Object to print.
#' @param ... Additional arguments passed to print methods.
#'
#' @details
#' The reported value is an empirical bootstrap exceedance probability, i.e. a
#' sample proportion over finite bootstrap draws. It is useful for threshold
#' questions, but it is not a Bayesian posterior probability.
#'
#' For multi-group bootstrap objects with \code{total = TRUE},
#' \code{exceedance_popsize()} reconstructs the total bootstrap distribution by
#' summing the per-replicate group estimates in \code{object$t}. This preserves
#' the within-replicate dependence across groups.
#'
#' @return An object of class \code{"uncounted_popsize_exceedance"} with
#'   components \code{group}, \code{threshold}, \code{direction},
#'   \code{n_boot}, \code{n_finite}, \code{estimate}, \code{count}, and
#'   \code{distribution_summary}.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n_obs <- 15
#' sim_data <- data.frame(
#'   country = rep(paste0("C", 1:5), each = 3),
#'   year    = rep(2018:2020, 5),
#'   N       = round(exp(rnorm(n_obs, mean = 13, sd = 0.2)))
#' )
#' sim_data$n <- rpois(n_obs, lambda = pmax(1, 0.003 * sim_data$N))
#' sim_data$m <- rpois(n_obs, lambda = sim_data$N^0.6 *
#'   (0.005 + sim_data$n / sim_data$N)^0.8)
#'
#' fit <- estimate_hidden_pop(
#'   data = sim_data, observed = ~m, auxiliary = ~n,
#'   reference_pop = ~N, method = "poisson",
#'   countries = ~country
#' )
#'
#' if (requireNamespace("fwb", quietly = TRUE)) {
#'   boot <- bootstrap_popsize(
#'     fit, R = 49, cluster = ~country, seed = 123,
#'     total = TRUE, verbose = FALSE
#'   )
#'   exceedance_popsize(boot, threshold = 2000)
#' }
#' }
#'
#' @export
exceedance_popsize <- function(object, threshold, group = NULL,
                               direction = c("above", "below")) {
  direction <- match.arg(direction)

  if (!inherits(object, "uncounted_boot")) {
    stop("'object' must inherit from class 'uncounted_boot'.", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) ||
      !is.finite(threshold)) {
    stop("'threshold' must be a single finite numeric value.", call. = FALSE)
  }
  if (!is.null(group) &&
      (!is.character(group) || length(group) != 1L || is.na(group))) {
    stop("'group' must be NULL or a single character string.", call. = FALSE)
  }

  dist_info <- .bootstrap_popsize_distribution(object, group = group)
  dist <- dist_info$distribution
  n_boot <- length(dist)
  finite <- is.finite(dist)
  n_finite <- sum(finite)

  if (n_finite < 10L) {
    stop("At least 10 finite bootstrap draws are required for exceedance calculations.",
         call. = FALSE)
  }

  dist_finite <- dist[finite]
  count <- if (identical(direction, "above")) {
    sum(dist_finite > threshold)
  } else {
    sum(dist_finite < threshold)
  }
  estimate <- count / n_finite

  distribution_summary <- c(
    mean = mean(dist_finite),
    median = stats::median(dist_finite),
    sd = stats::sd(dist_finite),
    q025 = as.numeric(stats::quantile(dist_finite, 0.025)),
    q975 = as.numeric(stats::quantile(dist_finite, 0.975))
  )

  out <- list(
    group = dist_info$group,
    threshold = threshold,
    direction = direction,
    n_boot = n_boot,
    n_finite = n_finite,
    estimate = estimate,
    count = count,
    distribution_summary = distribution_summary
  )
  class(out) <- "uncounted_popsize_exceedance"
  out
}

#' @rdname exceedance_popsize
#' @export
print.uncounted_popsize_exceedance <- function(x, ...) {
  relation <- if (identical(x$direction, "above")) ">" else "<"

  cat("Bootstrap exceedance probability\n")
  cat("Group:", x$group, "\n")
  cat("P*(xi ", relation, " ", format(x$threshold, trim = TRUE), ") = ",
      format(round(x$estimate, 4), nsmall = 4), " (",
      x$count, "/", x$n_finite, " finite draws)\n", sep = "")
  cat("Bootstrap draws:", x$n_boot, "\n")
  cat("Distribution summary:\n")
  print(round(x$distribution_summary, 3))

  invisible(x)
}


# ---- S3 methods ----

#' Print Bootstrap Population Size Results
#'
#' Prints a formatted table showing all point estimate types and the
#' bootstrap confidence interval. The output columns are:
#' \describe{
#'   \item{Plugin}{Plug-in estimate \eqn{\hat{\xi} = \sum N^{\hat{\alpha}}}.}
#'   \item{Plugin (BC)}{Analytical bias-corrected plug-in estimate from
#'     \code{popsize()}.}
#'   \item{Boot median}{Median of the bootstrap distribution (recommended).}
#'   \item{Boot mean}{Mean of the bootstrap distribution.}
#'   \item{CI lower / CI upper}{Bootstrap confidence interval bounds.}
#' }
#' A Total row is added when there are multiple groups. The header
#' indicates which point estimate type and CI method were selected.
#'
#' @param x An `"uncounted_boot"` object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.uncounted_boot <- function(x, ...) {
  cat("Bootstrap population size estimation\n")
  cat("R =", x$R, "| CI type:", x$ci_type,
      "| Point estimate:", x$point_estimate,
      "| Converged:", x$n_converged, "/", x$R, "\n")
  if (x$cluster) cat("Cluster bootstrap\n")
  cat(sprintf("%.0f%% CI\n\n", x$level * 100))

  .print_popsize_table_boot(x$popsize_full, x$ci_type, x$point_estimate, x$total,
                            show_total = !is.null(x$total))

  invisible(x)
}

#' Summary of Bootstrap Population Size Results
#'
#' Prints the same table as \code{\link{print.uncounted_boot}}, followed by
#' a per-group summary of the bootstrap distribution: mean, standard deviation,
#' 2.5\% quantile, median (50\%), and 97.5\% quantile. This is useful for
#' assessing the shape of the bootstrap distribution (e.g., skewness) and
#' comparing the spread across groups.
#'
#' @param object An `"uncounted_boot"` object.
#' @param ... Additional arguments (ignored).
#'
#' @export
summary.uncounted_boot <- function(object, ...) {
  print(object)

  cat("\nBootstrap distribution summary:\n")
  boot_summary <- apply(object$t, 2, function(col) {
    col <- col[is.finite(col)]
    c(mean = mean(col), sd = sd(col),
      q025 = quantile(col, 0.025),
      q50 = quantile(col, 0.5),
      q975 = quantile(col, 0.975))
  })
  group_names <- object$popsize$group[object$popsize$group != "Total"]
  colnames(boot_summary) <- group_names
  print(round(t(boot_summary)))

  if (!is.null(object$boot_params)) {
    cat("\nParameter bootstrap summary:\n")
    param_summary <- apply(object$boot_params, 2, function(col) {
      col <- col[is.finite(col)]
      if (length(col) == 0) return(rep(NA_real_, 5))
      c(mean = mean(col), sd = sd(col),
        q025 = quantile(col, 0.025),
        q50 = quantile(col, 0.5),
        q975 = quantile(col, 0.975))
    })
    print(round(t(param_summary), 6))
  }

  invisible(object)
}

#' Format bootstrap popsize table for printing
#' @noRd
.print_popsize_table_boot <- function(ps_full, ci_type = "perc",
                                      point_estimate = "median",
                                      total = NULL,
                                      show_total = FALSE) {
  fmt <- function(x) {
    ifelse(is.finite(x) & abs(x) < 1e12,
           format(round(x), big.mark = ",", trim = TRUE),
           ifelse(is.finite(x),
                  formatC(x, format = "e", digits = 2),
                  "NA"))
  }

  tab <- data.frame(
    `Plugin` = fmt(ps_full$plugin),
    `Plugin (BC)` = fmt(ps_full$plugin_bc),
    `Boot median` = fmt(ps_full$boot_median),
    `Boot mean` = fmt(ps_full$boot_mean),
    `CI lower` = fmt(ps_full$lower),
    `CI upper` = fmt(ps_full$upper),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  rownames(tab) <- ps_full$group

  # Total row: already present in ps_full when total=TRUE was used,

  # otherwise fall back to manual summing for display
  if ("Total" %in% ps_full$group) {
    # Already in ps_full — no extra action needed (it's in tab already)
  } else if (!is.null(total)) {
    tab <- rbind(tab, data.frame(
      `Plugin` = fmt(total$plugin),
      `Plugin (BC)` = fmt(total$plugin_bc),
      `Boot median` = fmt(total$median),
      `Boot mean` = fmt(total$mean),
      `CI lower` = fmt(total$lower),
      `CI upper` = fmt(total$upper),
      check.names = FALSE, stringsAsFactors = FALSE,
      row.names = "Total"
    ))
  } else if (show_total && nrow(ps_full) > 1) {
    tab <- rbind(tab, data.frame(
      `Plugin` = fmt(sum(ps_full$plugin)),
      `Plugin (BC)` = fmt(sum(ps_full$plugin_bc)),
      `Boot median` = fmt(sum(ps_full$boot_median)),
      `Boot mean` = fmt(sum(ps_full$boot_mean)),
      `CI lower` = fmt(sum(ps_full$lower)),
      `CI upper` = fmt(sum(ps_full$upper)),
      check.names = FALSE, stringsAsFactors = FALSE,
      row.names = "Total"
    ))
  }

  ci_label <- if (ci_type == "bc") "bias-corrected percentile" else "percentile"
  pe_label <- switch(point_estimate,
    median = "bootstrap median (recommended)",
    plugin = "plug-in",
    mean = "bootstrap mean"
  )
  cat(sprintf("  Point estimate: %s | CI: bootstrap %s\n", pe_label, ci_label))
  print(tab)
}

#' Extract a bootstrap population-size distribution for exceedance summaries
#' @noRd
.bootstrap_popsize_distribution <- function(object, group = NULL) {
  valid_groups <- setdiff(object$popsize$group, "Total")

  if (is.null(group)) {
    if (!is.null(object$total) || "Total" %in% object$popsize$group) {
      return(list(group = "Total", distribution = rowSums(object$t)))
    }
    if (length(valid_groups) == 1L) {
      return(list(group = valid_groups[1], distribution = object$t[, 1]))
    }
    stop("When multiple groups are present, use bootstrap_popsize(..., total = TRUE) ",
         "or specify 'group'.", call. = FALSE)
  }

  if (!(group %in% valid_groups)) {
    stop("'group' must match one of the non-total bootstrap groups.",
         call. = FALSE)
  }

  idx <- match(group, valid_groups)
  list(group = group, distribution = object$t[, idx])
}

#' Add a Total metadata row when a bootstrap popsize table has one
#' @noRd
.bootstrap_popsize_groups_with_total <- function(groups, labels) {
  if (is.null(groups)) {
    groups <- data.frame(.group = labels[labels != "Total"],
                         group = labels[labels != "Total"],
                         stringsAsFactors = FALSE)
  }
  out <- groups
  if ("Total" %in% labels && !("Total" %in% out$.group)) {
    total_row <- out[NA_integer_, , drop = FALSE][1, , drop = FALSE]
    total_row$.group <- "Total"
    if ("group" %in% names(total_row)) {
      if (is.factor(out$group)) {
        out$group <- factor(as.character(out$group),
                            levels = unique(c(levels(out$group), "Total")))
        total_row$group <- factor("Total", levels = levels(out$group))
      } else {
        total_row$group <- "Total"
      }
    }
    out <- rbind(out, total_row)
  }
  out
}
