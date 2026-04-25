#' Dependence Bounds for Total Population Size
#'
#' Quantifies how the estimated total hidden population changes when the
#' identifying separability assumption behind
#' \eqn{E(m_i \mid N_i, n_i) = \xi_i \rho_i} is relaxed.
#'
#' @param object An \code{"uncounted"} object.
#' @param Gamma Numeric vector of sensitivity parameters. Values must be
#'   finite and at least 1. Duplicates are removed and the remaining values are
#'   sorted increasingly before the bounds are computed.
#' @param bias_correction Logical; use the bias-corrected total from
#'   \code{\link{popsize}}? Default \code{TRUE}. When \code{FALSE}, the
#'   plug-in total is used instead.
#' @param threshold Optional positive numeric threshold used to compute a
#'   tipping-point \code{gamma_star}. Default \code{NULL}.
#' @param threshold_side Which bound should be compared against
#'   \code{threshold}? \code{"lower"} (default) finds the first
#'   \code{Gamma} value with \code{lower <= threshold}; \code{"upper"} finds
#'   the first value with \code{upper >= threshold}.
#' @param x Object to print.
#' @param ... Additional arguments passed to print methods.
#'
#' @details
#' The fitted baseline model writes the expected observed count as
#'
#' \deqn{\mu_i = \xi_i \rho_i,}
#'
#' where \eqn{\xi_i} is the latent population component and \eqn{\rho_i} is the
#' detection component. This function relaxes the separability assumption by
#' introducing a multiplicative distortion term:
#'
#' \deqn{\mu_i = \xi_i \rho_i \kappa_i.}
#'
#' The sensitivity model imposes the bound
#'
#' \deqn{1 / \Gamma \le \kappa_i \le \Gamma, \qquad \Gamma \ge 1.}
#'
#' When \eqn{\Gamma = 1}, the distortion disappears and the returned bounds
#' collapse to the baseline total. Larger values of \eqn{\Gamma} widen the
#' identification envelope around the baseline estimate.
#'
#' The reported bounds are a sensitivity analysis for an untestable identifying
#' assumption. They are not confidence intervals and should not be merged with
#' the sampling uncertainty returned by \code{\link{popsize}}.
#'
#' @return An object of class \code{"uncounted_dependence_bounds"}, a list with:
#' \describe{
#'   \item{\code{table}}{Data frame with columns \code{Gamma},
#'     \code{estimate}, \code{lower}, \code{upper},
#'     \code{pct_change_lower}, and \code{pct_change_upper}.}
#'   \item{\code{baseline}}{One-row data frame containing the baseline total
#'     columns \code{estimate}, \code{estimate_bc}, \code{se},
#'     \code{lower}, and \code{upper}.}
#'   \item{\code{bias_correction}}{Logical; whether the baseline used the
#'     analytical bias correction.}
#'   \item{\code{threshold}}{The supplied threshold, or \code{NULL}.}
#'   \item{\code{threshold_side}}{Which bound was compared against the
#'     threshold.}
#'   \item{\code{gamma_star}}{\code{NULL} when no threshold was supplied;
#'     otherwise the first \code{Gamma} on the supplied grid that crosses the
#'     threshold, or \code{NA_real_} if the grid never crosses it.}
#' }
#'
#' @examples
#' set.seed(123)
#' d <- data.frame(
#'   N = round(exp(rnorm(20, mean = 7, sd = 0.35)))
#' )
#' d$n <- rpois(20, lambda = pmax(1, 0.03 * d$N))
#' d$m <- rpois(20, lambda = d$N^0.5 * (0.01 + d$n / d$N)^0.8)
#'
#' fit <- estimate_hidden_pop(
#'   data = d, observed = ~m, auxiliary = ~n, reference_pop = ~N,
#'   method = "poisson", gamma = 0.01
#' )
#'
#' dependence_bounds(fit, Gamma = c(1, 1.1, 1.25))
#' dependence_bounds(
#'   fit,
#'   Gamma = c(1, 1.1, 1.25),
#'   threshold = 900,
#'   threshold_side = "lower"
#' )
#'
#' @export
dependence_bounds <- function(object,
# nolint start: object_name_linter.
                              Gamma = c(1, 1.05, 1.10, 1.25, 1.50, 2.00),
# nolint end: object_name_linter.
                              bias_correction = TRUE,
                              threshold = NULL,
                              threshold_side = c("lower", "upper")) {
  threshold_side <- match.arg(threshold_side)

  if (!inherits(object, "uncounted")) {
    stop("'object' must inherit from class 'uncounted'.", call. = FALSE)
  }
  if (!is.logical(bias_correction) || length(bias_correction) != 1L ||
      is.na(bias_correction)) {
    stop("'bias_correction' must be either TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(Gamma) || length(Gamma) == 0L || anyNA(Gamma) ||
      any(!is.finite(Gamma))) {
    stop("'Gamma' must be a non-empty numeric vector of finite values >= 1.",
         call. = FALSE)
  }
  if (any(Gamma < 1)) {
    stop("'Gamma' values must all be >= 1.", call. = FALSE)
  }
  if (!is.null(threshold) &&
      (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) ||
       !is.finite(threshold) || threshold <= 0)) {
    stop("'threshold' must be NULL or a single positive finite number.",
         call. = FALSE)
  }

  Gamma <- sort(unique(as.numeric(Gamma)))

  ps <- popsize(object, bias_correction = bias_correction, total = TRUE)
  baseline <- .popsize_total_baseline(ps)
  estimate_col <- if (isTRUE(bias_correction)) "estimate_bc" else "estimate"
  estimate <- baseline[[estimate_col]][1]

  if (!is.finite(estimate) || estimate <= 0) {
    stop("Baseline total estimate must be finite and positive.", call. = FALSE)
  }

  tab <- data.frame(
    Gamma = Gamma,
    estimate = rep.int(estimate, length(Gamma)),
    lower = estimate / Gamma,
    upper = estimate * Gamma,
    pct_change_lower = 100 * ((estimate / Gamma) / estimate - 1),
    pct_change_upper = 100 * ((estimate * Gamma) / estimate - 1)
  )

  gamma_star <- NULL
  if (!is.null(threshold)) {
    bound <- if (identical(threshold_side, "lower")) tab$lower else tab$upper
    hit <- if (identical(threshold_side, "lower")) {
      which(bound <= threshold)
    } else {
      which(bound >= threshold)
    }
    gamma_star <- if (length(hit) > 0L) tab$Gamma[min(hit)] else NA_real_
  }

  out <- list(
    table = tab,
    baseline = baseline,
    bias_correction = bias_correction,
    threshold = threshold,
    threshold_side = threshold_side,
    gamma_star = gamma_star
  )
  class(out) <- c("uncounted_dependence_bounds", "uncounted_sensitivity")
  out
}

#' @rdname dependence_bounds
#' @export
sensitivity_dependence <- function(object,
# nolint start: object_name_linter.
                                   Gamma = c(1, 1.05, 1.10, 1.25, 1.50, 2.00),
# nolint end: object_name_linter.
                                   bias_correction = TRUE,
                                   threshold = NULL,
                                   threshold_side = c("lower", "upper")) {
  .Deprecated("dependence_bounds")
  dependence_bounds(
    object = object,
    Gamma = Gamma,
    bias_correction = bias_correction,
    threshold = threshold,
    threshold_side = threshold_side
  )
}

#' @rdname dependence_bounds
#' @export
print.uncounted_dependence_bounds <- function(x, ...) {
  baseline <- x$table$estimate[1]
  cat("Dependence bounds analysis\n")
  cat("Baseline total:", format(round(baseline), big.mark = ","), "\n")
  cat("Interpretation: bounds reflect identification sensitivity, not sampling uncertainty.\n\n")
  print.data.frame(x$table, row.names = FALSE, ...)

  if (!is.null(x$threshold)) {
    gamma_label <- if (is.na(x$gamma_star)) {
      "not crossed on supplied Gamma grid"
    } else {
      format(x$gamma_star, trim = TRUE)
    }
    cat("\nThreshold (", x$threshold_side, " bound): ",
        format(x$threshold, trim = TRUE), "\n", sep = "")
    cat("gamma_star: ", gamma_label, "\n", sep = "")
  }

  invisible(x)
}

#' @rdname dependence_bounds
#' @export
print.uncounted_sensitivity <- function(x, ...) {
  print.uncounted_dependence_bounds(x, ...)
}

#' Extract a total-like baseline row from popsize() output
#' @noRd
.popsize_total_baseline <- function(ps) {
  tot <- attr(ps, "total")
  if (!is.null(tot)) {
    return(data.frame(
      estimate = tot$estimate,
      estimate_bc = tot$estimate_bc,
      se = tot$se,
      lower = tot$lower,
      upper = tot$upper
    ))
  }

  data.frame(
    estimate = ps$estimate[1],
    estimate_bc = ps$estimate_bc[1],
    se = NA_real_,
    lower = ps$lower[1],
    upper = ps$upper[1]
  )
}
