#' Profile Dependence Sensitivity for Population Size
#'
#' Refits the model over a grid of fixed dependence offsets and records how the
#' plug-in total population-size estimate changes.
#'
#' @param object An \code{"uncounted"} object fitted with
#'   \code{estimator = "mle"} and \code{method = "poisson"} or
#'   \code{method = "nb"}.
#' @param delta_grid Numeric vector of dependence offsets to evaluate.
#'   The baseline no-dependence case is \code{delta = 0}.
#' @param plot Logical; produce the two-panel profile plot? Default
#'   \code{TRUE}.
#' @param ... Additional arguments passed to \code{plot()}.
#'
#' @details
#' The bounded helper \code{\link{dependence_bounds}} keeps the baseline
#' point estimate fixed and only widens the identification envelope. In
#' contrast, \code{profile_dependence()} imposes a parametric dependence model
#' and refits the mean function over a fixed grid:
#'
#' \deqn{\mu_i(\delta) = \exp(\delta)\,\xi_i\rho_i,}
#'
#' or equivalently
#'
#' \deqn{\log \mu_i(\delta) = \delta + \alpha_i \log N_i + \log(\rho_i).}
#'
#' Because the model is refitted for each \code{delta}, the reported
#' population-size estimate \eqn{\hat\xi(\delta)} can move with the
#' sensitivity parameter. This is a model-based sensitivity analysis, not a
#' partial-identification bound.
#'
#' @return Invisibly, a data frame with columns:
#' \describe{
#'   \item{\code{delta}}{Fixed dependence offset used in the refit.}
#'   \item{\code{kappa}}{Multiplicative distortion \code{exp(delta)}.}
#'   \item{\code{xi}}{Plug-in total population-size estimate
#'     \code{sum(popsize(fit_i)$estimate)}.}
#'   \item{\code{loglik}}{Refitted log-likelihood, or \code{NA} when the
#'     refit failed.}
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
#' profile_dependence(fit, delta_grid = c(-0.25, 0, 0.25), plot = FALSE)
#'
#' @seealso \code{\link{dependence_bounds}}, \code{\link{robustness_dependence}}
#'
#' @export
profile_dependence <- function(object,
                               delta_grid = seq(-1, 1, length.out = 41),
                               plot = TRUE, ...) {
  if (!inherits(object, "uncounted")) {
    stop("'object' must inherit from class 'uncounted'.", call. = FALSE)
  }
  if (!identical(object$estimator, "mle")) {
    stop("profile_dependence() is only available for models fitted with ",
         "estimator = 'mle'.", call. = FALSE)
  }
  if (!object$method %in% c("poisson", "nb")) {
    stop("profile_dependence() is only available for Poisson and NB models.",
         call. = FALSE)
  }
  if (!is.numeric(delta_grid) || length(delta_grid) == 0L ||
      anyNA(delta_grid) || any(!is.finite(delta_grid))) {
    stop("'delta_grid' must be a non-empty numeric vector of finite values.",
         call. = FALSE)
  }
  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("'plot' must be either TRUE or FALSE.", call. = FALSE)
  }

  refit_args <- .refit_args_from_object(object)
  results <- data.frame(
    delta = as.numeric(delta_grid),
    kappa = exp(as.numeric(delta_grid)),
    xi = NA_real_,
    loglik = NA_real_
  )

  for (i in seq_along(delta_grid)) {
    fit_i <- tryCatch(
      .refit_uncounted_with_dependence(object, delta_grid[i], refit_args = refit_args),
      error = function(e) NULL
    )

    if (!is.null(fit_i) && isTRUE(fit_i$convergence == 0)) {
      ps <- popsize(fit_i)
      results$xi[i] <- sum(ps$estimate)
      results$loglik[i] <- fit_i$loglik
    }
  }

  if (plot) {
    op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    on.exit(par(op))

    valid_xi <- is.finite(results$xi)
    if (any(valid_xi)) {
      plot(results$delta[valid_xi], results$xi[valid_xi],
           type = "b", pch = 19, cex = 0.6,
           xlab = expression(delta), ylab = expression(hat(xi)),
           main = expression("Population size" ~ hat(xi)(delta)),
           ...)
      abline(v = 0, lty = 2, col = "red")
    }

    valid_ll <- is.finite(results$loglik)
    if (any(valid_ll)) {
      plot(results$delta[valid_ll], results$loglik[valid_ll],
           type = "b", pch = 19, cex = 0.6,
           xlab = expression(delta), ylab = "Log-likelihood",
           main = expression("Log-likelihood" ~ ell(delta)),
           ...)
      abline(v = 0, lty = 2, col = "red")
    }
  }

  invisible(results)
}

#' Dependence Robustness Value for Population Size
#'
#' Summarises the smallest dependence strength needed for the profiled
#' population-size estimate to cross a target value.
#'
#' @param object Either an \code{"uncounted"} object or a data frame returned
#'   by \code{\link{profile_dependence}} containing at least \code{delta},
#'   \code{kappa}, and \code{xi}.
#' @param q Relative change target used when \code{threshold = NULL}. For
#'   \code{direction = "decrease"}, the target is
#'   \code{baseline_xi * (1 - q)}. For \code{direction = "increase"}, the
#'   target is \code{baseline_xi * (1 + q)}.
#' @param direction Which target crossing to report: \code{"decrease"} or
#'   \code{"increase"}.
#' @param threshold Optional positive numeric target for \eqn{\xi}. When
#'   supplied, \code{q} is ignored for target construction.
#' @param delta_grid Numeric grid passed to \code{\link{profile_dependence}}
#'   when \code{object} is an \code{"uncounted"} fit.
#' @param x Object to print.
#' @param ... Additional arguments passed to print methods.
#'
#' @details
#' This helper adapts the robustness-value idea to the one-dimensional
#' dependence profile. It searches the profiled \eqn{\hat\xi(\delta)} values
#' and returns the smallest absolute \code{delta} whose profile crosses the
#' requested target. The corresponding multiplicative scales are
#' \code{rv_kappa = exp(rv_delta)} and \code{rv_Gamma = exp(abs(rv_delta))}.
#'
#' @return An object of class \code{"uncounted_dependence_robustness"} with
#'   components \code{baseline_xi}, \code{target_xi}, \code{rv_delta},
#'   \code{rv_kappa}, \code{rv_Gamma}, \code{direction}, \code{q},
#'   \code{threshold}, and \code{reached}.
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
#' robustness_dependence(fit, q = 0.10, delta_grid = seq(-0.5, 0.5, length.out = 9))
#'
#' @seealso \code{\link{profile_dependence}}, \code{\link{dependence_bounds}}
#'
#' @export
robustness_dependence <- function(object,
                                  q = 0.25,
                                  direction = c("decrease", "increase"),
                                  threshold = NULL,
                                  delta_grid = seq(-1, 1, length.out = 81)) {
  direction <- match.arg(direction)

  if (!is.numeric(q) || length(q) != 1L || is.na(q) || !is.finite(q) || q < 0) {
    stop("'q' must be a single finite number >= 0.", call. = FALSE)
  }
  if (!is.null(threshold) &&
      (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) ||
       !is.finite(threshold) || threshold <= 0)) {
    stop("'threshold' must be NULL or a single positive finite number.",
         call. = FALSE)
  }

  profile_df <- if (inherits(object, "uncounted")) {
    profile_dependence(object, delta_grid = delta_grid, plot = FALSE)
  } else if (is.data.frame(object)) {
    object
  } else {
    stop("'object' must be an 'uncounted' fit or a dependence-profile data frame.",
         call. = FALSE)
  }

  profile_df <- .validate_dependence_profile(profile_df)
  baseline_idx <- which(abs(profile_df$delta) <= sqrt(.Machine$double.eps))
  baseline_xi <- profile_df$xi[baseline_idx[1]]

  target_xi <- if (!is.null(threshold)) {
    threshold
  } else if (identical(direction, "decrease")) {
    baseline_xi * (1 - q)
  } else {
    baseline_xi * (1 + q)
  }

  hits <- if (identical(direction, "decrease")) {
    which(is.finite(profile_df$xi) & profile_df$xi <= target_xi)
  } else {
    which(is.finite(profile_df$xi) & profile_df$xi >= target_xi)
  }

  reached <- length(hits) > 0L
  if (reached) {
    ord <- order(abs(profile_df$delta[hits]), profile_df$delta[hits])
    hit <- hits[ord[1]]
    rv_delta <- profile_df$delta[hit]
    rv_kappa <- exp(rv_delta)
    rv_Gamma <- exp(abs(rv_delta))
  } else {
    rv_delta <- NA_real_
    rv_kappa <- NA_real_
    rv_Gamma <- NA_real_
  }

  out <- list(
    baseline_xi = baseline_xi,
    target_xi = target_xi,
    rv_delta = rv_delta,
    rv_kappa = rv_kappa,
    rv_Gamma = rv_Gamma,
    direction = direction,
    q = q,
    threshold = threshold,
    reached = reached
  )
  class(out) <- "uncounted_dependence_robustness"
  out
}

#' @rdname robustness_dependence
#' @export
print.uncounted_dependence_robustness <- function(x, ...) {
  cat("Dependence robustness analysis\n")
  cat("Baseline xi:", format(round(x$baseline_xi), big.mark = ","), "\n")
  cat("Target xi:", format(round(x$target_xi), big.mark = ","), "\n")

  if (isTRUE(x$reached)) {
    cat("rv_delta:", format(x$rv_delta, digits = 4), "\n")
    cat("rv_kappa:", format(x$rv_kappa, digits = 4), "\n")
    cat("rv_Gamma:", format(x$rv_Gamma, digits = 4), "\n")
  } else {
    cat("Target not reached on supplied delta grid.\n")
  }

  invisible(x)
}

#' Build refit arguments from a fitted object
#' @noRd
.refit_args_from_object <- function(object) {
  cl <- object$call

  list(
    data = object$data,
    observed = object$call_args$observed,
    auxiliary = object$call_args$auxiliary,
    reference_pop = object$call_args$reference_pop,
    method = object$method,
    cov_alpha = if (!is.null(cl$cov_alpha)) cl$cov_alpha else NULL,
    cov_beta = if (!is.null(cl$cov_beta)) cl$cov_beta else NULL,
    gamma = if (isTRUE(object$gamma_estimated)) {
      "estimate"
    } else if (!is.null(object$gamma)) {
      object$gamma
    } else {
      NULL
    },
    cov_gamma = if (!is.null(cl$cov_gamma)) cl$cov_gamma else NULL,
    gamma_bounds = if (!is.null(cl$gamma_bounds)) cl$gamma_bounds else c(1e-10, 0.5),
    theta_start = if (!is.null(object$theta)) object$theta else 1,
    link_rho = object$link_rho,
    estimator = object$estimator,
    vcov = object$vcov_type,
    weights = object$obs_weights,
    constrained = isTRUE(object$constrained),
    countries = if (!is.null(cl$countries)) cl$countries else NULL,
    cluster = if (!is.null(cl$cluster)) cl$cluster else NULL
  )
}

#' Refit a model under a fixed dependence offset
#' @noRd
.refit_uncounted_with_dependence <- function(object, delta, refit_args = NULL) {
  if (is.null(refit_args)) {
    refit_args <- .refit_args_from_object(object)
  }
  .with_dependence_offset(delta, do.call(estimate_hidden_pop, refit_args))
}

#' Validate a dependence-profile data frame
#' @noRd
.validate_dependence_profile <- function(profile_df) {
  required_cols <- c("delta", "kappa", "xi")
  if (!all(required_cols %in% names(profile_df))) {
    stop("Profile data must contain columns 'delta', 'kappa', and 'xi'.",
         call. = FALSE)
  }
  if (!is.numeric(profile_df$delta) || !is.numeric(profile_df$kappa) ||
      !is.numeric(profile_df$xi) ||
      anyNA(profile_df$delta) || any(!is.finite(profile_df$delta))) {
    stop("Profile columns 'delta', 'kappa', and 'xi' must be numeric, ",
         "with finite 'delta' values.",
         call. = FALSE)
  }
  if (!any(abs(profile_df$delta) <= sqrt(.Machine$double.eps) &
           is.finite(profile_df$xi))) {
    stop("Profile data must include a finite baseline row with delta = 0.",
         call. = FALSE)
  }
  profile_df
}
