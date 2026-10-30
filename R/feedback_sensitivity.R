#' Detection-Feedback Sensitivity Analysis for Population Size
#'
#' Relaxes the "no feedback" assumption of the main model, under which detection
#' does not depend on the size of the unauthorized population, and reports the
#' population size consistent with each degree of feedback. The analysis is a
#' deterministic transform of the fitted cell-level estimates: nothing is refitted.
#'
#' @param fit An \code{"uncounted"} object.
#' @param anchor How to choose \eqn{K_i}, the size at which the main model's
#'   detection function is assumed correct. There is no default, because the
#'   anchor can reverse the direction of the result (see Details).
#'   \code{"geomean"} uses the geometric mean of \eqn{\xi_{0i}} over all cells
#'   (or within \code{anchor_by}), so "large" means large in absolute terms.
#'   \code{"share"} uses \eqn{K_i = \bar{s} N_i} with
#'   \eqn{\bar{s}} the geometric mean of \eqn{\xi_{0i} / N_i}, so "large" means
#'   large relative to the reference population. \code{"value"} uses \code{K}.
#' @param delta_grid Numeric vector of \eqn{\delta} values. Must contain 0, and
#'   every value must exceed \eqn{-1}.
#' @param K Anchor for \code{anchor = "value"}: a positive scalar, or a positive
#'   vector with one value per observation.
#' @param anchor_by Optional one-sided formula: compute the \code{"geomean"} or
#'   \code{"share"} anchor within groups. With \code{~ country_code} the levels of
#'   countries barely move, so the analysis concerns changes within countries.
#' @param by Optional one-sided formula for the reported groups, with the same
#'   meaning as in \code{\link{popsize}}.
#' @param total Logical; add a \code{"Total"} row when there are several groups?
#' @param boot Optional \code{\link{bootstrap_popsize}} result for the same fit.
#'   Its replicate coefficients are transformed; the bootstrap is not re-run. The
#'   anchor stays at its point-estimate value in every replicate, because it is
#'   part of the assumption rather than an estimate.
#' @param level Level of the percentile bootstrap intervals.
#' @param sublinear Logical; apply constraint A1 (see Details)?
#' @param p_max Upper limit for the implied detection in constraint A2.
#' @param keep_cells Logical; return cell-level matrices in \code{$cells}?
#' @param delta_prior Optional Monte Carlo sensitivity analysis (Greenland, 2005):
#'   a function of \code{R} returning \code{R} draws of \eqn{\delta}, or
#'   \code{TRUE} for a uniform distribution on the admissible range. One draw is
#'   used per bootstrap replicate, so \code{boot} is required. The data do not
#'   update this distribution, because \eqn{\delta} is not identified: the result
#'   depends completely on the distribution chosen.
#' @param x Object to print or plot.
#' @param group For \code{plot()}: labels of the groups to draw. Defaults to
#'   \code{"Total"} when available, otherwise the first six groups.
#' @param plausible For \code{plot()}: optional length-2 range of \eqn{\delta} to
#'   mark as plausible.
#' @param ... Additional arguments passed to \code{plot()}.
#'
#' @details
#' \strong{Model.} The main model assumes
#' \eqn{E(p_i \mid M_i, N_i, n_i) = \rho_i}. Here \eqn{\delta} is the elasticity of
#' detection with respect to the unauthorized population \eqn{M_i}:
#'
#' \deqn{E(p_i \mid M_i, N_i, n_i) = \rho_i (M_i / K_i)^\delta, \qquad \delta > -1,}
#'
#' which implies \eqn{\mu_i \approx \rho_i \xi_i (\xi_i / K_i)^\delta}. The
#' approximation treats \eqn{M_i} as equal to its mean. Holding the fitted
#' \eqn{\hat\mu_i} and \eqn{\hat\rho_i} fixed, the population size consistent with
#' \eqn{\delta} is
#'
#' \deqn{\xi_i(\delta) = K_i (\xi_{0i} / K_i)^{1 / (1 + \delta)},
#'   \qquad \xi_{0i} = N_i^{\hat\alpha_i}.}
#'
#' \strong{Fit preservation.} Fitted values, likelihood and residuals are the same
#' for every \eqn{\delta}: only the split of \eqn{\hat\mu_i} into size and detection
#' changes. So \eqn{\delta} is not identified and the data cannot choose among its
#' values (compare Link, 2003; Barker et al., 2018). Refitting with an offset would
#' not help. The offset model is
#' \eqn{\log\mu_i = (1+\delta)\alpha_i \log N_i + \beta_i \log r_i - \delta \log K},
#' and because the model has no intercept, the likelihood varies with \eqn{\delta}
#' only through the normalization \eqn{-\delta \log K}: with \eqn{K = 1} it does not
#' vary at all.
#'
#' \strong{The anchor is an assumption.} \eqn{K_i} is a pivot: for \eqn{\delta > 0},
#' cells above it shrink and cells below it grow. Totals are dominated by large
#' cells, so the direction depends on which side of the pivot they lie. Under
#' \code{"geomean"} large countries are far above \eqn{K} and the total falls with
#' \eqn{\delta}. Under \code{"share"} they usually have a below-typical
#' \eqn{\xi_{0i} / N_i}, so they lie below \eqn{K_i} and the total can rise. With
#' \code{"geomean"} computed over exactly the summed cells, the total is
#' non-increasing in \eqn{\delta}. Report results under more than one anchor.
#'
#' \strong{Admissibility}, from point estimates only:
#' \itemize{
#'   \item A1 (\code{sublinear = TRUE}, constant anchors only):
#'     \eqn{\delta > \max_i \hat\alpha_i - 1}, so that size stays sublinear in
#'     \eqn{N}. It does not apply to \code{"share"} or to a \code{K} that varies
#'     across cells. With \code{anchor_by}, the anchor is constant within groups and
#'     A1 is applied.
#'   \item A2: implied detection \eqn{\hat\mu_i / \xi_i(\delta) \le} \code{p_max}.
#'     Cells that already exceed \code{p_max} at \eqn{\delta = 0} (possible because
#'     the power link does not bound \eqn{\rho_i}) would make every \eqn{\delta}
#'     inadmissible, so they are left out of A2 with a warning and listed in
#'     \code{$a2_excluded}.
#'   \item A3: for every group, \eqn{\sum_{i \in g} \xi_i(\delta) \ge
#'     \sum_{i \in g} m_i}.
#' }
#' \code{$admissible_range} is the largest run of consecutive admissible grid points
#' containing 0, so its ends are limited by \code{delta_grid}. The range can hinge on
#' very small cells; \code{$diagnostics} shows which constraint binds and where.
#'
#' \strong{Point estimates.} \code{estimate} uses the plug-in \eqn{\xi_{0i}}, because
#' the analytical bias correction in \code{\link{popsize}} was derived for
#' \eqn{\delta = 0}. With a bootstrap, \code{boot_median} is the recommended point
#' estimate, as in \code{\link{bootstrap_popsize}}.
#'
#' @return An object of class \code{"uncounted_feedback"}, a list with:
#' \describe{
#'   \item{\code{profile}}{Data frame with columns \code{delta}, \code{group},
#'     \code{estimate}, \code{boot_median}, \code{lower}, \code{upper} and
#'     \code{admissible}.}
#'   \item{\code{admissible_range}}{\code{c(delta_L, delta_U)}, or \code{NA}s when
#'     \eqn{\delta = 0} is inadmissible.}
#'   \item{\code{anchor}}{List with the \code{type}, the anchor \code{K} per cell, and
#'     \code{by}.}
#'   \item{\code{diagnostics}}{Data frame per \eqn{\delta}: \code{A1}, \code{A2},
#'     \code{A3}, \code{binding}, \code{max_p}, \code{max_p_row}, and
#'     \code{d_total}, the derivative of the total with respect to \eqn{\delta}.}
#'   \item{\code{a2_excluded}}{Data frame of the cells left out of A2.}
#'   \item{\code{cells}}{With \code{keep_cells = TRUE}: matrices \code{xi},
#'     \code{detection} and \code{alpha} (cells by \eqn{\delta}). Implied
#'     \code{alpha} is \eqn{\hat\alpha_i / (1 + \delta)} for constant anchors,
#'     \eqn{1 - (1 - \hat\alpha_i) / (1 + \delta)} for \code{"share"}, and \code{NA}
#'     for \code{"value"}.}
#'   \item{\code{monte_carlo}}{With \code{delta_prior}: data frame with \code{group},
#'     \code{median}, \code{lower} and \code{upper}, and the draws in
#'     \code{attr(, "delta")}.}
#' }
#'
#' @references
#' Barker, R. J., Schofield, M. R., Link, W. A., and Sauer, J. R. (2018). On the
#' reliability of N-mixture models for count data. \emph{Biometrics}, 74(1),
#' 369--377. \doi{10.1111/biom.12734}
#'
#' Greenland, S. (2005). Multiple-bias modelling for analysis of observational data.
#' \emph{Journal of the Royal Statistical Society: Series A}, 168(2), 267--306.
#' \doi{10.1111/j.1467-985X.2004.00349.x}
#'
#' Link, W. A. (2003). Nonidentifiability of population size from capture-recapture
#' data with heterogeneous detection probabilities. \emph{Biometrics}, 59(4),
#' 1123--1130. \doi{10.1111/j.0006-341X.2003.00129.x}
#'
#' @examples
#' set.seed(123)
#' d <- data.frame(N = round(exp(rnorm(40, mean = 7, sd = 1))))
#' d$n <- rpois(40, lambda = pmax(1, 0.03 * d$N))
#' d$m <- rpois(40, lambda = d$N^0.6 * (0.01 + d$n / d$N)^0.6)
#' fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.01)
#'
#' feedback_sensitivity(fit, anchor = "geomean")
#' feedback_sensitivity(fit, anchor = "share")
#'
#' @seealso \code{\link{feedback_identified_set}}, \code{\link{feedback_tipping}},
#'   \code{\link{feedback_calibrate}}, \code{\link{dependence_bounds}},
#'   \code{\link{profile_dependence}} (where \code{delta} is a different quantity, a
#'   constant shift of the mean)
#'
#' @export
feedback_sensitivity <- function(fit,
                                 anchor,
                                 delta_grid = seq(-0.3, 0.3, by = 0.01),
                                 K = NULL,
                                 anchor_by = NULL,
                                 by = NULL,
                                 total = FALSE,
                                 boot = NULL,
                                 level = 0.95,
                                 sublinear = TRUE,
                                 p_max = 1,
                                 keep_cells = FALSE,
                                 delta_prior = NULL) {
  if (!inherits(fit, "uncounted")) {
    stop("'fit' must inherit from class 'uncounted'.", call. = FALSE)
  }
  if (missing(anchor)) {
    stop("'anchor' must be supplied: \"geomean\", \"share\" or \"value\". It has no ",
         "default because it can reverse the direction of the result; see ",
         "?feedback_sensitivity.", call. = FALSE)
  }
  anchor <- match.arg(anchor, c("geomean", "share", "value"))
  delta_grid <- .feedback_delta_grid(delta_grid)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(p_max) || length(p_max) != 1L || is.na(p_max) || p_max <= 0) {
    stop("'p_max' must be a single positive number.", call. = FALSE)
  }

  log_N <- log(fit$N)
  alpha <- fit$alpha_values
  log_xi0 <- alpha * log_N
  mu <- fit$fitted.values
  n_cells <- length(log_xi0)

  log_K <- .feedback_log_anchor(fit, anchor, log_xi0, log_N, K, anchor_by)
  constant_anchor <- anchor == "geomean" ||
    (anchor == "value" && length(unique(log_K)) == 1L)

  info <- .popsize_group_info(fit, by = by)
  group_id <- integer(n_cells)
  for (g in seq_along(info$index)) group_id[info$index[[g]]] <- g
  labels <- names(info$index)
  add_total <- isTRUE(total) && length(labels) > 1L

  replicates <- if (!is.null(boot)) .feedback_replicates(fit, boot, log_N) else NULL

  # ---- Admissibility, from point estimates only ----
  a2_excluded <- .feedback_a2_excluded(fit, mu, log_xi0, p_max)
  a2_cells <- setdiff(seq_len(n_cells), a2_excluded$row)
  m_group <- as.numeric(rowsum(fit$m, group_id))
  L <- log_xi0 - log_K

  diagnostics <- do.call(rbind, lapply(delta_grid, function(d) {
    xi <- exp(log_K + L / (1 + d))
    p <- mu / xi
    A1 <- if (isTRUE(sublinear) && constant_anchor) d > max(alpha) - 1 else NA
    p_kept <- p[a2_cells]
    A2 <- length(p_kept) == 0L || max(p_kept) <= p_max
    A3 <- all(as.numeric(rowsum(xi, group_id)) >= m_group)
    failed <- c(if (isFALSE(A1)) "A1", if (!A2) "A2", if (!A3) "A3")
    data.frame(
      delta = d, A1 = A1, A2 = A2, A3 = A3,
      binding = if (length(failed)) paste(failed, collapse = "+") else "none",
      max_p = if (length(p_kept)) max(p_kept) else NA_real_,
      max_p_row = if (length(p_kept)) a2_cells[which.max(p_kept)] else NA_integer_,
      d_total = -sum(xi * L) / (1 + d)^2
    )
  }))
  admissible <- diagnostics$binding == "none"
  admissible_range <- .feedback_admissible_range(delta_grid, admissible)

  # ---- Profile ----
  profile <- do.call(rbind, lapply(seq_along(delta_grid), function(j) {
    d <- delta_grid[j]
    est <- .feedback_totals(log_xi0, log_K, d, group_id, add_total)
    out <- data.frame(delta = d, group = c(labels, if (add_total) "Total"),
                      estimate = est, boot_median = NA_real_,
                      lower = NA_real_, upper = NA_real_,
                      admissible = admissible[j], stringsAsFactors = FALSE)
    if (!is.null(replicates)) {
      draws <- .feedback_totals(replicates, log_K, d, group_id, add_total)
      for (g in seq_len(nrow(out))) {
        x_g <- draws[g, ]
        x_g <- x_g[is.finite(x_g)]
        if (length(x_g) >= 10L) {
          ci <- .bootstrap_ci(x_g, est[g], level, "perc")
          out$boot_median[g] <- stats::median(x_g)
          out$lower[g] <- ci[1]
          out$upper[g] <- ci[2]
        }
      }
    }
    out
  }))
  rownames(profile) <- NULL

  # ---- Optional outputs ----
  cells <- NULL
  if (isTRUE(keep_cells)) {
    xi_mat <- exp(outer(L, 1 / (1 + delta_grid)) + log_K)
    alpha_mat <- switch(anchor,
      share = 1 - outer(1 - alpha, 1 / (1 + delta_grid)),
      value = matrix(NA_real_, n_cells, length(delta_grid)),
      outer(alpha, 1 / (1 + delta_grid))
    )
    colnames(xi_mat) <- colnames(alpha_mat) <- format(delta_grid, trim = TRUE)
    cells <- list(xi = xi_mat, detection = mu / xi_mat, alpha = alpha_mat)
  }

  monte_carlo <- NULL
  if (!is.null(delta_prior) && !isFALSE(delta_prior)) {
    monte_carlo <- .feedback_monte_carlo(delta_prior, replicates, log_K, group_id,
                                         add_total, labels, level, admissible_range)
  }

  out <- list(
    profile = profile,
    cells = cells,
    admissible_range = admissible_range,
    anchor = list(type = anchor, K = exp(log_K), by = anchor_by),
    diagnostics = diagnostics,
    a2_excluded = a2_excluded,
    monte_carlo = monte_carlo,
    delta_grid = delta_grid,
    level = level,
    total = add_total,
    labels = labels,
    groups = info$groups,
    group_id = group_id,
    log_xi0 = log_xi0,
    log_K = log_K,
    replicates = replicates
  )
  class(out) <- "uncounted_feedback"
  out
}

#' @rdname feedback_sensitivity
#' @export
print.uncounted_feedback <- function(x, ...) {
  cat("Detection-feedback sensitivity analysis\n")
  K <- x$anchor$K
  cat("Anchor: ", x$anchor$type,
      if (!is.null(x$anchor$by)) paste0(" within ", deparse(x$anchor$by)),
      if (length(unique(K)) == 1L) paste0(" (K = ", signif(K[1], 4), ")")
      else paste0(" (K from ", signif(min(K), 4), " to ", signif(max(K), 4), ")"),
      "\n", sep = "")
  cat("The anchor is an assumption and can change the direction of the result.\n")
  rng <- x$admissible_range
  if (anyNA(rng)) {
    cat("Admissible range: none (delta = 0 is inadmissible)\n")
  } else {
    cat("Admissible range: [", format(rng[1]), ", ", format(rng[2]), "]\n", sep = "")
  }
  if (nrow(x$a2_excluded) > 0L) {
    cat("Cells left out of A2:", nrow(x$a2_excluded), "(see $a2_excluded)\n")
  }

  show <- c(-0.2, -0.1, 0, 0.1, 0.2)
  show <- x$delta_grid[vapply(x$delta_grid, function(d) any(abs(d - show) < 1e-8), logical(1))]
  if (!anyNA(rng)) show <- show[show >= rng[1] & show <= rng[2]]
  tab <- x$profile[x$profile$delta %in% show, , drop = FALSE]
  has_boot <- any(is.finite(tab$boot_median))
  groups <- unique(tab$group)
  n_show <- min(length(groups), 15L)
  if (nrow(tab) > 0L) {
    cat("\nPlug-in estimates", if (has_boot) " [bootstrap percentile interval]", ":\n", sep = "")
    cells <- ifelse(
      is.finite(tab$lower),
      sprintf("%s [%s, %s]", .feedback_fmt(tab$estimate), .feedback_fmt(tab$lower),
              .feedback_fmt(tab$upper)),
      .feedback_fmt(tab$estimate)
    )
    mat <- matrix(cells, nrow = length(groups), byrow = FALSE,
                  dimnames = list(groups, paste0("d=", format(show, trim = TRUE))))
    print(noquote(mat[seq_len(n_show), , drop = FALSE]))
    if (length(groups) > n_show) cat("... and", length(groups) - n_show, "more groups\n")
  }
  if (!is.null(x$monte_carlo)) {
    cat("\nMonte Carlo sensitivity analysis over delta:\n")
    mc <- x$monte_carlo
    mc[c("median", "lower", "upper")] <- lapply(mc[c("median", "lower", "upper")], round)
    print(as.data.frame(mc), row.names = FALSE)
  }
  invisible(x)
}

#' @rdname feedback_sensitivity
#' @export
plot.uncounted_feedback <- function(x, group = NULL, plausible = NULL, ...) {
  labels <- unique(x$profile$group)
  if (is.null(group)) {
    group <- if ("Total" %in% labels) "Total" else utils::head(labels, 6L)
  }
  unknown <- setdiff(group, labels)
  if (length(unknown) > 0L) {
    stop("Unknown group(s): ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  n <- length(group)
  op <- par(mfrow = c(ceiling(n / min(n, 3L)), min(n, 3L)), mar = c(4, 4, 3, 1))
  on.exit(par(op))

  for (g in group) {
    tab <- x$profile[x$profile$group == g, , drop = FALSE]
    adm <- tab$admissible
    keep <- if (any(adm)) adm else rep(TRUE, nrow(tab))
    ## Scale to the admissible part, which inadmissible totals would otherwise dwarf
    ylim <- range(tab$estimate[keep], tab$lower[keep], tab$upper[keep], na.rm = TRUE)
    plot(tab$delta, tab$estimate, type = "n", ylim = ylim, xlab = expression(delta),
         ylab = expression(hat(xi)), main = g, ...)
    usr <- par("usr")
    if (!is.null(plausible)) {
      rect(plausible[1], usr[3], plausible[2], usr[4], border = NA,
           col = adjustcolor("forestgreen", alpha.f = 0.12))
    }
    step <- if (nrow(tab) > 1L) min(diff(tab$delta)) / 2 else 0
    for (j in which(!adm)) {
      rect(tab$delta[j] - step, usr[3], tab$delta[j] + step, usr[4], border = NA,
           col = adjustcolor("grey50", alpha.f = 0.35))
    }
    if (any(is.finite(tab$lower))) {
      ok <- is.finite(tab$lower)
      graphics::polygon(c(tab$delta[ok], rev(tab$delta[ok])), c(tab$lower[ok], rev(tab$upper[ok])),
              border = NA, col = adjustcolor("steelblue", alpha.f = 0.25))
    }
    lines(tab$delta, tab$estimate, lwd = 2)
    abline(v = 0, lty = 2, col = "red")
  }
  invisible(x)
}


# ---- Internal helpers ----

#' Validate delta_grid; values within rounding of 0 become exactly 0
#' @noRd
.feedback_delta_grid <- function(delta_grid) {
  if (!is.numeric(delta_grid) || length(delta_grid) == 0L ||
      anyNA(delta_grid) || any(!is.finite(delta_grid))) {
    stop("'delta_grid' must be a non-empty numeric vector of finite values.",
         call. = FALSE)
  }
  if (any(delta_grid <= -1)) {
    stop("All 'delta_grid' values must be greater than -1.", call. = FALSE)
  }
  delta_grid[abs(delta_grid) < sqrt(.Machine$double.eps)] <- 0
  if (!any(delta_grid == 0)) {
    stop("'delta_grid' must contain 0, the main model.", call. = FALSE)
  }
  sort(unique(delta_grid))
}

#' Log anchor per cell
#' @noRd
.feedback_log_anchor <- function(fit, anchor, log_xi0, log_N, K, anchor_by) {
  n_cells <- length(log_xi0)
  if (anchor == "value") {
    if (is.null(K)) {
      stop("anchor = \"value\" requires 'K'.", call. = FALSE)
    }
    if (!is.null(anchor_by)) {
      stop("'anchor_by' applies to \"geomean\" and \"share\", not \"value\".",
           call. = FALSE)
    }
    if (!is.numeric(K) || anyNA(K) || any(!is.finite(K)) || any(K <= 0)) {
      stop("'K' must contain finite positive values.", call. = FALSE)
    }
    if (!length(K) %in% c(1L, n_cells)) {
      stop("'K' must be a scalar or have one value per observation.", call. = FALSE)
    }
    return(rep_len(log(K), n_cells))
  }
  if (!is.null(K)) {
    stop("'K' is only used with anchor = \"value\".", call. = FALSE)
  }
  group_id <- rep(1L, n_cells)
  if (!is.null(anchor_by)) {
    index <- .popsize_group_info(fit, by = anchor_by)$index
    for (g in seq_along(index)) group_id[index[[g]]] <- g
  }
  if (anchor == "geomean") {
    stats::ave(log_xi0, group_id)
  } else {
    stats::ave(log_xi0 - log_N, group_id) + log_N
  }
}

#' Cells above p_max already at delta = 0, which A2 must leave out
#' @noRd
.feedback_a2_excluded <- function(fit, mu, log_xi0, p_max) {
  rows <- which(mu / exp(log_xi0) > p_max)
  out <- data.frame(row = rows, N = fit$N[rows], n = fit$n_aux[rows], m = fit$m[rows],
                    fitted = mu[rows], detection = (mu / exp(log_xi0))[rows])
  if (!is.null(fit$countries_var)) {
    out$country <- as.character(fit$countries_var)[rows]
  }
  if (length(rows) > 0L) {
    shown <- if (is.null(out$country)) out$row else paste0(out$country, " (row ", out$row, ")")
    warning(length(rows), " cell(s) already have implied detection above p_max = ", p_max,
            " at delta = 0 and are left out of constraint A2: ",
            paste(utils::head(shown, 5L), collapse = ", "),
            if (length(rows) > 5L) ", ...", ". See $a2_excluded.", call. = FALSE)
  }
  out
}

#' Largest run of consecutive admissible grid points containing 0
#' @noRd
.feedback_admissible_range <- function(delta_grid, admissible) {
  zero <- which(delta_grid == 0)
  if (!admissible[zero]) {
    warning("delta = 0 is inadmissible, so there is no admissible range. Check ",
            "$diagnostics; A1 fails when some alpha is at least 1.", call. = FALSE)
    return(c(NA_real_, NA_real_))
  }
  lo <- zero
  while (lo > 1L && admissible[lo - 1L]) lo <- lo - 1L
  hi <- zero
  while (hi < length(delta_grid) && admissible[hi + 1L]) hi <- hi + 1L
  if (sum(admissible) > hi - lo + 1L) {
    warning("The admissible set is not contiguous; $admissible_range is the run ",
            "containing 0.", call. = FALSE)
  }
  delta_grid[c(lo, hi)]
}

#' Replicate log xi0 (cells by replicates) from the bootstrap's coefficients
#' @noRd
.feedback_replicates <- function(fit, boot, log_N) {
  if (!inherits(boot, "uncounted_boot")) {
    stop("'boot' must come from bootstrap_popsize().", call. = FALSE)
  }
  coef_names <- colnames(fit$X_alpha)
  params <- boot$boot_params
  if (is.null(params) || !all(coef_names %in% colnames(params))) {
    stop("'boot' does not hold the alpha coefficients of 'fit'.", call. = FALSE)
  }
  xi_fit <- sum(fit$N^fit$alpha_values)
  if (!isTRUE(all.equal(sum(boot$t0), xi_fit, tolerance = 1e-6))) {
    stop("'boot' does not come from the same fit: its plug-in total is ",
         signif(sum(boot$t0), 6), ", the fit's is ", signif(xi_fit, 6), ".", call. = FALSE)
  }
  params <- params[, coef_names, drop = FALSE]
  params <- params[stats::complete.cases(params), , drop = FALSE]
  alpha <- fit$X_alpha %*% t(params)
  if (isTRUE(fit$constrained)) alpha[] <- .inv_logit(alpha)
  alpha * log_N
}

#' Group totals at one delta, for a vector or a cells-by-replicates matrix of log xi0
#' @noRd
.feedback_totals <- function(log_xi0, log_K, delta, group_id, add_total) {
  xi <- exp(log_K + (log_xi0 - log_K) / (1 + delta))
  by_group <- rowsum(xi, group_id)
  if (is.matrix(log_xi0)) {
    dimnames(by_group) <- NULL
    return(if (add_total) rbind(by_group, colSums(by_group)) else by_group)
  }
  by_group <- as.numeric(by_group)
  if (add_total) c(by_group, sum(by_group)) else by_group
}

#' One delta per replicate, so the interval reflects both sources of uncertainty
#' @noRd
.feedback_monte_carlo <- function(delta_prior, replicates, log_K, group_id, add_total,
                                  labels, level, admissible_range) {
  if (is.null(replicates)) {
    stop("'delta_prior' requires 'boot'.", call. = FALSE)
  }
  n_rep <- ncol(replicates)
  if (isTRUE(delta_prior)) {
    if (anyNA(admissible_range)) {
      stop("There is no admissible range to draw delta from; supply a function.",
           call. = FALSE)
    }
    delta_prior <- function(R) stats::runif(R, admissible_range[1], admissible_range[2])
  }
  if (!is.function(delta_prior)) {
    stop("'delta_prior' must be a function, TRUE, or NULL.", call. = FALSE)
  }
  draws <- delta_prior(n_rep)
  if (!is.numeric(draws) || length(draws) != n_rep || anyNA(draws) || any(draws <= -1)) {
    stop("'delta_prior' must return ", n_rep, " numeric values above -1.", call. = FALSE)
  }
  xi <- exp(sweep(replicates - log_K, 2, 1 / (1 + draws), "*") + log_K)
  totals <- rowsum(xi, group_id)
  dimnames(totals) <- NULL
  if (add_total) totals <- rbind(totals, colSums(totals))
  a <- (1 - level) / 2
  out <- data.frame(
    group = c(labels, if (add_total) "Total"),
    median = apply(totals, 1, stats::median),
    lower = apply(totals, 1, stats::quantile, probs = a, names = FALSE),
    upper = apply(totals, 1, stats::quantile, probs = 1 - a, names = FALSE),
    stringsAsFactors = FALSE
  )
  attr(out, "delta") <- draws
  out
}

#' @noRd
.feedback_fmt <- function(x) format(round(x), big.mark = ",", trim = TRUE, scientific = FALSE)


#' Identified Set for Population Size under Detection Feedback
#'
#' Bounds the population size over a range of \eqn{\delta}, with a confidence
#' interval for the true value or for the whole set.
#'
#' @param x An \code{"uncounted_feedback"} object.
#' @param delta_range Length-2 range of \eqn{\delta}. Defaults to
#'   \code{x$admissible_range}.
#' @param level Confidence level.
#' @param method \code{"imbens_manski"} covers the true \eqn{\xi_g} (Imbens and
#'   Manski, 2004). \code{"union"} covers the whole set, so it is wider.
#' @param se How to measure the spread of the replicate bounds on the log scale:
#'   \code{"sd"}, or \code{"iqr"} for \eqn{\mathrm{IQR} / 1.349}.
#' @param ... Ignored.
#'
#' @details
#' With \eqn{D} = \code{delta_range}, the bounds are
#' \eqn{L_g = \min_{\delta \in D} \xi_g(\delta)} and
#' \eqn{U_g = \max_{\delta \in D} \xi_g(\delta)}. They are found on a grid and refined
#' with \code{optimize()}, not read from the ends of \eqn{D}, because group totals
#' need not be monotone in \eqn{\delta}. Replicate bounds are replicate-wise minima
#' and maxima over \eqn{D}.
#'
#' \code{"imbens_manski"}, on the log scale with
#' \eqn{\Delta = \log U_g - \log L_g}, solves for \eqn{c} in
#' \deqn{\Phi\left(c + \frac{\Delta}{\max(s_L, s_U)}\right) - \Phi(-c) = \mathrm{level},}
#' which lies between \code{qnorm(level)} and \code{qnorm((1 + level) / 2)}, and
#' returns \eqn{\exp(\log L_g - c s_L)} and \eqn{\exp(\log U_g + c s_U)}.
#' \code{"union"} returns the \eqn{(1 - \mathrm{level}) / 2} quantile of the replicate
#' \eqn{L_g} and the \eqn{(1 + \mathrm{level}) / 2} quantile of the replicate \eqn{U_g}.
#'
#' \strong{Validity.} The Imbens--Manski interval requires superefficient estimation
#' of \eqn{\Delta} (Stoye, 2009). Stoye's Lemma 3 gives this when the bounds are jointly
#' asymptotically normal and \eqn{U_g \ge L_g} by construction, as here. Normality
#' can fail when the minimum or maximum over \eqn{D} is attained at more than one
#' \eqn{\delta}, so treat the intervals with caution if \code{delta_at_min} or
#' \code{delta_at_max} is interior to \eqn{D}. \eqn{D} itself is treated as fixed.
#'
#' @return A data frame of class \code{"uncounted_feedback_set"} with one row per group
#'   and columns \code{group}, \code{L}, \code{U}, \code{delta_at_min},
#'   \code{delta_at_max}, \code{ci_lower}, \code{ci_upper}, \code{c}, \code{Delta},
#'   \code{se_L} and \code{se_U}. Without a bootstrap in \code{x} the interval columns
#'   are \code{NA}.
#'
#' @references
#' Imbens, G. W. and Manski, C. F. (2004). Confidence intervals for partially
#' identified parameters. \emph{Econometrica}, 72(6), 1845--1857.
#' \doi{10.1111/j.1468-0262.2004.00555.x}
#'
#' Stoye, J. (2009). More on confidence intervals for partially identified
#' parameters. \emph{Econometrica}, 77(4), 1299--1315. \doi{10.3982/ECTA7347}
#'
#' @seealso \code{\link{feedback_sensitivity}}
#'
#' @export
feedback_identified_set <- function(x, delta_range = NULL, level = 0.95,
                                    method = c("imbens_manski", "union"),
                                    se = c("sd", "iqr")) {
  if (!inherits(x, "uncounted_feedback")) {
    stop("'x' must come from feedback_sensitivity().", call. = FALSE)
  }
  method <- match.arg(method)
  se <- match.arg(se)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number between 0 and 1.", call. = FALSE)
  }
  D <- if (is.null(delta_range)) x$admissible_range else delta_range
  if (!is.numeric(D) || length(D) != 2L || anyNA(D) || D[1] > D[2] || D[1] <= -1) {
    stop("The delta range must be two increasing values above -1. x has no ",
         "admissible range if delta = 0 is inadmissible; then supply 'delta_range'.",
         call. = FALSE)
  }
  grid <- sort(unique(c(seq(D[1], D[2], length.out = 41L),
                        x$delta_grid[x$delta_grid >= D[1] & x$delta_grid <= D[2]])))

  members <- split(seq_along(x$group_id), x$group_id)
  labels <- x$labels
  if (isTRUE(x$total)) {
    members <- c(members, list(seq_along(x$group_id)))
    labels <- c(labels, "Total")
  }

  rows <- lapply(seq_along(members), function(g) {
    idx <- members[[g]]
    lo <- .feedback_extreme(x$log_xi0[idx], x$log_K[idx], grid, minimum = TRUE)
    hi <- .feedback_extreme(x$log_xi0[idx], x$log_K[idx], grid, minimum = FALSE)
    out <- data.frame(group = labels[g], L = lo$value, U = hi$value,
                      delta_at_min = lo$delta, delta_at_max = hi$delta,
                      ci_lower = NA_real_, ci_upper = NA_real_, c = NA_real_,
                      Delta = log(hi$value) - log(lo$value),
                      se_L = NA_real_, se_U = NA_real_, stringsAsFactors = FALSE)
    if (is.null(x$replicates)) return(out)

    L_b <- apply(x$replicates[idx, , drop = FALSE], 2, function(a) {
      .feedback_extreme(a, x$log_K[idx], grid, minimum = TRUE)$value
    })
    U_b <- apply(x$replicates[idx, , drop = FALSE], 2, function(a) {
      .feedback_extreme(a, x$log_K[idx], grid, minimum = FALSE)$value
    })
    spread <- if (se == "sd") stats::sd else function(v) stats::IQR(v) / 1.349
    out$se_L <- spread(log(L_b))
    out$se_U <- spread(log(U_b))
    if (method == "imbens_manski") {
      out$c <- .im_critical(out$Delta, max(out$se_L, out$se_U), level)
      out$ci_lower <- exp(log(out$L) - out$c * out$se_L)
      out$ci_upper <- exp(log(out$U) + out$c * out$se_U)
    } else {
      out$ci_lower <- unname(stats::quantile(L_b, (1 - level) / 2))
      out$ci_upper <- unname(stats::quantile(U_b, (1 + level) / 2))
    }
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "delta_range") <- D
  attr(out, "level") <- level
  attr(out, "method") <- method
  class(out) <- c("uncounted_feedback_set", "data.frame")
  out
}

#' @rdname feedback_identified_set
#' @export
print.uncounted_feedback_set <- function(x, ...) {
  D <- attr(x, "delta_range")
  ## Subsetting drops the attributes and possibly columns: show it as plain data
  if (is.null(D) || !all(c("group", "L", "U") %in% names(x))) {
    print(as.data.frame(unclass(x), stringsAsFactors = FALSE), ...)
    return(invisible(x))
  }
  cat("Identified set over delta in [", signif(D[1], 3), ", ", signif(D[2], 3), "]\n", sep = "")
  has_ci <- any(is.finite(x$ci_lower))
  if (has_ci) {
    cat(100 * attr(x, "level"), "% interval: ", attr(x, "method"), "\n", sep = "")
  }
  tab <- data.frame(group = x$group, L = .feedback_fmt(x$L), U = .feedback_fmt(x$U),
                    delta_at_min = signif(x$delta_at_min, 3),
                    delta_at_max = signif(x$delta_at_max, 3))
  if (has_ci) {
    tab$ci_lower <- .feedback_fmt(x$ci_lower)
    tab$ci_upper <- .feedback_fmt(x$ci_upper)
  }
  print(tab, row.names = FALSE)
  invisible(x)
}

#' Minimum or maximum of a total over delta: grid, then refined when interior
#' @noRd
.feedback_extreme <- function(log_xi0, log_K, grid, minimum = TRUE) {
  total <- function(d) sum(exp(log_K + (log_xi0 - log_K) / (1 + d)))
  values <- vapply(grid, total, numeric(1))
  best <- if (minimum) which.min(values) else which.max(values)
  if (best == 1L || best == length(grid)) {
    return(list(delta = grid[best], value = values[best]))
  }
  opt <- stats::optimize(total, grid[c(best - 1L, best + 1L)], maximum = !minimum,
                         tol = 1e-8)
  list(delta = if (minimum) opt$minimum else opt$maximum, value = opt$objective)
}

#' Imbens-Manski critical value for an identified set of length delta
#' @noRd
.im_critical <- function(delta, se_max, level) {
  upper <- stats::qnorm((1 + level) / 2)
  lower <- stats::qnorm(level)
  if (!is.finite(se_max) || se_max <= 0) {
    return(if (delta > 0) lower else upper)
  }
  ratio <- delta / se_max
  f <- function(c) stats::pnorm(c + ratio) - stats::pnorm(-c) - level
  if (f(lower) >= 0) return(lower)
  if (f(upper) <= 0) return(upper)
  stats::uniroot(f, c(lower, upper), tol = 1e-10)$root
}


#' Tipping Points under Detection Feedback
#'
#' Finds how much detection feedback it takes to overturn a statement about the
#' population size.
#'
#' @param x An \code{"uncounted_feedback"} object.
#' @param threshold A single value for a group total. Use with \code{group}. The tested
#'   direction is the side of the threshold on which the estimate lies at
#'   \eqn{\delta = 0}.
#' @param group Group label for \code{threshold}. Defaults to \code{"Total"} when
#'   available, or the only group.
#' @param hypothesis A comparison in the \code{xi[...]} syntax of
#'   \code{\link{hypotheses_popsize}}, such as
#'   \code{"xi[year == 2024] - xi[year == 2019] > 0"}. The filters refer to the
#'   \code{by} variables of \code{x}.
#' @param level Level of the one-sided bootstrap test: the statement is supported at
#'   \eqn{\delta} when the \code{1 - level} quantile of the replicate contrasts is
#'   above 0 (for \code{>}), or the \code{level} quantile is below 0 (for \code{<}).
#'   For \code{=}, both tails use \code{(1 - level) / 2}.
#' @param ... Ignored.
#'
#' @details
#' Two kinds of tipping point are reported on each side of 0, within the
#' \code{delta_grid} of \code{x}:
#' \itemize{
#'   \item \code{point}: the \eqn{\delta} nearest 0 at which the point estimate of the
#'     contrast changes sign, found exactly between grid points.
#'   \item \code{test}: with a bootstrap, the grid value nearest 0 at which the test
#'     stops supporting the statement.
#' }
#' If the test does not support the statement at \eqn{\delta = 0}, there is nothing to
#' overturn: \code{test} is \code{NA} and \code{supported_at_zero} is \code{FALSE}.
#' Tipping points are reported whether or not they are admissible, with the admissible
#' range shown alongside.
#'
#' @return An object of class \code{"uncounted_feedback_tipping"}, a list with
#'   \code{statement}, \code{alternative}, \code{baseline} (the contrast at
#'   \eqn{\delta = 0}), \code{point} and \code{test} (each \code{c(negative, positive)},
#'   \code{NA} when not reached), \code{supported_at_zero} (\code{NA} without a
#'   bootstrap), \code{admissible_range}, and \code{curve}, a data frame with
#'   \code{delta}, \code{contrast}, \code{lower}, \code{upper} and \code{supported}.
#'
#' @seealso \code{\link{feedback_sensitivity}}, \code{\link{hypotheses_popsize}}
#'
#' @export
feedback_tipping <- function(x, threshold = NULL, group = NULL, hypothesis = NULL,
                             level = 0.95) {
  if (!inherits(x, "uncounted_feedback")) {
    stop("'x' must come from feedback_sensitivity().", call. = FALSE)
  }
  if (is.null(threshold) == is.null(hypothesis)) {
    stop("Supply either 'threshold' or 'hypothesis'.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number between 0 and 1.", call. = FALSE)
  }

  if (!is.null(hypothesis)) {
    if (!is.character(hypothesis) || length(hypothesis) != 1L) {
      stop("'hypothesis' must be a single character string.", call. = FALSE)
    }
    parsed <- .pop_hyp_parse_comparison(hypothesis)
    contrast <- function(v) {
      .pop_hyp_eval_numeric(parsed$lhs, v, x$groups) -
        .pop_hyp_eval_numeric(parsed$rhs, v, x$groups)
    }
    statement <- hypothesis
    alternative <- parsed$alternative
  } else {
    if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold)) {
      stop("'threshold' must be a single finite number.", call. = FALSE)
    }
    if (is.null(group)) {
      group <- if (isTRUE(x$total)) "Total" else if (length(x$labels) == 1L) x$labels
    }
    if (is.null(group) || !group %in% c(x$labels, if (isTRUE(x$total)) "Total")) {
      stop("'group' must name one group: ",
           paste(c(x$labels, if (isTRUE(x$total)) "Total"), collapse = ", "), ".",
           call. = FALSE)
    }
    pick <- if (group == "Total") seq_along(x$labels) else match(group, x$labels)
    contrast <- function(v) sum(v[pick]) - threshold
    at_zero <- contrast(.feedback_totals(x$log_xi0, x$log_K, 0, x$group_id, FALSE))
    alternative <- if (at_zero > 0) "greater" else "less"
    statement <- paste0("xi[", group, "] ", if (at_zero > 0) ">" else "<", " ",
                        format(threshold, big.mark = ",", scientific = FALSE))
  }

  point_at <- function(d) {
    contrast(.feedback_totals(x$log_xi0, x$log_K, d, x$group_id, FALSE))
  }
  grid <- x$delta_grid
  values <- vapply(grid, point_at, numeric(1))
  zero <- which(grid == 0)
  baseline <- values[zero]

  curve <- data.frame(delta = grid, contrast = values, lower = NA_real_,
                      upper = NA_real_, supported = NA)
  a <- if (alternative == "two.sided") (1 - level) / 2 else 1 - level
  if (!is.null(x$replicates)) {
    for (j in seq_along(grid)) {
      totals <- .feedback_totals(x$replicates, x$log_K, grid[j], x$group_id, FALSE)
      draws <- apply(totals, 2, contrast)
      draws <- draws[is.finite(draws)]
      curve$lower[j] <- unname(stats::quantile(draws, a))
      curve$upper[j] <- unname(stats::quantile(draws, 1 - a))
      curve$supported[j] <- switch(alternative,
        greater = curve$lower[j] > 0,
        less = curve$upper[j] < 0,
        two.sided = curve$lower[j] > 0 || curve$upper[j] < 0
      )
    }
  }

  ## Sign change nearest 0 on one side; an exact 0 at baseline has no side to cross from
  crossing <- function(side) {
    if (baseline == 0) return(NA_real_)
    path <- if (side < 0) rev(seq_len(zero)) else seq(zero, length(grid))
    for (k in seq_along(path)[-1L]) {
      if (sign(values[path[k]]) != sign(baseline)) {
        ends <- sort(grid[path[c(k - 1L, k)]])
        return(stats::uniroot(point_at, ends, tol = 1e-8)$root)
      }
    }
    NA_real_
  }
  stops <- function(side) {
    if (is.null(x$replicates) || !isTRUE(curve$supported[zero])) return(NA_real_)
    path <- if (side < 0) rev(seq_len(zero)) else seq(zero, length(grid))
    hit <- path[!curve$supported[path]]
    if (length(hit)) grid[hit[1]] else NA_real_
  }

  out <- list(
    statement = statement,
    alternative = alternative,
    baseline = baseline,
    point = c(negative = crossing(-1), positive = crossing(1)),
    test = c(negative = stops(-1), positive = stops(1)),
    supported_at_zero = if (is.null(x$replicates)) NA else isTRUE(curve$supported[zero]),
    admissible_range = x$admissible_range,
    level = level,
    curve = curve
  )
  class(out) <- "uncounted_feedback_tipping"
  out
}

#' @rdname feedback_tipping
#' @export
print.uncounted_feedback_tipping <- function(x, ...) {
  fmt <- function(v) {
    ifelse(is.na(v), "not reached", format(round(v, 4), trim = TRUE))
  }
  cat("Detection-feedback tipping points\n")
  cat("Statement:", x$statement, "\n")
  cat("Contrast at delta = 0:", .feedback_fmt(x$baseline), "\n")
  rng <- x$admissible_range
  cat("Admissible range: ",
      if (anyNA(rng)) "none" else paste0("[", format(rng[1]), ", ", format(rng[2]), "]"),
      "\n", sep = "")
  cat("Point estimate changes sign at delta: ", fmt(x$point["negative"]), " (below 0), ",
      fmt(x$point["positive"]), " (above 0)\n", sep = "")
  if (isTRUE(x$supported_at_zero)) {
    cat("Test at ", 100 * x$level, "% stops supporting it at delta: ",
        fmt(x$test["negative"]), " (below 0), ", fmt(x$test["positive"]),
        " (above 0)\n", sep = "")
  } else if (isFALSE(x$supported_at_zero)) {
    cat("The bootstrap test at ", 100 * x$level, "% does not support the statement at ",
        "delta = 0, so there is no significant result to overturn.\n", sep = "")
  }
  invisible(x)
}


#' Calibrate the Feedback Parameter
#'
#' Summarizes \eqn{s_i = \xi_{0i} / (N_i + \xi_{0i})}, the share of the unauthorized
#' population in the foreign population of a cell. If inspections are spread over the
#' whole foreign population, the implied mechanical congestion elasticity of detection
#' is \eqn{-s}, which helps to choose a \code{delta_grid} for
#' \code{\link{feedback_sensitivity}}.
#'
#' @param fit An \code{"uncounted"} object.
#' @param x Object to print.
#' @param ... Ignored.
#'
#' @return An object of class \code{"uncounted_feedback_calibration"}, a list with
#'   \code{s} (one value per cell), \code{mean_s}, \code{weighted_mean_s} (weighted by
#'   \eqn{\xi_{0i}}) and \code{congestion_elasticity}, the negative of
#'   \code{weighted_mean_s}.
#'
#' @seealso \code{\link{feedback_sensitivity}}
#'
#' @export
feedback_calibrate <- function(fit) {
  if (!inherits(fit, "uncounted")) {
    stop("'fit' must inherit from class 'uncounted'.", call. = FALSE)
  }
  xi0 <- fit$N^fit$alpha_values
  s <- xi0 / (fit$N + xi0)
  weighted <- stats::weighted.mean(s, xi0)
  out <- list(s = s, mean_s = mean(s), weighted_mean_s = weighted,
              congestion_elasticity = -weighted)
  class(out) <- "uncounted_feedback_calibration"
  out
}

#' @rdname feedback_calibrate
#' @export
print.uncounted_feedback_calibration <- function(x, ...) {
  cat("Feedback calibration: s = xi0 / (N + xi0)\n")
  cat("Cell mean of s:", round(x$mean_s, 3), "\n")
  cat("xi0-weighted mean of s:", round(x$weighted_mean_s, 3), "\n")
  cat("Implied mechanical congestion elasticity:", round(x$congestion_elasticity, 3), "\n")
  invisible(x)
}
