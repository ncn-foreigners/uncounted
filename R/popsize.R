#' Estimated Population Size
#'
#' Computes the estimated total unauthorized population
#' \eqn{\hat{\xi} = \sum_{i} N_i^{\hat{\alpha}}} for each group
#' defined by the covariates in alpha. Includes bias correction via
#' a multiplicative lognormal adjustment for unconstrained models (and a
#' second-order Taylor approximation for constrained models), plus log-normal
#' delta-method confidence intervals on the population-size scale.
#'
#' @details
#' **Point estimate.** For a group \eqn{g} with reference populations
#' \eqn{N_1, \ldots, N_n} and estimated exponent \eqn{\hat{\alpha}_g},
#' the plug-in population size is
#' \deqn{\hat{\xi}_g = \sum_{i=1}^{n} N_i^{\hat{\alpha}_g}.}
#'
#' **Bias correction.** Because \eqn{\xi} is a nonlinear
#' function of \eqn{\hat{\alpha}}, the plug-in estimate is biased upward
#' (by Jensen's inequality, \eqn{E[N^{\hat{\alpha}}] \geq N^{E[\hat{\alpha}]}}
#' when \eqn{N > 1}).
#' A second-order Taylor expansion of \eqn{h(\alpha) = \sum_i N_i^{\alpha}}
#' around \eqn{\alpha_0} gives the approximate bias
#' \deqn{
#'   \mathrm{Bias}(\hat{\xi}_g) \approx
#'     \frac{1}{2} \sum_{i=1}^{n} N_i^{\alpha_g} (\log N_i)^2 \;
#'     \mathbf{x}_g' \mathbf{V} \mathbf{x}_g,
#' }
#' where \eqn{\mathbf{x}_g} is the design vector for group \eqn{g} and
#' \eqn{\mathbf{V}} is the variance-covariance matrix of \eqn{\hat{\alpha}}.
#' For unconstrained models, the exact multiplicative correction is used:
#' \deqn{
#'   \hat{\xi}^{BC}_g = \sum_{i=1}^{n} N_i^{\hat{\alpha}_g}
#'     \exp\!\left(-\frac{1}{2} (\log N_i)^2 \mathbf{x}_g' \mathbf{V}
#'     \mathbf{x}_g\right),
#' }
#' which is exact under normality of \eqn{\hat{\alpha}} and always positive.
#' For constrained models the subtractive Taylor correction is used instead
#' (the logistic-normal integral has no closed form), and the bias can be
#' positive or negative depending on \eqn{\alpha}.
#' Model-based variance is used for bias correction rather than HC-robust
#' variance. For Poisson and NB count models this means the Fisher-style
#' inverse information for the mean-model parameters, evaluated at the fitted
#' coefficients. This remains true when the coefficients were obtained by
#' \code{estimator = "gmm"} or \code{estimator = "el"}: the robust HC or FWB
#' covariance is still used for confidence intervals, but the bias correction
#' uses the same Fisher-style model variance as in the MLE case, evaluated at
#' the non-MLE estimate. For iOLS/GPML, the model-based variance is
#' \eqn{(\mathbf{Z}'\mathbf{Z})^{-1}} (Gamma Fisher information, no
#' dispersion scaling).
#'
#' When \code{constrained = TRUE}, the delta method accounts for the
#' logit link: \eqn{\mathrm{Var}(\alpha) = \mathrm{Var}(\eta) \cdot
#' [\sigma'(\eta)]^2} where \eqn{\sigma'(\eta) = \alpha(1-\alpha)}.
#'
#' **Confidence intervals via the delta method.** Let
#' \eqn{\mathbf{g}_g = \partial \hat{\xi}_g / \partial \boldsymbol{\alpha}}
#' denote the gradient of the plug-in estimator with respect to the alpha
#' coefficients. For unconstrained models,
#' \eqn{\mathbf{g}_g = \sum_i N_i^{\hat{\alpha}_i} \log(N_i)\mathbf{x}_i;}
#' for constrained models the same expression is multiplied by the derivative
#' of the inverse-logit map. Using the HC-robust covariance
#' \eqn{\mathbf{V}}, the package computes
#' \deqn{
#'   \widehat{\mathrm{Var}}(\hat{\xi}_g) =
#'     \mathbf{g}_g^\top \mathbf{V}\mathbf{g}_g.
#' }
#' To preserve positivity, the subgroup interval is then reported on a
#' log-normal scale:
#' \deqn{
#'   \hat{\xi}_g \exp\!\left(
#'     \pm z_{\alpha/2}
#'     \frac{\widehat{\mathrm{se}}(\hat{\xi}_g)}{\hat{\xi}_g}
#'   \right).
#' }
#' When \code{bias_correction = TRUE}, the lower and upper bounds are
#' rescaled by \eqn{\hat{\xi}^{BC}_g / \hat{\xi}_g}, matching the returned
#' bias-corrected point estimate.
#'
#' **Total across groups.** When multiple groups exist, the total
#' \eqn{\hat{\xi} = \sum_g \hat{\xi}_g} has its own delta-method CI
#' computed via the gradient \eqn{\nabla_\alpha \xi} and a log-normal
#' approximation for positivity. This is stored in \code{attr(result, "total")}.
#'
#' @param object An `"uncounted"` object.
#' @param level Confidence level for intervals (default 0.95).
#' @param bias_correction Logical; apply analytical bias correction?
#'   Default TRUE. Uses model-based variance (not HC-robust) to avoid
#'   overcorrection from inflated leverage-driven standard errors.
#' @param total Logical; if \code{TRUE} and multiple groups exist, compute a
#'   delta-method total with SE and CI, stored in \code{attr(result, "total")}.
#'   Default \code{FALSE}. \strong{Warning}: for panel data where groups are
#'   defined by time periods (e.g., \code{by = ~ year}), the total sums
#'   population estimates across years. This is generally not meaningful because
#'   the same individuals may appear in multiple years. The total is only
#'   appropriate when groups represent non-overlapping subpopulations
#'   (e.g., \code{by = ~ sex} within a single year).
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{group}{Group label derived from alpha covariates, or \code{"(all)"}
#'     when no alpha covariates are specified.}
#'   \item{estimate}{Plug-in estimate \eqn{\hat{\xi}_g = \sum N_i^{\hat{\alpha}_g}}.}
#'   \item{estimate_bc}{Bias-corrected estimate \eqn{\hat{\xi}^{BC}_g}; \code{NA} when \code{bias_correction = FALSE}.}
#'   \item{lower}{Lower bound of the (bias-corrected) confidence interval.}
#'   \item{upper}{Upper bound of the (bias-corrected) confidence interval.}
#'   \item{share_pct}{Group share as percentage of total \eqn{\hat{\xi}}.}
#' }
#' When multiple groups exist, \code{attr(result, "total")} contains
#' the total estimate, bias-corrected estimate, standard error, and CI.
#'
#' @examples
#' # Simulate synthetic data for 5 countries, 3 years each
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
#' # Fit a Poisson model
#' fit <- estimate_hidden_pop(
#'   data = sim_data, observed = ~m, auxiliary = ~n,
#'   reference_pop = ~N, method = "poisson"
#' )
#'
#' # Population size with bias correction and 95% CI
#' popsize(fit)
#'
#' # Without bias correction
#' popsize(fit, bias_correction = FALSE)
#'
#' # 90% confidence interval
#' popsize(fit, level = 0.90)
#'
#' @export
popsize <- function(object, ...) {
  UseMethod("popsize")
}

#' @rdname popsize
#' @param by Optional formula specifying grouping variables for stratified
#'   population size estimation (e.g., \code{~ year}, \code{~ country},
#'   \code{~ year + sex}). The variables must exist in the data used to fit
#'   the model. When \code{NULL} (default), groups are defined by the alpha
#'   covariate pattern from \code{cov_alpha}. When provided, the same fitted
#'   alpha values are used but summed over the \code{by}-groups instead.
#'   This is analogous to \code{singleRcapture::stratifyPopsize()}.
#'
#' @export
popsize.uncounted <- function(object, by = NULL, level = 0.95,
                               bias_correction = TRUE, total = FALSE, ...) {
  alpha_coefs <- object$alpha_coefs
  V <- object$vcov
  V_model <- object$vcov_model
  N <- object$N
  X_alpha <- object$X_alpha
  p_alpha <- object$p_alpha
  z_crit <- qnorm(1 - (1 - level) / 2)
  is_constr <- isTRUE(object$constrained)

  # Robust vcov for CI (HC3 or user-selected)
  V_alpha <- V[seq_len(p_alpha), seq_len(p_alpha), drop = FALSE]

  # Model-based vcov for bias correction (homoscedastic)
  V_alpha_model <- if (!is.null(V_model)) {
    V_model[seq_len(p_alpha), seq_len(p_alpha), drop = FALSE]
  } else {
    V_alpha
  }

  # Compute alpha_i for every observation
  alpha_all_eta <- as.numeric(X_alpha %*% alpha_coefs)
  alpha_all <- if (is_constr) .inv_logit(alpha_all_eta) else alpha_all_eta

  # Define groups
  group_info <- .popsize_group_info(object, by = by)
  group_idx <- group_info$index

  # Compute per-group estimates
  results <- vector("list", length(group_idx))
  gradients <- matrix(NA_real_, nrow = length(group_idx), ncol = p_alpha,
                      dimnames = list(names(group_idx), colnames(X_alpha)))
  bc_scale <- rep(NA_real_, length(group_idx))

  for (k in seq_along(group_idx)) {
    idx <- group_idx[[k]]
    label <- names(group_idx)[k]
    N_g <- N[idx]
    alpha_g <- alpha_all[idx]
    m_g <- sum(object$m[idx])

    # Point estimate: sum(N_i^alpha_i)
    est <- sum(N_g^alpha_g)

    # --- Delta-method CI ---
    # Gradient w.r.t. alpha coefficients: sum_i N_i^alpha_i * log(N_i) * x_i
    if (is_constr) {
      dsig <- alpha_g * (1 - alpha_g)
      g_vec <- colSums(X_alpha[idx, , drop = FALSE] *
                         as.numeric(N_g^alpha_g * log(N_g) * dsig))
    } else {
      g_vec <- colSums(X_alpha[idx, , drop = FALSE] *
                         as.numeric(N_g^alpha_g * log(N_g)))
    }

    var_xi <- as.numeric(t(g_vec) %*% V_alpha %*% g_vec)
    se_xi <- sqrt(max(var_xi, 0))
    gradients[k, ] <- g_vec

    # Log-normal CI for positivity
    if (isTRUE(se_xi > 0) && isTRUE(est > 0)) {
      log_se <- se_xi / est
      ci_lower <- est * exp(-z_crit * log_se)
      ci_upper <- est * exp(z_crit * log_se)
    } else {
      ci_lower <- ci_upper <- NA_real_
    }

    # --- Bias correction ---
    est_bc <- NA_real_
    ci_lower_bc <- ci_lower
    ci_upper_bc <- ci_upper
    if (bias_correction) {
      log_N2 <- (log(N_g))^2
      x_sub <- X_alpha[idx, , drop = FALSE]
      xVx <- rowSums((x_sub %*% V_alpha_model) * x_sub)

      if (is_constr) {
        # Constrained (logit link): subtractive correction with sigma' and sigma''
        # NOTE: h''(eta) is NOT always positive under logit — bias can go either
        # direction. Near alpha=1, the correction can INCREASE the estimate.
        dsig <- alpha_g * (1 - alpha_g)
        dsig2 <- dsig^2
        dsig_dd <- dsig * (1 - 2 * alpha_g)
        bias <- 0.5 * sum(N_g^alpha_g * (log_N2 * dsig2 + log(N_g) * dsig_dd) * xVx)
        est_bc <- est - bias
        # Clamp only to positive (bias can go either direction for constrained)
        if (est_bc <= 0) est_bc <- est
        if (isTRUE(est > 0) && isTRUE(est_bc > 0)) {
          bc_ratio <- est_bc / est
          bc_scale[k] <- bc_ratio
          if (!is.na(ci_lower)) {
            ci_lower_bc <- ci_lower * bc_ratio
            ci_upper_bc <- ci_upper * bc_ratio
          }
        }
      } else {
        # Unconstrained: multiplicative lognormal correction (exact under normality)
        # E[N^alpha_hat] = N^alpha_0 * exp(0.5 * (log N)^2 * x'Vx)
        # So xi_BC = sum N^alpha_hat * exp(-0.5 * (log N)^2 * x'Vx)
        # Always positive. More accurate than subtractive Taylor approx.
        correction <- exp(-0.5 * log_N2 * xVx)
        est_bc <- sum(N_g^alpha_g * correction)
        if (isTRUE(est > 0)) {
          bc_ratio <- est_bc / est
          bc_scale[k] <- bc_ratio
          if (!is.na(ci_lower)) {
            ci_lower_bc <- ci_lower * bc_ratio
            ci_upper_bc <- ci_upper * bc_ratio
          }
        }
      }
    }

    results[[k]] <- data.frame(
      group = label,
      observed = m_g,
      estimate = est,
      estimate_bc = est_bc,
      lower = ci_lower_bc,
      upper = ci_upper_bc,
      stringsAsFactors = FALSE
    )
  }

  ps <- do.call(rbind, results)
  rownames(ps) <- NULL
  ps$share_pct <- ps$estimate / sum(ps$estimate) * 100

  vcov_estimate <- gradients %*% V_alpha %*% t(gradients)
  rownames(vcov_estimate) <- colnames(vcov_estimate) <- ps$group
  vcov_estimate_bc <- vcov_estimate
  if (bias_correction && all(is.finite(bc_scale))) {
    vcov_estimate_bc <- vcov_estimate * tcrossprod(bc_scale, bc_scale)
  } else {
    vcov_estimate_bc[,] <- NA_real_
  }

  # Total with delta-method CI
  if (total && nrow(ps) > 1) {
    xi_total <- sum(ps$estimate)
    xi_total_bc <- if (bias_correction) sum(ps$estimate_bc, na.rm = TRUE) else NA_real_

    # Total gradient over all observations
    if (is_constr) {
      dsig_all <- alpha_all * (1 - alpha_all)
      g_total <- colSums(X_alpha * as.numeric(N^alpha_all * log(N) * dsig_all))
    } else {
      g_total <- colSums(X_alpha * as.numeric(N^alpha_all * log(N)))
    }
    var_xi_total <- as.numeric(t(g_total) %*% V_alpha %*% g_total)
    se_xi_total <- sqrt(max(var_xi_total, 0))

    if (isTRUE(se_xi_total > 0) && isTRUE(xi_total > 0)) {
      log_se <- se_xi_total / xi_total
      total_lower <- xi_total * exp(-z_crit * log_se)
      total_upper <- xi_total * exp(z_crit * log_se)
    } else {
      total_lower <- total_upper <- NA_real_
    }

    attr(ps, "total") <- list(
      estimate = xi_total,
      estimate_bc = xi_total_bc,
      se = se_xi_total,
      lower = total_lower,
      upper = total_upper
    )
  }

  attr(ps, "groups") <- group_info$groups
  attr(ps, "vcov_estimate") <- vcov_estimate
  attr(ps, "vcov_estimate_bc") <- vcov_estimate_bc
  class(ps) <- c("uncounted_popsize", "data.frame")
  ps
}


#' @export
print.uncounted_popsize <- function(x, ...) {
  # Print the data frame
  print.data.frame(x, ...)
  # Print total if present
  tot <- attr(x, "total")
  if (!is.null(tot)) {
    cat("\nTotal: estimate =", format(round(tot$estimate), big.mark = ","),
        " BC =", format(round(tot$estimate_bc), big.mark = ","),
        " [", format(round(tot$lower), big.mark = ","),
        ",", format(round(tot$upper), big.mark = ","), "]\n")
  }
  invisible(x)
}


#' Plot Population Size Estimates
#'
#' Visualizes population size estimates by group with confidence intervals.
#'
#' @param x An \code{"uncounted_popsize"} object from \code{\link{popsize}}.
#' @param type \code{"estimate"} (default) shows bias-corrected estimates with
#'   CIs; \code{"compare"} shows plug-in and bias-corrected side by side.
#' @param ... Additional graphical arguments passed to \code{\link{plot}}.
#'
#' @export
plot.uncounted_popsize <- function(x, type = c("estimate", "compare"), ...) {
  type <- match.arg(type)
  ng <- nrow(x)
  labs <- x$group
  x_pos <- seq_len(ng)

  if (type == "compare") {
    ylim <- range(c(x$lower, x$upper, x$estimate, x$estimate_bc), na.rm = TRUE)
    ylim[1] <- max(ylim[1], 0)
    offset <- 0.15

    plot(NULL, xlim = c(0.5, ng + 0.5), ylim = ylim,
         xaxt = "n", xlab = "", ylab = "Population size",
         main = "Plug-in vs bias-corrected estimates", ...)
    axis(1, at = x_pos, labels = labs, las = 2, cex.axis = 0.8)
    abline(v = x_pos, col = "gray90")

    points(x_pos - offset, x$estimate, pch = 16, col = "steelblue")
    points(x_pos + offset, x$estimate_bc, pch = 17, col = "firebrick")
    segments(x_pos + offset, x$lower, x_pos + offset, x$upper,
             col = "firebrick", lwd = 1.5)

    legend("topright", legend = c("Plug-in", "Bias-corrected"),
           pch = c(16, 17), col = c("steelblue", "firebrick"),
           bty = "n", cex = 0.9)
  } else {
    ylim <- range(c(x$lower, x$upper), na.rm = TRUE)
    ylim[1] <- max(ylim[1], 0)

    plot(x_pos, x$estimate_bc, ylim = ylim, xlim = c(0.5, ng + 0.5),
         pch = 16, xaxt = "n", xlab = "", ylab = "Population size",
         main = "Population size estimates (bias-corrected)", ...)
    axis(1, at = x_pos, labels = labs, las = 2, cex.axis = 0.8)
    segments(x_pos, x$lower, x_pos, x$upper, lwd = 1.5)
    abline(v = x_pos, col = "gray90", lty = 3)
    points(x_pos, x$estimate_bc, pch = 16)
  }
  invisible(x)
}


#' Compare Population Size Across Models
#'
#' Produces a side-by-side comparison of population size estimates from
#' multiple fitted models.
#'
#' @param ... Fitted \code{uncounted} objects to compare.
#' @param labels Character vector of model labels. Defaults to model method.
#' @param by Optional formula for grouping (passed to \code{\link{popsize}}).
#' @param bias_correction Logical; apply bias correction? Default TRUE.
#' @return An object of class \code{"uncounted_popsize_compare"} with
#'   components \code{table} (data frame) and \code{labels} (model labels).
#'
#' @export
compare_popsize <- function(..., labels = NULL, by = NULL,
                            bias_correction = TRUE) {
  models <- list(...)
  n_models <- length(models)
  if (n_models < 2) stop("Need at least 2 models to compare.")

  if (is.null(labels)) {
    labels <- vapply(models, function(m) toupper(m$method), character(1))
  }

  tabs <- lapply(seq_len(n_models), function(i) {
    ps <- popsize(models[[i]], by = by, bias_correction = bias_correction)
    ps$model <- labels[i]
    ps
  })

  combined <- do.call(rbind, tabs)
  rownames(combined) <- NULL

  result <- list(table = combined, labels = labels, n_models = n_models)
  class(result) <- "uncounted_popsize_compare"
  result
}


#' @export
print.uncounted_popsize_compare <- function(x, ...) {
  cat("Population size comparison:", paste(x$labels, collapse = " vs "), "\n\n")
  print(x$table[, c("model", "group", "estimate", "estimate_bc", "lower", "upper")])
  invisible(x)
}


#' @export
plot.uncounted_popsize_compare <- function(x, ...) {
  tab <- x$table
  labels <- x$labels
  groups <- unique(tab$group)
  ng <- length(groups)
  nm <- x$n_models

  cols <- c("steelblue", "firebrick", "forestgreen", "darkorange",
            "purple", "brown")[seq_len(nm)]

  ylim <- range(c(tab$lower, tab$upper), na.rm = TRUE)
  ylim[1] <- max(ylim[1], 0)
  x_pos <- seq_len(ng)
  spacing <- 0.8 / nm

  plot(NULL, xlim = c(0.5, ng + 0.5), ylim = ylim,
       xaxt = "n", xlab = "", ylab = "Population size",
       main = "Population size comparison", ...)
  axis(1, at = x_pos, labels = groups, las = 2, cex.axis = 0.8)
  abline(v = x_pos, col = "gray90", lty = 3)

  for (j in seq_len(nm)) {
    sub <- tab[tab$model == labels[j], ]
    offset <- (j - (nm + 1) / 2) * spacing
    x_j <- match(sub$group, groups) + offset
    points(x_j, sub$estimate_bc, pch = 15 + j, col = cols[j])
    segments(x_j, sub$lower, x_j, sub$upper, col = cols[j], lwd = 1.5)
  }

  legend("topright", legend = labels, pch = 15 + seq_len(nm),
         col = cols, bty = "n", cex = 0.9)
  invisible(x)
}


#' @rdname popsize
#' @usage NULL
#' @export
xi <- function(object, ...) {
  .Deprecated("popsize")
  UseMethod("popsize")
}

#' Build a readable group label from covariate data
#' @noRd
.make_group_label <- function(row_idx, cov_data) {
  if (is.null(cov_data) || ncol(cov_data) == 0) return("(all)")
  vals <- cov_data[row_idx, , drop = FALSE]
  parts <- vapply(seq_len(ncol(vals)), function(j) {
    v <- vals[[j]]
    nm <- names(vals)[j]
    paste0(nm, "=", as.character(v))
  }, character(1))
  paste(parts, collapse = ", ")
}

#' Resolve population-size grouping indices
#' @noRd
.popsize_group_index <- function(object, by = NULL) {
  .popsize_group_info(object, by = by)$index
}

#' Resolve population-size grouping indices and metadata
#' @noRd
.popsize_group_info <- function(object, by = NULL) {
  X_alpha <- object$X_alpha

  if (!is.null(by)) {
    by_vars <- all.vars(by)
    if (length(by_vars) == 0L) {
      out <- list(seq_len(nrow(object$data)))
      names(out) <- "(all)"
      groups <- data.frame(.group = "(all)", group = "(all)",
                           stringsAsFactors = FALSE)
      return(list(index = out, groups = groups))
    }
    if (!all(by_vars %in% names(object$data))) {
      missing <- setdiff(by_vars, names(object$data))
      stop("Variables not found in data: ", paste(missing, collapse = ", "))
    }
    by_data <- object$data[, by_vars, drop = FALSE]
    group_factor <- interaction(by_data, drop = TRUE, sep = ", ")
    group_labels <- levels(group_factor)
    out <- lapply(group_labels, function(lev) which(group_factor == lev))
    names(out) <- group_labels
    groups <- do.call(rbind, lapply(seq_along(out), function(i) {
      row <- by_data[out[[i]][1], , drop = FALSE]
      row$.group <- group_labels[i]
      row
    }))
    rownames(groups) <- NULL
    groups <- groups[, c(".group", by_vars), drop = FALSE]
    if (!("group" %in% names(groups))) {
      groups$group <- groups$.group
    }
    return(list(index = out, groups = groups))
  }

  X_key <- apply(X_alpha, 1, paste, collapse = "|")
  unique_keys <- unique(X_key)
  out <- lapply(unique_keys, function(key) which(X_key == key))
  group_labels <- vapply(out, function(idx) {
    .make_group_label(idx[1], object$cov_alpha_vars)
  }, character(1))
  names(out) <- group_labels
  if (!is.null(object$cov_alpha_vars) && ncol(object$cov_alpha_vars) > 0) {
    groups <- do.call(rbind, lapply(seq_along(out), function(i) {
      row <- object$cov_alpha_vars[out[[i]][1], , drop = FALSE]
      row$.group <- group_labels[i]
      row
    }))
    rownames(groups) <- NULL
    groups <- groups[, c(".group", names(object$cov_alpha_vars)), drop = FALSE]
  } else {
    groups <- data.frame(.group = group_labels, stringsAsFactors = FALSE)
  }
  if (!("group" %in% names(groups))) {
    groups$group <- groups$.group
  }
  list(index = out, groups = groups)
}
