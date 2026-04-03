#' Estimated Population Size
#'
#' Computes the estimated total unauthorized population
#' \eqn{\hat{\xi} = \sum_{i} N_i^{\hat{\alpha}}} for each group
#' defined by the covariates in alpha. Includes bias correction via
#' second-order Taylor expansion and confidence intervals via monotone
#' transformation of the Wald interval on the link scale.
#'
#' @details
#' **Point estimate.** For a group \eqn{g} with reference populations
#' \eqn{N_1, \ldots, N_n} and estimated exponent \eqn{\hat{\alpha}_g},
#' the plug-in population size is
#' \deqn{\hat{\xi}_g = \sum_{i=1}^{n} N_i^{\hat{\alpha}_g}.}
#'
#' **Bias correction via Taylor expansion.** Because \eqn{\xi} is a nonlinear
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
#' The bias-corrected estimate is
#' \eqn{\hat{\xi}^{BC}_g = \hat{\xi}_g - \widehat{\mathrm{Bias}}}.
#' Model-based (homoscedastic) variance is used for bias correction rather
#' than HC-robust variance, because HC3 can be inflated by high-leverage
#' observations in skewed data, leading to overcorrection.
#'
#' When \code{constrained = TRUE}, the delta method accounts for the
#' logit link: \eqn{\mathrm{Var}(\alpha) = \mathrm{Var}(\eta) \cdot
#' [\sigma'(\eta)]^2} where \eqn{\sigma'(\eta) = \alpha(1-\alpha)}.
#'
#' **Confidence intervals via monotone transformation.** A Wald interval is
#' first constructed on the link scale for the linear predictor:
#' \deqn{\hat{\eta}_g \pm z_{\alpha/2} \cdot \mathrm{se}(\hat{\eta}_g),}
#' where the standard error uses the HC-robust variance.
#' The interval endpoints are then mapped through the monotone transformation
#' \eqn{g(\alpha) = \sum_i N_i^{\alpha}} (increasing for \eqn{N_i \geq 1})
#' to obtain \eqn{[\hat{\xi}_L, \hat{\xi}_U]}. When \code{constrained = TRUE},
#' the logit link is applied before exponentiation. Bias correction
#' is also applied to the CI bounds at their respective alpha values.
#'
#' **Total across groups.** When multiple groups exist, the total
#' \eqn{\hat{\xi} = \sum_g \hat{\xi}_g} has its own delta-method CI
#' computed via the gradient \eqn{\nabla_\alpha \xi} and a log-normal
#' approximation for positivity. This is stored in \code{attr(result, "total")}.
#'
#' @param object An `"uncounted"` object.
#' @param level Confidence level for intervals (default 0.95).
#' @param bias_correction Logical; apply Taylor-expansion bias correction?
#'   Default TRUE. Uses model-based variance (not HC-robust) to avoid
#'   overcorrection from inflated leverage-driven standard errors.
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{group}{Group label derived from alpha covariates, or \code{"(all)"}
#'     when no alpha covariates are specified.}
#'   \item{estimate}{Plug-in estimate \eqn{\hat{\xi}_g = \sum N_i^{\hat{\alpha}_g}}.}
#'   \item{estimate_bc}{Bias-corrected estimate \eqn{\hat{\xi}^{BC}_g}.}
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
  if (!is.null(by)) {
    # Stratified by user-specified variables
    by_vars <- all.vars(by)
    if (!all(by_vars %in% names(object$data))) {
      missing <- setdiff(by_vars, names(object$data))
      stop("Variables not found in data: ", paste(missing, collapse = ", "))
    }
    by_data <- object$data[, by_vars, drop = FALSE]
    group_factor <- interaction(by_data, drop = TRUE, sep = ", ")
    group_labels <- levels(group_factor)
    group_idx <- lapply(group_labels, function(lev) which(group_factor == lev))
    names(group_idx) <- group_labels
  } else {
    # Default: group by unique cov_alpha patterns
    X_key <- apply(X_alpha, 1, paste, collapse = "|")
    unique_keys <- unique(X_key)
    group_idx <- lapply(unique_keys, function(key) which(X_key == key))
    # Labels from covariate values
    group_labels <- vapply(group_idx, function(idx) {
      .make_group_label(idx[1], object$cov_alpha_vars)
    }, character(1))
    names(group_idx) <- group_labels
  }

  # Compute per-group estimates
  results <- vector("list", length(group_idx))

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

    # Log-normal CI for positivity
    if (se_xi > 0 && est > 0) {
      log_se <- se_xi / est
      ci_lower <- est * exp(-z_crit * log_se)
      ci_upper <- est * exp(z_crit * log_se)
    } else {
      ci_lower <- ci_upper <- NA_real_
    }

    # --- Bias correction ---
    est_bc <- est
    ci_lower_bc <- ci_lower
    ci_upper_bc <- ci_upper
    if (bias_correction) {
      # Per-observation bias: 0.5 * N_i^alpha_i * (log N_i)^2 * x_i' V_model x_i
      log_N2 <- (log(N_g))^2
      x_sub <- X_alpha[idx, , drop = FALSE]
      # x_i' V x_i for each obs
      xVx <- rowSums((x_sub %*% V_alpha_model) * x_sub)
      if (is_constr) {
        dsig2 <- (alpha_g * (1 - alpha_g))^2
        xVx <- xVx * dsig2
      }
      bias <- 0.5 * sum(N_g^alpha_g * log_N2 * xVx)
      est_bc <- est - bias
      # Also correct CI bounds (approximate: same relative bias)
      if (!is.na(ci_lower) && est > 0) {
        ci_lower_bc <- ci_lower - bias * (ci_lower / est)
        ci_upper_bc <- ci_upper - bias * (ci_upper / est)
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

  # Total with delta-method CI
  if (total && nrow(ps) > 1) {
    xi_total <- sum(ps$estimate)
    xi_total_bc <- sum(ps$estimate_bc)

    # Total gradient over all observations
    if (is_constr) {
      dsig_all <- alpha_all * (1 - alpha_all)
      g_total <- colSums(X_alpha * as.numeric(N^alpha_all * log(N) * dsig_all))
    } else {
      g_total <- colSums(X_alpha * as.numeric(N^alpha_all * log(N)))
    }
    var_xi_total <- as.numeric(t(g_total) %*% V_alpha %*% g_total)
    se_xi_total <- sqrt(max(var_xi_total, 0))

    if (se_xi_total > 0 && xi_total > 0) {
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

  ps
}

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
