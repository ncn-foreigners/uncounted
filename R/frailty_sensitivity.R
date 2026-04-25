#' Omitted-Frailty Sensitivity Analysis for Hidden Population Size
#'
#' Computes a `sensemakr`-style omitted-frailty sensitivity analysis for the
#' hidden-population estimand \eqn{\Xi}. The method linearizes the fitted
#' Poisson or Negative Binomial mean model with one IRLS step, treats the
#' omitted shared frailty as an unobserved weighted confounder, and reports
#' how large the resulting first-order bias would have to be to change the
#' plug-in hidden-population estimate.
#'
#' @param object An \code{"uncounted"} object.
#' @param by Optional one-sided formula defining the \eqn{\Xi} targets. This
#'   uses the same grouping logic as \code{\link{popsize}}.
#' @param r2_d Numeric vector of weighted partial-\eqn{R^2} values describing
#'   how strongly the omitted frailty aligns with the targeted alpha contrast.
#'   Values must satisfy \code{0 <= r2_d < 1}. Duplicates are removed and the
#'   grid is sorted increasingly.
#' @param r2_y Numeric vector of weighted partial-\eqn{R^2} values describing
#'   how strongly the omitted frailty predicts the IRLS working outcome after
#'   conditioning on the fitted model. Values must satisfy
#'   \code{0 <= r2_y <= 1}. Duplicates are removed and the grid is sorted
#'   increasingly.
#' @param q Numeric vector of relative target changes used when
#'   \code{threshold = NULL}. For \code{direction = "decrease"}, each target is
#'   \code{Xi_hat * (1 - q)} and values must lie in \code{[0, 1]}. For
#'   \code{direction = "increase"}, the targets are
#'   \code{Xi_hat * (1 + q)} and values must be non-negative.
#' @param threshold Optional positive numeric threshold for \eqn{\Xi}. When
#'   supplied, \code{q} is ignored for target construction and the same
#'   threshold is applied to every reported group.
#' @param direction Which target crossing to summarize:
#'   \code{"decrease"} or \code{"increase"}.
#' @param df_method Degrees-of-freedom rule for the omitted-frailty formulas.
#'   \code{"model"} uses \code{n_obs - rank(X_work)}. \code{"cluster"} uses
#'   the number of observed clusters minus the same rank when cluster metadata
#'   is present on the fitted object. \code{"auto"} prefers cluster degrees of
#'   freedom when available and valid, otherwise falls back to the model-based
#'   rule.
#' @param plot Logical; if \code{TRUE} (default), draw a contour plot for the
#'   first reported group.
#' @param x Object to print or plot.
#' @param group Optional group label for the plot method. When \code{NULL},
#'   the first available group is used.
#' @param type Plot scale: \code{"ratio"} for values relative to baseline
#'   \eqn{\Xi}, or \code{"xi"} for raw \eqn{\Xi} values.
#' @param side Which sensitivity surface to draw: \code{"lower"} or
#'   \code{"upper"}.
#' @param ... Additional graphical arguments passed to
#'   \code{\link{plot.uncounted_frailty_sensitivity}} when
#'   \code{plot = TRUE}.
#'
#' @details
#' The baseline \pkg{uncounted} mean model factorizes the observed-count mean as
#'
#' \deqn{\mu_i = E(m_i \mid N_i, n_i) = \xi_i \rho_i,}
#'
#' where \eqn{\xi_i} is the latent population component and \eqn{\rho_i} is the
#' detection component. Under the default power-link specification,
#'
#' \deqn{
#' \log \mu_i = \alpha_i \log N_i + \beta_i \log(\gamma + n_i / N_i).
#' }
#'
#' This helper studies violations of the identifying assumption by introducing
#' an omitted shared frailty \eqn{U_i} into the log-mean:
#'
#' \deqn{
#' \log \mu_i^\dagger =
#'   \alpha_i \log N_i +
#'   \beta_i \log(\gamma + n_i / N_i) +
#'   \delta U_i.
#' }
#'
#' The frailty is not observed or estimated directly. Instead, its strength is
#' indexed by two weighted partial-\eqn{R^2} quantities:
#'
#' \describe{
#'   \item{\code{r2_d}}{The weighted partial association between the omitted
#'     frailty and the targeted alpha contrast after conditioning on the fitted
#'     nuisance regressors.}
#'   \item{\code{r2_y}}{The weighted partial association between the omitted
#'     frailty and the working IRLS outcome after conditioning on the targeted
#'     contrast and the same nuisance regressors.}
#' }
#'
#' For a fitted Poisson or NB model, the method linearizes the mean model with
#' the working response
#'
#' \deqn{
#' Y_i^* = \hat\eta_i + \frac{m_i - \hat\mu_i}{\hat\mu_i},
#' }
#'
#' using working weights
#'
#' \deqn{
#' w_i = \hat\mu_i \qquad \text{(Poisson)}
#' }
#'
#' and
#'
#' \deqn{
#' w_i = \frac{\hat\mu_i}{1 + \hat\mu_i / \hat\theta} \qquad \text{(NB2)}.
#' }
#'
#' The alpha block of the weighted linearized design is
#' \eqn{X_{\alpha,i} \log N_i}; the beta block is
#' \eqn{X_{\beta,i} \log(\gamma + n_i/N_i)}. When scalar gamma is estimated,
#' the nuisance design also includes the derivative column
#'
#' \deqn{
#' \frac{\partial \eta_i}{\partial \gamma} =
#' \frac{\hat\beta_i}{\hat\gamma + n_i / N_i}.
#' }
#'
#' The primary estimand is not a raw alpha coefficient but the grouped
#' hidden-population estimate
#'
#' \deqn{
#' \hat\Xi_g = \sum_{i \in g} N_i^{\hat\alpha_i}.
#' }
#'
#' For each requested group, the helper computes the gradient of
#' \eqn{\hat\Xi_g} with respect to the alpha-coefficient vector and uses that
#' gradient to define a one-dimensional alpha contrast. This makes the omitted
#' frailty analysis applicable to vector-valued \code{cov_alpha} models such as
#' \code{~ year * ukr + sex}. The reported contrast is an internal device: it
#' aligns the sensitivity calculation with the direction in coefficient space
#' that matters most for the chosen \eqn{\Xi_g}.
#'
#' Let \eqn{\theta_g} denote the targeted alpha contrast, with estimated
#' standard error \eqn{\widehat{\mathrm{se}}(\hat\theta_g)} and degrees of
#' freedom \eqn{\nu}. The first-order omitted-frailty bias is
#'
#' \deqn{
#' B_{\theta,g}(r2_d, r2_y) =
#' \widehat{\mathrm{se}}(\hat\theta_g)
#' \sqrt{
#'   \nu \frac{r2_y \, r2_d}{1 - r2_d}
#' }.
#' }
#'
#' The hidden-population effect is then approximated by a delta method:
#'
#' \deqn{
#' \hat\Xi_{g,\mathrm{lower}} \approx
#' \hat\Xi_g - \frac{\partial \Xi_g}{\partial \theta_g}
#' B_{\theta,g}(r2_d, r2_y),
#' }
#'
#' \deqn{
#' \hat\Xi_{g,\mathrm{upper}} \approx
#' \hat\Xi_g + \frac{\partial \Xi_g}{\partial \theta_g}
#' B_{\theta,g}(r2_d, r2_y).
#' }
#'
#' The returned \code{surface} table stores these lower/upper approximations as
#' both raw \eqn{\Xi} values and ratios relative to the baseline estimate.
#'
#' The \code{robustness} table summarizes two robustness values for threshold
#' questions:
#'
#' \describe{
#'   \item{\code{rv_equal}}{The minimum equal-strength partial-\eqn{R^2}
#'     required for the omitted frailty to move the targeted \eqn{\Xi} to the
#'     requested threshold under \eqn{r2_d = r2_y}.}
#'   \item{\code{rv_extreme}}{The corresponding treatment-side strength when the
#'     omitted frailty is given the extreme outcome-side scenario
#'     \eqn{r2_y = 1}.}
#' }
#'
#' These robustness values are derived from the same first-order bias mapping.
#' They quantify how much residual alignment between the omitted frailty and the
#' targeted alpha contrast would be required to overturn a substantive claim.
#'
#' \strong{Interpretation.} This is an identification-sensitivity analysis,
#' distinct from the package's existing distributional sensitivity
#' (\code{method = "nb"}) and bootstrap uncertainty summaries. It is inspired by
#' omitted-variable-bias diagnostics such as \pkg{sensemakr}, but it is not a
#' literal OLS port: the count model is first linearized with its IRLS working
#' representation, and the hidden-population estimand is handled through a
#' grouped delta-method contrast.
#'
#' \strong{Limitations.} The method is first-order and deliberately restricted
#' to \code{link_rho = "power"}, unconstrained Poisson/NB MLE fits without
#' \code{cov_gamma}. It does not yet implement benchmark covariates or a
#' simulation-based calibration layer.
#'
#' @return An object of class \code{"uncounted_frailty_sensitivity"} with:
#' \describe{
#'   \item{\code{baseline}}{A \code{\link{popsize}} table computed with
#'     \code{bias_correction = FALSE} and the requested \code{by =} grouping.}
#'   \item{\code{working}}{Data frame with one row per reported target and
#'     columns including \code{group}, \code{xi_hat}, \code{theta_hat},
#'     \code{se_theta}, \code{d_xi_dtheta}, and \code{df}.}
#'   \item{\code{surface}}{Tidy grid over \code{group}, \code{r2_d}, and
#'     \code{r2_y} containing \code{bias_theta}, \code{xi_lower},
#'     \code{xi_upper}, \code{xi_ratio_lower}, and \code{xi_ratio_upper}.}
#'   \item{\code{robustness}}{Data frame with robustness values for each
#'     requested target.}
#'   \item{\code{settings}}{List of resolved inputs and defaults.}
#' }
#'
#' @references
#' Cinelli, C., & Hazlett, C. (2020). Making sense of sensitivity:
#' Extending omitted variable bias. \emph{Journal of the Royal Statistical
#' Society: Series B}, 82(1), 39--67.
#'
#' Beręsewicz, M., & Pawlukiewicz, K. (2020). Estimation of the number of
#' irregular foreigners in Poland using non-linear count regression models.
#' \emph{arXiv preprint} arXiv:2008.09407.
#'
#' Zhang, L.-C. (2008). Developing methods for determining the number of
#' unauthorized foreigners in Norway. \emph{Documents} 2008/11, Statistics
#' Norway.
#'
#' @examples
#' data(irregular_migration)
#' d <- irregular_migration[
#'   irregular_migration$m > 0 & irregular_migration$n > 0,
#' ]
#' keep <- unique(d$country)[1:8]
#' d <- droplevels(d[d$country %in% keep, ])
#'
#' fit <- estimate_hidden_pop(
#'   data = d,
#'   observed = ~m,
#'   auxiliary = ~n,
#'   reference_pop = ~N,
#'   method = "poisson",
#'   cov_alpha = ~ sex,
#'   gamma = 0.005
#' )
#'
#' frailty_sensitivity(
#'   fit,
#'   by = ~ year,
#'   r2_d = c(0, 0.05, 0.10),
#'   r2_y = c(0, 0.05, 0.10),
#'   plot = FALSE
#' )
#'
#' @seealso \code{\link{dependence_bounds}}, \code{\link{profile_dependence}},
#'   \code{\link{robustness_dependence}}, \code{\link{popsize}}
#'
#' @export
frailty_sensitivity <- function(object, by = NULL,
                                r2_d = seq(0, 0.25, by = 0.01),
                                r2_y = seq(0, 0.25, by = 0.01),
                                q = c(0.10, 0.25, 0.50),
                                threshold = NULL,
                                direction = c("decrease", "increase"),
                                df_method = c("auto", "model", "cluster"),
                                plot = TRUE, ...) {
  direction <- match.arg(direction)
  df_method <- match.arg(df_method)

  .validate_frailty_object(object)
  r2_d <- .validate_r2_grid(r2_d, lower = 0, upper = 1, inclusive_upper = FALSE,
                            arg = "r2_d")
  r2_y <- .validate_r2_grid(r2_y, lower = 0, upper = 1, inclusive_upper = TRUE,
                            arg = "r2_y")
  q <- .validate_q_grid(q, direction = direction)

  if (!is.null(threshold) &&
      (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) ||
       !is.finite(threshold) || threshold <= 0)) {
    stop("'threshold' must be NULL or a single positive finite number.",
         call. = FALSE)
  }
  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("'plot' must be either TRUE or FALSE.", call. = FALSE)
  }

  working_model <- .frailty_working_model(object)
  df_info <- .frailty_degrees_of_freedom(object, working_model$x_work,
                                         df_method = df_method)
  group_idx <- .popsize_group_index(object, by = by)
  baseline <- popsize(object, by = by, bias_correction = FALSE,
                      total = length(group_idx) > 1L)

  v_alpha <- object$vcov[seq_len(object$p_alpha), seq_len(object$p_alpha),
                         drop = FALSE]

  target_details <- lapply(seq_along(group_idx), function(i) {
    .frailty_target_stats(
      object = object,
      idx = group_idx[[i]],
      label = names(group_idx)[i],
      v_alpha = v_alpha,
      df = df_info$df
    )
  })
  if (length(group_idx) > 1L) {
    target_details[[length(target_details) + 1L]] <- .frailty_target_stats(
      object = object,
      idx = seq_len(object$n_obs),
      label = "Total",
      v_alpha = v_alpha,
      df = df_info$df
    )
  }

  working <- do.call(rbind, lapply(target_details, function(x) x$summary))
  surface <- .frailty_surface(working, r2_d = r2_d, r2_y = r2_y)
  robustness <- .frailty_robustness(
    working,
    q = q,
    threshold = threshold,
    direction = direction
  )

  out <- list(
    baseline = baseline,
    working = working,
    surface = surface,
    robustness = robustness,
    settings = list(
      by = by,
      r2_d = r2_d,
      r2_y = r2_y,
      q = q,
      threshold = threshold,
      direction = direction,
      df_method = df_method,
      resolved_df_method = df_info$resolved_method
    )
  )
  class(out) <- "uncounted_frailty_sensitivity"

  if (plot) {
    plot(out, ...)
  }

  out
}

#' @rdname frailty_sensitivity
#' @export
print.uncounted_frailty_sensitivity <- function(x, ...) {
  cat("Omitted-frailty sensitivity analysis\n")
  cat("Targets:", nrow(x$working), "| DF method:",
      x$settings$resolved_df_method, "\n\n")

  base_tab <- x$baseline[, c("group", "estimate"), drop = FALSE]
  names(base_tab)[2] <- "xi_hat"
  print.data.frame(base_tab, row.names = FALSE, ...)

  total_info <- attr(x$baseline, "total")
  if (!is.null(total_info)) {
    cat("\nTotal xi:", format(round(total_info$estimate), big.mark = ","), "\n")
  }

  cat("\nWorking contrasts:\n")
  print.data.frame(
    x$working[, c("group", "theta_hat", "se_theta", "d_xi_dtheta", "df"),
              drop = FALSE],
    row.names = FALSE, ...
  )

  cat("\nRobustness values:\n")
  print.data.frame(x$robustness, row.names = FALSE, ...)

  invisible(x)
}

#' @rdname frailty_sensitivity
#' @export
plot.uncounted_frailty_sensitivity <- function(x, group = NULL,
                                               type = c("ratio", "xi"),
                                               side = c("lower", "upper"),
                                               ...) {
  type <- match.arg(type)
  side <- match.arg(side)

  if (is.null(group)) {
    group <- unique(x$working$group)[1]
  }
  if (!group %in% unique(x$working$group)) {
    stop("'group' must match one of the groups in the frailty-sensitivity output.",
         call. = FALSE)
  }

  sdf <- x$surface[x$surface$group == group, , drop = FALSE]
  if (nrow(sdf) == 0L) {
    stop("No sensitivity surface available for the selected group.",
         call. = FALSE)
  }

  x_vals <- sort(unique(sdf$r2_d))
  y_vals <- sort(unique(sdf$r2_y))
  z_col <- if (identical(type, "ratio")) {
    if (identical(side, "lower")) "xi_ratio_lower" else "xi_ratio_upper"
  } else {
    if (identical(side, "lower")) "xi_lower" else "xi_upper"
  }
  z_mat <- matrix(sdf[[z_col]], nrow = length(x_vals), ncol = length(y_vals))
  z_lab <- if (identical(type, "ratio")) "Xi / Xi_hat" else "Xi"
  side_lab <- if (identical(side, "lower")) "Lower" else "Upper"

  graphics::contour(
    x = x_vals,
    y = y_vals,
    z = z_mat,
    xlab = expression(r[D]^2),
    ylab = expression(r[Y]^2),
    main = paste(side_lab, z_lab, "surface -", group),
    ...
  )

  invisible(x)
}

#' Validate frailty-sensitivity support for a fitted object
#' @noRd
.validate_frailty_object <- function(object) {
  if (!inherits(object, "uncounted")) {
    stop("'object' must inherit from class 'uncounted'.", call. = FALSE)
  }
  if (!identical(object$estimator, "mle")) {
    stop("frailty_sensitivity() is only available for models fitted with ",
         "estimator = 'mle'.", call. = FALSE)
  }
  if (!object$method %in% c("poisson", "nb")) {
    stop("frailty_sensitivity() is only available for Poisson and NB models.",
         call. = FALSE)
  }
  if (!identical(object$link_rho, "power")) {
    stop("frailty_sensitivity() currently supports only link_rho = 'power'.",
         call. = FALSE)
  }
  if (isTRUE(object$constrained)) {
    stop("frailty_sensitivity() is not available for constrained fits.",
         call. = FALSE)
  }
  if (isTRUE(object$has_cov_gamma)) {
    stop("frailty_sensitivity() is not available for models with cov_gamma.",
         call. = FALSE)
  }
  invisible(object)
}

#' Validate a partial-R2 grid
#' @noRd
.validate_r2_grid <- function(x, lower, upper, inclusive_upper = TRUE, arg) {
  if (!is.numeric(x) || length(x) == 0L || anyNA(x) || any(!is.finite(x))) {
    stop("'", arg, "' must be a non-empty numeric vector of finite values.",
         call. = FALSE)
  }
  if (any(x < lower) || any(if (inclusive_upper) x > upper else x >= upper)) {
    upper_txt <- if (inclusive_upper) paste0("<= ", upper) else paste0("< ", upper)
    stop("'", arg, "' values must satisfy ", lower, " <= ", arg, " ", upper_txt, ".",
         call. = FALSE)
  }
  sort(unique(as.numeric(x)))
}

#' Validate the q grid
#' @noRd
.validate_q_grid <- function(q, direction) {
  if (!is.numeric(q) || length(q) == 0L || anyNA(q) || any(!is.finite(q)) ||
      any(q < 0)) {
    stop("'q' must be a non-empty numeric vector of finite values >= 0.",
         call. = FALSE)
  }
  if (identical(direction, "decrease") && any(q > 1)) {
    stop("For direction = 'decrease', 'q' values must be <= 1.", call. = FALSE)
  }
  sort(unique(as.numeric(q)))
}

#' Build the IRLS working model used by frailty sensitivity
#' @noRd
.frailty_working_model <- function(object) {
  ratio <- object$n_aux / object$N
  rate <- if (!is.null(object$gamma)) object$gamma + ratio else ratio
  log_rate <- log(rate)

  y_star <- object$log_mu + (object$m - object$fitted.values) / object$fitted.values
  weights <- if (identical(object$method, "poisson")) {
    object$fitted.values
  } else {
    object$fitted.values / (1 + object$fitted.values / object$theta)
  }

  z_alpha <- object$X_alpha * log(object$N)
  z_beta <- object$X_beta * log_rate
  gamma_deriv <- if (isTRUE(object$gamma_estimated) && !is.null(object$gamma)) {
    matrix(
      object$beta_values / rate,
      ncol = 1,
      dimnames = list(NULL, "gamma")
    )
  } else {
    NULL
  }

  x_work <- cbind(z_alpha, z_beta)
  if (!is.null(gamma_deriv)) {
    x_work <- cbind(x_work, gamma_deriv)
  }

  list(
    y_star = y_star,
    weights = weights,
    z_alpha = z_alpha,
    z_beta = z_beta,
    gamma_deriv = gamma_deriv,
    x_work = x_work
  )
}

#' Resolve degrees of freedom for frailty sensitivity
#' @noRd
.frailty_degrees_of_freedom <- function(object, x_work,
                                        df_method = c("auto", "model", "cluster")) {
  df_method <- match.arg(df_method)

  rank_x <- qr(x_work)$rank
  df_model <- object$n_obs - rank_x
  if (!isTRUE(df_model > 0)) {
    stop("Model-based frailty-sensitivity degrees of freedom are not positive.",
         call. = FALSE)
  }

  df_cluster <- NA_real_
  if (!is.null(object$cluster_var)) {
    n_clusters <- length(unique(object$cluster_var[!is.na(object$cluster_var)]))
    df_cluster <- n_clusters - rank_x
  }

  if (identical(df_method, "model")) {
    return(list(df = df_model, resolved_method = "model"))
  }
  if (identical(df_method, "cluster")) {
    if (!is.finite(df_cluster) || !isTRUE(df_cluster > 0)) {
      warning("Cluster-based degrees of freedom are not positive; falling back ",
              "to model-based degrees of freedom.", call. = FALSE)
      return(list(df = df_model, resolved_method = "model"))
    }
    return(list(df = df_cluster, resolved_method = "cluster"))
  }

  if (is.finite(df_cluster) && isTRUE(df_cluster > 0)) {
    list(df = df_cluster, resolved_method = "cluster")
  } else {
    list(df = df_model, resolved_method = "model")
  }
}

#' Compute group-specific frailty target statistics
#' @noRd
.frailty_target_stats <- function(object, idx, label, v_alpha, df) {
  x_alpha <- object$X_alpha[idx, , drop = FALSE]
  n_g <- object$N[idx]
  alpha_g <- object$alpha_values[idx]
  xi_hat <- sum(n_g^alpha_g)
  grad_alpha <- colSums(
    x_alpha * as.numeric(n_g^alpha_g * log(n_g))
  )

  grad_norm <- sqrt(sum(grad_alpha^2))
  if (is.finite(grad_norm) && grad_norm > sqrt(.Machine$double.eps)) {
    contrast <- grad_alpha / grad_norm
    d_xi_dtheta <- sum(grad_alpha * contrast)
  } else {
    contrast <- rep(0, object$p_alpha)
    contrast[1] <- 1
    d_xi_dtheta <- 0
  }

  theta_hat <- sum(contrast * object$alpha_coefs)
  se_theta <- sqrt(max(as.numeric(t(contrast) %*% v_alpha %*% contrast), 0))

  list(
    summary = data.frame(
      group = label,
      xi_hat = xi_hat,
      theta_hat = theta_hat,
      se_theta = se_theta,
      d_xi_dtheta = d_xi_dtheta,
      df = df,
      stringsAsFactors = FALSE
    ),
    contrast = contrast,
    gradient = grad_alpha
  )
}

#' Build the frailty-sensitivity surface
#' @noRd
.frailty_surface <- function(working, r2_d, r2_y) {
  grid <- expand.grid(
    r2_d = r2_d,
    r2_y = r2_y,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  out <- lapply(seq_len(nrow(working)), function(i) {
    wrow <- working[i, , drop = FALSE]
    bias_theta <- wrow$se_theta * sqrt(wrow$df * grid$r2_y * grid$r2_d /
                                         pmax(1 - grid$r2_d, .Machine$double.eps))
    xi_shift <- wrow$d_xi_dtheta * bias_theta
    xi_lower <- pmax(wrow$xi_hat - xi_shift, 0)
    xi_upper <- wrow$xi_hat + xi_shift

    data.frame(
      group = wrow$group,
      r2_d = grid$r2_d,
      r2_y = grid$r2_y,
      bias_theta = bias_theta,
      xi_lower = xi_lower,
      xi_upper = xi_upper,
      xi_ratio_lower = if (isTRUE(wrow$xi_hat > 0)) xi_lower / wrow$xi_hat else NA_real_,
      xi_ratio_upper = if (isTRUE(wrow$xi_hat > 0)) xi_upper / wrow$xi_hat else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

#' Compute robustness values for frailty sensitivity
#' @noRd
.frailty_robustness <- function(working, q, threshold = NULL,
                                direction = c("decrease", "increase")) {
  direction <- match.arg(direction)

  out <- lapply(seq_len(nrow(working)), function(i) {
    wrow <- working[i, , drop = FALSE]
    xi_hat <- wrow$xi_hat
    se_theta <- wrow$se_theta
    d_xi_dtheta <- wrow$d_xi_dtheta
    df <- wrow$df

    targets <- if (is.null(threshold)) {
      if (identical(direction, "decrease")) xi_hat * (1 - q) else xi_hat * (1 + q)
    } else {
      threshold
    }
    q_vals <- if (is.null(threshold)) q else NA_real_

    do.call(rbind, lapply(seq_along(targets), function(j) {
      target <- targets[j]
      already_crossed <- if (identical(direction, "decrease")) {
        xi_hat <= target
      } else {
        xi_hat >= target
      }

      if (already_crossed) {
        rv_equal <- 0
        rv_extreme <- 0
        status <- "already_crossed"
      } else if (!is.finite(d_xi_dtheta) || d_xi_dtheta <= sqrt(.Machine$double.eps) ||
                 !is.finite(se_theta) || se_theta <= sqrt(.Machine$double.eps) ||
                 !is.finite(df) || df <= 0) {
        rv_equal <- 1
        rv_extreme <- 1
        status <- "flat_gradient"
      } else {
        delta_xi <- abs(target - xi_hat)
        delta_theta <- delta_xi / d_xi_dtheta
        f_val <- delta_theta / (se_theta * sqrt(df))

        if (!is.finite(f_val)) {
          rv_equal <- 1
          rv_extreme <- 1
          status <- "flat_gradient"
        } else {
          rv_equal <- 0.5 * (sqrt(f_val^4 + 4 * f_val^2) - f_val^2)
          rv_extreme <- f_val^2 / (1 + f_val^2)
          rv_equal <- min(max(rv_equal, 0), 1)
          rv_extreme <- min(max(rv_extreme, 0), 1)
          status <- "ok"
        }
      }

      data.frame(
        group = wrow$group,
        direction = direction,
        q = q_vals[j],
        threshold = if (is.null(threshold)) NA_real_ else threshold,
        xi_hat = xi_hat,
        target_xi = target,
        rv_equal = rv_equal,
        rv_extreme = rv_extreme,
        status = status,
        stringsAsFactors = FALSE
      )
    }))
  })

  do.call(rbind, out)
}
