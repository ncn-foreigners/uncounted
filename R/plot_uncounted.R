#' Diagnostic plots for uncounted models
#'
#' Produces a 4-panel diagnostic display following the recommendations of
#' Zhang (2008) for the power-law hidden-population model. The panels are:
#'
#' \enumerate{
#'   \item \strong{Fitted vs Observed (sqrt scale):}
#'     Scatterplot of \eqn{\sqrt{m}} against \eqn{\sqrt{\hat\mu}} with a
#'     45-degree equality line. Points should cluster around the line if the
#'     mean function is correctly specified.
#'   \item \strong{Anscombe residuals vs \eqn{\sqrt{\hat\mu}}:}
#'     Anscombe (variance-stabilised) residuals plotted against
#'     \eqn{\sqrt{\hat\mu}} with a running-mean smoother. A flat trend near
#'     zero indicates no systematic misfit.
#'   \item \strong{Scale-Location (Anscombe):}
#'     Absolute Anscombe residuals \eqn{|r_A|} against \eqn{\sqrt{\hat\mu}}
#'     with a running-mean smoother. An increasing trend suggests
#'     under-modelled variance (e.g., overdispersion).
#'   \item \strong{Normal Q-Q Plot:}
#'     Normal quantile-quantile plot of the Anscombe residuals. For
#'     well-specified count models the residuals should be approximately
#'     standard normal.
#' }
#'
#' @param x An object of class \code{"uncounted"}.
#' @param which Integer vector selecting which panels to display (default
#'   \code{1:4}).
#' @param ask Logical; if \code{TRUE}, the user is prompted before each new
#'   page. Defaults to \code{TRUE} in interactive sessions when more than one
#'   panel is requested, and \code{FALSE} otherwise.
#' @param ... Additional graphical arguments passed to \code{\link{plot}}.
#'
#' @references
#' Zhang, L.-C. (2008). Developing methods for determining the number of
#' unauthorized foreigners in Norway. \emph{Documents} 2008/11, Statistics
#' Norway. \url{https://www.ssb.no/a/english/publikasjoner/pdf/doc_200811_en/doc_200811_en.pdf}
#'
#' McCullagh, P. and Nelder, J. A. (1989). \emph{Generalized Linear Models},
#' 2nd ed. Chapman & Hall. \doi{10.1007/978-1-4899-3242-6}
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(
#'   N = rep(1000, 30),
#'   n = rpois(30, lambda = 50)
#' )
#' df$m <- rpois(30, lambda = df$N^0.5 * (df$n / df$N)^0.8)
#'
#' fit <- estimate_hidden_pop(data = df, observed = ~m, auxiliary = ~n,
#'                            reference_pop = ~N, method = "poisson")
#'
#' # All four diagnostic panels
#' plot(fit, ask = FALSE)
#'
#' # Only the Q-Q plot
#' plot(fit, which = 4)
#'
#' @export
plot.uncounted <- function(x, which = 1:4, ask = interactive() && length(which) > 1, ...) {
  m <- x$m
  mu <- x$fitted.values
  r_a <- residuals(x, type = "anscombe")
  sqrt_mu <- sqrt(mu)
  sqrt_m <- sqrt(m)

  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  running_mean <- .running_mean

  # Panel 1: sqrt(m) vs sqrt(mu)
  if (1 %in% which) {
    plot(sqrt_mu, sqrt_m, xlab = expression(sqrt(hat(mu))),
         ylab = expression(sqrt(m)),
         main = "Fitted vs Observed (sqrt scale)", pch = 1, ...)
    abline(0, 1, lty = 1)
  }

  # Panel 2: Anscombe residuals vs sqrt(mu)
  if (2 %in% which) {
    plot(sqrt_mu, r_a, xlab = expression(sqrt(hat(mu))),
         ylab = "Anscombe residual",
         main = "Anscombe residuals", pch = 1, ...)
    abline(h = 0, lty = 2)
    rm <- running_mean(sqrt_mu, r_a)
    lines(rm$x, rm$y, lty = 1, lwd = 2)
  }

  # Panel 3: |Anscombe residuals| vs sqrt(mu)
  if (3 %in% which) {
    plot(sqrt_mu, abs(r_a), xlab = expression(sqrt(hat(mu))),
         ylab = expression("|" * r[A] * "|"),
         main = "Scale-Location (Anscombe)", pch = 1, ...)
    rm <- running_mean(sqrt_mu, abs(r_a))
    lines(rm$x, rm$y, lty = 2, lwd = 2)
  }

  # Panel 4: Normal Q-Q
  if (4 %in% which) {
    qqnorm(r_a, main = "Normal Q-Q Plot", pch = 1, ...)
    qqline(r_a, lty = 1)
  }
}

#' Exploratory plots for the power-law model
#'
#' Produces two side-by-side log-log scatterplots that visualise the marginal
#' relationships underlying the power-law model
#' \eqn{\log(\mu/N) = (\alpha - 1)\log N + \beta\log(n/N)}, following
#' Zhang (2008, Figures 5, 8, 11):
#'
#' \enumerate{
#'   \item \strong{log(m/N) vs log(N):} Explores the scaling of the
#'     apprehension rate with total population size. The OLS slope approximates
#'     \eqn{\alpha - 1}. A slope near zero (\eqn{\alpha \approx 1}) means the
#'     rate is independent of population size.
#'   \item \strong{log(m/N) vs log(n/N):} Explores the relationship between
#'     the apprehension rate and the auxiliary-to-population ratio. The OLS
#'     slope approximates \eqn{\beta}.
#' }
#'
#' Observations with \eqn{m = 0} or \eqn{n = 0} are excluded (undefined on the
#' log scale). OLS regression lines are overlaid as visual guides.
#'
#' @param object An object of class \code{"uncounted"}.
#' @param ... Additional graphical arguments passed to \code{\link{plot}}.
#'
#' @references
#' Zhang, L.-C. (2008). Developing methods for determining the number of
#' unauthorized foreigners in Norway. \emph{Documents} 2008/11, Statistics
#' Norway. \url{https://www.ssb.no/a/english/publikasjoner/pdf/doc_200811_en/doc_200811_en.pdf}
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(
#'   N = round(exp(rnorm(50, 6, 1))),
#'   n = rpois(50, lambda = 30)
#' )
#' df$m <- rpois(50, lambda = df$N^0.5 * (df$n / df$N)^0.8)
#'
#' fit <- estimate_hidden_pop(data = df, observed = ~m, auxiliary = ~n,
#'                            reference_pop = ~N, method = "poisson")
#' plot_explore(fit)
#'
#' @export
plot_explore <- function(object, ...) {
  if (!inherits(object, "uncounted")) stop("object must be of class 'uncounted'")

  m <- object$m; N <- object$N; n_aux <- object$n_aux
  pos <- m > 0 & n_aux > 0
  log_mN <- log(m[pos] / N[pos])
  log_N <- log(N[pos])
  log_nN <- log(n_aux[pos] / N[pos])

  old_par <- par(mfrow = c(1, 2))
  on.exit(par(old_par))

  plot(log_N, log_mN, xlab = "log(N)", ylab = "log(m/N)",
       main = "log(m/N) vs log(N)", pch = 1, ...)
  fit1 <- lm(log_mN ~ log_N)
  abline(fit1, lty = 2)

  plot(log_nN, log_mN, xlab = "log(n/N)", ylab = "log(m/N)",
       main = "log(m/N) vs log(n/N)", pch = 1, ...)
  fit2 <- lm(log_mN ~ log_nN)
  abline(fit2, lty = 2)
}

#' Rootogram for count model diagnostics
#'
#' Compares observed and fitted count frequencies using a rootogram
#' (Kleiber & Zeileis, 2016). Only available for Poisson and Negative Binomial
#' models. Three display styles are supported:
#'
#' \describe{
#'   \item{\strong{hanging} (default)}{Bars representing the observed
#'     frequencies are "hung" from the fitted curve. The bottom of each bar
#'     should touch the zero line when the fit is good. Deviations below zero
#'     indicate underfitting; above zero, overfitting at that count.}
#'   \item{\strong{standing}}{Observed bars stand on the x-axis and the fitted
#'     curve is overlaid. Harder to judge discrepancies than the hanging style
#'     because the eye must compare bar height to the curve.}
#'   \item{\strong{suspended}}{Plots the difference
#'     \eqn{\sqrt{f_{\mathrm{obs}}} - \sqrt{f_{\mathrm{exp}}}} directly. Bars
#'     that cross zero indicate counts where the model over- or under-predicts.}
#' }
#'
#' By default the y-axis is on the square-root scale (\code{sqrt_scale = TRUE}),
#' which stabilises the visual variability of the bars and makes it easier to
#' judge fit in the tails.
#'
#' @param object An object of class \code{"uncounted"} fitted with
#'   \code{method = "poisson"} or \code{method = "nb"}.
#' @param style Character: \code{"hanging"} (default), \code{"standing"}, or
#'   \code{"suspended"}.
#' @param max_count Maximum count value to display. If \code{NULL}
#'   (default), automatically determined from the data (capped at 100).
#' @param sqrt_scale Logical; use the square-root scale for frequencies
#'   (default \code{TRUE}).
#' @param ... Additional arguments passed to \code{\link{barplot}}.
#'
#' @return Invisibly, a list with components \code{observed} (observed
#'   frequencies), \code{expected} (expected frequencies from the fitted model),
#'   \code{counts} (count values), \code{style}, and \code{sqrt_scale}.
#'
#' @references
#' Kleiber, C. and Zeileis, A. (2016). Visualizing Count Data Regressions
#' Using Rootograms. \emph{The American Statistician}, 70(3), 296--303.
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(
#'   N = rep(1000, 50),
#'   n = rpois(50, lambda = 50)
#' )
#' df$m <- rpois(50, lambda = df$N^0.5 * (df$n / df$N)^0.8)
#'
#' fit <- estimate_hidden_pop(data = df, observed = ~m, auxiliary = ~n,
#'                            reference_pop = ~N, method = "poisson")
#'
#' # Hanging rootogram (default)
#' rootogram(fit)
#'
#' # Standing rootogram
#' rootogram(fit, style = "standing")
#'
#' # Suspended rootogram
#' rootogram(fit, style = "suspended")
#'
#' @export
rootogram <- function(object, ...) UseMethod("rootogram")

#' @rdname rootogram
#' @export
rootogram.uncounted <- function(object, style = c("hanging", "standing", "suspended"),
                                 max_count = NULL, sqrt_scale = TRUE, ...) {
  style <- match.arg(style)

  if (!object$method %in% c("poisson", "nb")) {
    stop("Rootogram is only available for Poisson and NB models.")
  }
  if (isTRUE(object$pseudo_loglik)) {
    stop("Rootogram needs whole-number counts; 'observed' has non-integer values.",
         call. = FALSE)
  }

  m <- object$m
  mu <- object$fitted.values
  theta <- object$theta

  if (is.null(max_count)) {
    max_count <- min(max(m), quantile(m, 0.99) + 10)
    if (max_count > 100) max_count <- 100
  }
  max_count <- as.integer(max_count)

  # Observed frequencies
  breaks <- 0:max_count
  obs <- tabulate(factor(pmin(m, max_count), levels = breaks), nbins = length(breaks))

  # Expected frequencies
  exp_freq <- numeric(length(breaks))
  for (k in seq_along(breaks)) {
    if (object$method == "poisson") {
      exp_freq[k] <- sum(dpois(breaks[k], mu))
    } else {
      exp_freq[k] <- sum(dnbinom(breaks[k], size = theta, mu = mu))
    }
  }
  # Last bin is cumulative (>= max_count)
  if (max_count < max(m)) {
    if (object$method == "poisson") {
      exp_freq[length(breaks)] <- sum(ppois(max_count - 1, mu, lower.tail = FALSE))
    } else {
      exp_freq[length(breaks)] <- sum(pnbinom(max_count - 1, size = theta, mu = mu, lower.tail = FALSE))
    }
  }

  # Transform
  obs_t <- if (sqrt_scale) sqrt(obs) else obs
  exp_t <- if (sqrt_scale) sqrt(exp_freq) else exp_freq

  # Plot
  ylab <- if (sqrt_scale) expression(sqrt(Frequency)) else "Frequency"

  x_pos <- barplot(
    switch(style,
      "hanging" = exp_t,
      "standing" = obs_t,
      "suspended" = obs_t - exp_t
    ),
    names.arg = breaks, xlab = "Count", ylab = ylab,
    main = paste0("Rootogram (", style, ") - ", toupper(object$method)),
    col = if (style == "hanging") "white" else "lightblue",
    border = "gray40", ...
  )

  if (style == "hanging") {
    # Draw observed bars hanging from top of expected
    rect(x_pos - 0.5, exp_t - obs_t, x_pos + 0.5, exp_t,
         col = "lightblue", border = "gray40")
    lines(x_pos, exp_t, type = "b", pch = 19, col = "red", cex = 0.5)
    abline(h = 0, lty = 2)
  } else if (style == "standing") {
    lines(x_pos, exp_t, type = "b", pch = 19, col = "red", cex = 0.5)
  } else {
    abline(h = 0, lty = 2, col = "red")
  }

  invisible(list(observed = obs, expected = exp_freq, counts = breaks,
                 style = style, sqrt_scale = sqrt_scale))
}


#' Gamma Profile: Population Size as a Function of Gamma
#'
#' Refits the model across a grid of fixed gamma values and plots the
#' estimated population size \eqn{\hat\xi(\gamma)} and log-likelihood
#' \eqn{\ell(\gamma)} as functions of gamma. Useful for assessing
#' identification strength of gamma and sensitivity of population
#' size estimates.
#'
#' @details
#' For each gamma value in the grid, the model is refitted with
#' \code{gamma = gamma_val} (fixed) and the total estimated population
#' size \eqn{\hat\xi = \sum_i N_i^{\hat\alpha}} is extracted. If the
#' original model estimated gamma, its estimate is marked on the plot.
#'
#' A flat \eqn{\xi(\gamma)} profile indicates that population size is
#' robust to gamma specification. A steep profile suggests weak
#' identification and sensitivity to the gamma assumption.
#'
#' @param object An \code{"uncounted"} object (Poisson or NB).
#' @param gamma_grid Numeric vector of gamma values to evaluate.
#'   Default is 20 points from 1e-4 to 0.5.
#' @param plot Logical; produce the plot? Default TRUE.
#' @param ... Additional arguments passed to \code{plot()}.
#'
#' @return Invisibly, a data frame with columns: \code{gamma}, \code{xi},
#'   \code{loglik}.
#'
#' @examples
#' set.seed(42)
#' df <- data.frame(
#'   N = round(exp(rnorm(40, 5, 1.5))),
#'   n = rpois(40, 20)
#' )
#' df$m <- rpois(40, df$N^0.7 * (0.01 + df$n / df$N)^0.5)
#' fit <- estimate_hidden_pop(df, ~m, ~n, ~N, method = "poisson")
#' profile_gamma(fit)
#'
#' @export
profile_gamma <- function(object, gamma_grid = seq(1e-4, 0.5, length.out = 20),
                          plot = TRUE, ...) {
  if (!identical(object$estimator, "mle")) {
    stop("profile_gamma() is only available for models fitted with ",
         "estimator = 'mle'.", call. = FALSE)
  }
  if (!object$method %in% c("poisson", "nb")) {
    stop("profile_gamma() is only available for Poisson and NB models.")
  }
  if (isTRUE(object$has_cov_gamma)) {
    stop("profile_gamma() is not available for models with covariate-varying ",
         "gamma (cov_gamma). Profiling a scalar gamma grid does not apply when ",
         "gamma varies by covariates.", call. = FALSE)
  }

  ## Refit with exactly the original specification, changing only gamma.
  ## Building this by hand previously dropped weights, cluster and countries,
  ## so a weighted or clustered fit was profiled by unweighted, unclustered
  ## refits and the profile did not belong to the model being profiled.
  refit_args <- .refit_args_from_object(object)
  ## cov_gamma requires gamma = "estimate" and is rejected above anyway, so it
  ## must not be carried into a fixed-gamma refit.
  refit_args$cov_gamma <- NULL

  results <- data.frame(gamma = gamma_grid, xi = NA_real_, loglik = NA_real_)

  for (i in seq_along(gamma_grid)) {
    args_i <- refit_args
    args_i$gamma <- gamma_grid[i]
    fit_i <- tryCatch(
      do.call(estimate_hidden_pop, args_i),
      error = function(e) NULL
    )
    if (!is.null(fit_i)) {
      ps <- popsize(fit_i)
      results$xi[i] <- sum(ps$estimate)
      results$loglik[i] <- fit_i$loglik
    }
  }

  if (plot) {
    op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    on.exit(par(op))

    # xi vs gamma
    plot(results$gamma, results$xi, type = "b", pch = 19, cex = 0.6,
         xlab = expression(gamma), ylab = expression(hat(xi)),
         main = expression("Population size" ~ hat(xi)(gamma)),
         ...)
    if (object$gamma_estimated && !is.null(object$gamma)) {
      abline(v = object$gamma, lty = 2, col = "red")
      mtext(sprintf("est. gamma = %.4f", object$gamma), side = 3,
            line = 0, cex = 0.7, col = "red")
    }

    # loglik vs gamma
    plot(results$gamma, results$loglik, type = "b", pch = 19, cex = 0.6,
         xlab = expression(gamma), ylab = "Log-likelihood",
         main = expression("Log-likelihood" ~ ell(gamma)),
         ...)
    if (object$gamma_estimated && !is.null(object$gamma)) {
      abline(v = object$gamma, lty = 2, col = "red")
    }
  }

  invisible(results)
}


# ---- Smoother for residual plots ----

#' Running mean over consecutive blocks covering all points, in x order
#' @noRd
.running_mean <- function(x_val, y_val, n_windows = 10) {
  ord <- order(x_val)
  n <- length(ord)
  ## Blocks of nearly equal size: integer division would leave out up to
  ## n_windows - 1 of the largest values
  block <- ceiling(seq_len(n) * min(n_windows, n) / n)
  list(x = as.numeric(tapply(x_val[ord], block, mean)),
       y = as.numeric(tapply(y_val[ord], block, mean)))
}


# ---- NLL helper for profiling ----

#' Build a negative log-likelihood closure from a fitted uncounted object.
#'
#' The criterion is the one the method optimises: the count likelihoods for
#' Poisson/NB, and a weighted Gaussian on the log scale for OLS and on the raw
#' scale for NLS. \code{nll} takes (alpha, beta) with gamma and theta held at the
#' fit; \code{nll_full} also takes the nuisance parameters, appended as
#' (log theta, gamma coefficients on the log scale) in \code{par_full}, with
#' their bounds in \code{lower} and \code{upper}.
#' @param object An uncounted object.
#' @return A list with nll, par0, compute_xi, nll_full, par_full, lower, upper.
#' @noRd
.build_nll_from_object <- function(object) {
  if (!identical(object$estimator, "mle")) {
    stop("Profile likelihood is only available for models fitted with ",
         "estimator = 'mle'.", call. = FALSE)
  }
  method <- object$method
  if (method == "iols") {
    stop("Profile likelihood is not available for method = 'iols', which ",
         "solves estimating equations rather than maximising a likelihood.",
         call. = FALSE)
  }
  m <- object$m
  N <- object$N
  n_aux <- object$n_aux
  ratio <- n_aux / N
  log_N <- log(N)
  n_obs <- object$n_obs
  X_alpha <- object$X_alpha
  X_beta <- object$X_beta
  p_alpha <- object$p_alpha
  p_beta <- object$p_beta
  p_ab <- p_alpha + p_beta
  weights <- if (!is.null(object$obs_weights)) object$obs_weights else rep(1, n_obs)
  is_constr <- isTRUE(object$constrained)
  link_rho <- object$link_rho
  theta_mle <- object$theta

  .alpha_fn <- if (is_constr) .inv_logit else identity
  .beta_fn  <- if (is_constr) function(x) exp(x) else identity

  # Use fitted gamma (scalar or vector)
  if (!is.null(object$gamma_values) && isTRUE(object$has_cov_gamma)) {
    rate_fit <- object$gamma_values + ratio
  } else if (!is.null(object$gamma)) {
    rate_fit <- object$gamma + ratio
  } else {
    rate_fit <- ratio
  }

  par0 <- c(object$alpha_coefs, object$beta_coefs)

  ## Profile the same criterion the fit optimised, including any dependence
  ## offset that was in force then.
  offset <- object$dependence_offset
  if (is.null(offset)) offset <- 0

  nll_at <- function(par, rate, theta) {
    a <- par[seq_len(p_alpha)]
    b <- par[p_alpha + seq_len(p_beta)]
    alpha_lin <- .alpha_fn(as.numeric(X_alpha %*% a))
    beta_lin <- .beta_fn(as.numeric(X_beta %*% b))
    log_mu <- pmin(.compute_log_mu(alpha_lin, log_N, beta_lin, rate,
                                   link_rho = link_rho,
                                   dependence_offset = offset), 20)
    mu <- pmax(exp(log_mu), 1e-300)

    val <- switch(method,
      nb = -sum(weights * dnbinom(m, size = theta, mu = mu, log = TRUE)),
      poisson = -sum(weights * (m * log(mu) - mu - lgamma(m + 1))),
      {
        resid <- if (method == "nls") m - mu else .log_m_transform(m) - log_mu
        sum_w <- sum(weights)
        sum_w / 2 * (log(2 * pi * sum(weights * resid^2) / sum_w) + 1)
      }
    )
    if (!is.finite(val)) return(1e20)
    val
  }

  nll_fn <- function(par) nll_at(par, rate_fit, theta_mle)

  # Nuisance parameters, on the scales the fitters optimise
  has_theta <- method == "nb"
  est_gamma <- isTRUE(object$gamma_estimated)
  gamma_par <- if (!est_gamma) {
    NULL
  } else if (!is.null(object$gamma_coefs)) {
    object$gamma_coefs
  } else {
    c(gamma = log(object$gamma))
  }
  p_gamma <- length(gamma_par)
  X_gamma <- if (p_gamma > 0L && !is.null(object$X_gamma)) {
    object$X_gamma
  } else if (p_gamma > 0L) {
    matrix(1, nrow = n_obs, ncol = 1L)
  }
  gamma_bounds <- object$call$gamma_bounds
  if (is.null(gamma_bounds)) gamma_bounds <- c(1e-10, 0.5)
  idx_theta <- p_ab + 1L
  idx_gamma <- p_ab + has_theta + seq_len(p_gamma)

  par_full <- c(par0, if (has_theta) c(theta = log(theta_mle)), gamma_par)
  lower <- rep(-Inf, length(par_full))
  upper <- rep(Inf, length(par_full))
  if (has_theta) {
    lower[idx_theta] <- log(0.01)
    upper[idx_theta] <- log(1e6)
  }
  if (p_gamma > 0L) {
    lower[idx_gamma[1]] <- log(gamma_bounds[1])
    upper[idx_gamma[1]] <- log(gamma_bounds[2])
  }

  nll_full <- function(par) {
    rate <- if (p_gamma > 0L) {
      exp(pmin(as.numeric(X_gamma %*% par[idx_gamma]), 10)) + ratio
    } else {
      rate_fit
    }
    nll_at(par, rate, if (has_theta) exp(par[idx_theta]) else theta_mle)
  }

  xi_fn <- function(par) {
    a <- par[seq_len(p_alpha)]
    alpha_lin <- .alpha_fn(as.numeric(X_alpha %*% a))
    sum(N^alpha_lin)
  }

  list(nll = nll_fn, par0 = par0, compute_xi = xi_fn, nll_full = nll_full,
       par_full = par_full, lower = lower, upper = upper)
}


# ---- Profile alpha / beta ----

#' Shared implementation for profiling alpha or beta coefficients.
#' @noRd
.profile_coef <- function(object, coef_index, block = c("alpha", "beta"),
                          grid = NULL, reoptimize = FALSE, plot = TRUE, ...) {
  if (!identical(object$estimator, "mle")) {
    stop("profile_alpha() and profile_beta() are only available for models ",
         "fitted with estimator = 'mle'.", call. = FALSE)
  }
  block <- match.arg(block)
  p_alpha <- object$p_alpha
  p_beta <- object$p_beta

  if (block == "alpha") {
    if (coef_index < 1 || coef_index > p_alpha)
      stop("coef_index must be between 1 and ", p_alpha, " for alpha.", call. = FALSE)
    full_idx <- coef_index
    coef_name <- names(object$alpha_coefs)[coef_index]
  } else {
    if (coef_index < 1 || coef_index > p_beta)
      stop("coef_index must be between 1 and ", p_beta, " for beta.", call. = FALSE)
    full_idx <- p_alpha + coef_index
    coef_name <- names(object$beta_coefs)[coef_index]
  }

  nll_info <- .build_nll_from_object(object)
  par0 <- nll_info$par0
  nll <- nll_info$nll
  xi_fn <- nll_info$compute_xi
  mle_val <- par0[full_idx]

  se_val <- tryCatch(
    sqrt(max(0, object$vcov[full_idx, full_idx])),
    error = function(e) abs(mle_val) * 0.1 + 0.01
  )
  if (!is.finite(se_val) || se_val <= 0) se_val <- abs(mle_val) * 0.1 + 0.01

  if (is.null(grid)) {
    grid <- seq(mle_val - 3 * se_val, mle_val + 3 * se_val, length.out = 30)
  }

  results <- data.frame(value = grid, xi = NA_real_, loglik = NA_real_)

  if (!reoptimize) {
    for (i in seq_along(grid)) {
      par_i <- par0
      par_i[full_idx] <- grid[i]
      results$loglik[i] <- -nll(par_i)
      results$xi[i] <- xi_fn(par_i)
    }
  } else {
    ## True profile: every other parameter is free, including gamma and theta
    par_full <- nll_info$par_full
    nll_full <- nll_info$nll_full
    free_idx <- setdiff(seq_along(par_full), full_idx)
    for (i in seq_along(grid)) {
      fixed_val <- grid[i]
      obj_reduced <- function(par_free) {
        par_i <- par_full
        par_i[free_idx] <- par_free
        par_i[full_idx] <- fixed_val
        nll_full(par_i)
      }
      fit_free <- function(start) {
        optim(start, obj_reduced, method = "L-BFGS-B",
              lower = nll_info$lower[free_idx], upper = nll_info$upper[free_idx],
              control = list(maxit = 500, factr = 1e5, ndeps = rep(1e-5, length(start))))
      }
      opt_i <- tryCatch({
        opt1 <- fit_free(par_full[free_idx])
        fit_free(opt1$par)  # restart: fresh curvature memory
      }, error = function(e) NULL)
      if (!is.null(opt_i)) {
        par_i <- par_full
        par_i[free_idx] <- opt_i$par
        par_i[full_idx] <- fixed_val
        results$loglik[i] <- -nll_full(par_i)
        results$xi[i] <- xi_fn(par_i)
      }
    }
  }

  if (plot) {
    op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    on.exit(par(op))

    plot(results$value, results$xi, type = "b", pch = 19, cex = 0.6,
         xlab = coef_name, ylab = expression(hat(xi)),
         main = bquote("Population size" ~ hat(xi)(.(coef_name))),
         ...)
    abline(v = mle_val, lty = 2, col = "red")

    valid_ll <- is.finite(results$loglik)
    if (any(valid_ll)) {
      plot(results$value[valid_ll], results$loglik[valid_ll],
           type = "b", pch = 19, cex = 0.6,
           xlab = coef_name, ylab = "Log-likelihood",
           main = bquote("Log-likelihood"(.(coef_name))))
      abline(v = mle_val, lty = 2, col = "red")
    }
  }

  invisible(results)
}


#' Profile Likelihood for Alpha Coefficients
#'
#' Evaluates the log-likelihood and population size over a grid of values
#' for a chosen alpha coefficient, holding other parameters at their MLE
#' (concentrated profile) or re-optimizing them (true profile).
#'
#' @param object An \code{"uncounted"} object.
#' @param coef_index Integer: which alpha coefficient to profile (1 = intercept).
#' @param grid Numeric vector of values to evaluate. If \code{NULL} (default),
#'   auto-generates 30 points spanning \eqn{\pm 3} SE around the MLE.
#' @param reoptimize Logical. If \code{TRUE}, optimize all other parameters at
#'   each grid point, including an estimated gamma and the NB theta (true
#'   profile likelihood). If \code{FALSE} (default), hold everything else at the
#'   MLE (concentrated profile, faster). The criterion is the one the method
#'   optimises: the count likelihood for Poisson/NB and a weighted Gaussian on
#'   the log scale for OLS and on the raw scale for NLS. Not available for
#'   \code{method = "iols"}.
#' @param plot Logical; produce the plot? Default \code{TRUE}.
#' @param ... Additional arguments passed to \code{plot()}.
#'
#' @return Invisibly, a data frame with columns: \code{value}, \code{xi},
#'   \code{loglik}.
#'
#' @seealso \code{\link{profile_beta}}, \code{\link{profile_gamma}},
#'   \code{\link{profile.uncounted}}
#'
#' @export
profile_alpha <- function(object, coef_index = 1, grid = NULL,
                          reoptimize = FALSE, plot = TRUE, ...) {
  .profile_coef(object, coef_index = coef_index, block = "alpha",
                grid = grid, reoptimize = reoptimize, plot = plot, ...)
}


#' Profile Likelihood for Beta Coefficients
#'
#' Evaluates the log-likelihood and population size over a grid of values
#' for a chosen beta coefficient. When \code{reoptimize = FALSE}, the
#' population size \eqn{\hat\xi} is constant because beta does not enter
#' \eqn{\xi = \sum N^{\alpha}} directly.
#'
#' @inheritParams profile_alpha
#'
#' @return Invisibly, a data frame with columns: \code{value}, \code{xi},
#'   \code{loglik}.
#'
#' @seealso \code{\link{profile_alpha}}, \code{\link{profile_gamma}},
#'   \code{\link{profile.uncounted}}
#'
#' @export
profile_beta <- function(object, coef_index = 1, grid = NULL,
                         reoptimize = FALSE, plot = TRUE, ...) {
  .profile_coef(object, coef_index = coef_index, block = "beta",
                grid = grid, reoptimize = reoptimize, plot = plot, ...)
}


#' Profile Likelihood for Uncounted Models
#'
#' S3 method for \code{\link[stats]{profile}} that dispatches to
#' \code{\link{profile_gamma}}, \code{\link{profile_alpha}},
#' \code{\link{profile_beta}}, or \code{\link{profile_dependence}}
#' depending on \code{param}.
#'
#' @param fitted An \code{"uncounted"} object.
#' @param param Character: which parameter to profile. One of
#'   \code{"gamma"}, \code{"alpha"}, \code{"beta"}, or
#'   \code{"dependence"}.
#' @param ... Additional arguments passed to the specific profile function.
#'
#' @return Invisibly, a data frame from the dispatched function.
#'
#' @importFrom stats profile
#' @export
profile.uncounted <- function(fitted,
                              param = c("gamma", "alpha", "beta", "dependence"),
                              ...) {
  param <- match.arg(param)
  switch(param,
    gamma = profile_gamma(fitted, ...),
    alpha = profile_alpha(fitted, ...),
    beta = profile_beta(fitted, ...),
    dependence = profile_dependence(fitted, ...)
  )
}
