#' Estimate the Size of an Unauthorized Migrant Population
#'
#' Fits a power-law model relating observed counts of unauthorized migrants
#' to reference population sizes and auxiliary registration counts. The model
#' is estimated using OLS, NLS, Poisson MLE, or Negative Binomial MLE, and
#' supports covariate-varying parameters through interaction formulas.
#'
#' @param data A data frame containing all variables.
#' @param observed One-sided formula for the observed count (e.g., \code{~ m}).
#' @param auxiliary One-sided formula for the auxiliary count (e.g., \code{~ n}).
#' @param reference_pop One-sided formula for the reference population (e.g., \code{~ N}).
#' @param method Estimation method: \code{"poisson"} (default), \code{"nb"},
#'   \code{"ols"}, or \code{"nls"}. See Details.
#' @param cov_alpha Formula for covariates in alpha (e.g., \code{~ sex + year}).
#'   Default \code{NULL} means a single alpha (intercept only).
#' @param cov_beta Formula for covariates in beta (e.g., \code{~ sex}).
#'   Default \code{NULL} means a single beta (intercept only).
#' @param gamma Controls the gamma offset in the rate term:
#'   \itemize{
#'     \item \code{"estimate"} (default): gamma is estimated from data.
#'     \item A numeric value: fixed gamma, uses \eqn{\log(\gamma + n_i / N_i)}.
#'     \item \code{NULL}: no gamma, uses \eqn{\log(n_i / N_i)} (requires \eqn{n_i > 0}).
#'   }
#' @param gamma_bounds Numeric vector of length 2: lower and upper bounds for
#'   the estimated gamma parameter. Default \code{c(1e-10, 0.5)}.
#' @param theta_start Starting value for the NB dispersion parameter (used
#'   only when \code{method = "nb"}). Default 1.
#' @param vcov Controls the variance-covariance estimator. Can be:
#'   \itemize{
#'     \item A character string specifying the HC type: \code{"HC0"} through
#'       \code{"HC5"}, or \code{"HC4m"}. Default \code{"HC3"}. When \code{cluster}
#'       is also provided, this type is passed to \code{sandwich::vcovCL()}.
#'     \item A function that takes the fitted \code{uncounted} object and returns
#'       a variance-covariance matrix. For example,
#'       \code{sandwich::vcovHC} or
#'       \code{function(x) sandwich::vcovCL(x, cluster = data$country_code)}.
#'   }
#' @param cluster Optional one-sided formula identifying a cluster variable
#'   for cluster-robust variance estimation (e.g., \code{~ country_code}).
#'   When provided and \code{vcov} is a character string, the variance is
#'   computed using \code{sandwich::vcovCL()} with the specified HC type.
#'   Ignored when \code{vcov} is a function.
#' @param weights Optional numeric vector of observation weights.
#' @param constrained Logical. If \code{TRUE}, applies link functions to
#'   ensure \eqn{\alpha \in (0, 1)} (logit) and \eqn{\beta > 0} (exp).
#'   Only available for \code{method = "poisson"} and \code{method = "nb"}.
#'   Default \code{FALSE}. See Details.
#' @param countries One-sided formula identifying a country or group variable
#'   used when reporting population size estimates (e.g., \code{~ country}).
#'   This does \strong{not} enable cluster-robust variance; use \code{cluster}
#'   for that.
#'
#' @details
#' \strong{Model specification.}
#' For observation \eqn{i}, the expected observed count \eqn{m_i} is modelled as:
#'
#' \deqn{E(m_i) = N_i^{\alpha_i} \cdot (\gamma + n_i / N_i)^{\beta_i}}{E(m_i) = N_i^alpha_i * (gamma + n_i / N_i)^beta_i}
#'
#' where \eqn{N_i} is the reference (total registered) population,
#' \eqn{n_i} is an auxiliary count (e.g., new registrations),
#' and \eqn{\gamma \ge 0}{gamma >= 0} is an intercept-like offset.
#' On the log scale, the model is linear:
#'
#' \deqn{\log E(m_i) = \alpha_i \log N_i + \beta_i \log(\gamma + n_i / N_i)}{log E(m_i) = alpha_i * log(N_i) + beta_i * log(gamma + n_i / N_i)}
#'
#' \strong{Parameter interpretation.}
#' \describe{
#'   \item{\eqn{\alpha}}{Elasticity of the observed count with respect to
#'     the reference population. When \eqn{\alpha < 1}, the total
#'     (unobserved) population \eqn{\xi = \sum_i N_i^\alpha}{xi = sum(N_i^alpha)} is smaller
#'     than \eqn{\sum_i N_i}{sum(N_i)}, reflecting incomplete coverage. Values
#'     of \eqn{\alpha} near 0 imply weak dependence on population size;
#'     values near 1 imply near-proportional scaling.}
#'   \item{\eqn{\beta}}{Elasticity with respect to the registration rate
#'     \eqn{\gamma + n_i / N_i}. A positive \eqn{\beta} means higher
#'     auxiliary rates are associated with more observed unauthorized
#'     migrants.}
#'   \item{\eqn{\gamma}}{Baseline registration-rate offset, ensuring that
#'     the rate term is positive even when \eqn{n_i = 0}. Typically small
#'     (close to zero). When \code{gamma = "estimate"}, it is profiled out
#'     (OLS) or jointly optimized (MLE methods).}
#' }
#'
#' \strong{Estimation methods.}
#' \describe{
#'   \item{OLS}{Ordinary least squares on the log-linearized model.
#'     Fast and transparent but ignores the count nature of \eqn{m_i}
#'     and may be inefficient under heteroscedasticity.
#'     Uses \eqn{\log(m_i)} as the response (or \eqn{\log(m_i + 1)} when
#'     zeros are present). Note: when the \eqn{\log(m_i + 1)} transformation
#'     is used, fitted values are \eqn{\exp(\hat{\mu})} where \eqn{\hat{\mu}}
#'     was estimated on the shifted scale, so response-scale residuals and
#'     diagnostics are approximate.}
#'   \item{NLS}{Nonlinear least squares on the original scale
#'     \eqn{m_i = N_i^{\alpha} (\gamma + n_i/N_i)^{\beta} + \varepsilon_i}{m_i = N_i^alpha * (gamma + n_i/N_i)^beta + epsilon_i}.
#'     Avoids the log-transformation bias of OLS but still treats \eqn{m_i}
#'     as continuous.}
#'   \item{Poisson}{Poisson pseudo-maximum likelihood (PPML). Consistent
#'     under heteroscedasticity as long as the conditional mean is correctly
#'     specified (Santos Silva and Tenreyro, 2006). Recommended as the
#'     default.}
#'   \item{NB}{Negative Binomial MLE. Adds a dispersion parameter
#'     \eqn{\theta} to accommodate overdispersion beyond what the Poisson
#'     allows. Standard errors for the regression coefficients are computed
#'     conditional on \eqn{\hat{\theta}}, consistent with the approach used
#'     in \code{MASS::glm.nb}. This means coefficient SEs do not account for
#'     uncertainty in \eqn{\theta} estimation.}
#' }
#'
#' \strong{Constrained vs. unconstrained estimation.}
#' When \code{constrained = FALSE} (default), \eqn{\alpha} and \eqn{\beta}
#' are estimated on the real line without restrictions. When
#' \code{constrained = TRUE}, link functions enforce parameter bounds:
#' \eqn{\alpha = \mathrm{logit}(\eta)}{alpha = logit(eta)} so that
#' \eqn{\alpha \in (0, 1)}{alpha in (0, 1)}, and
#' \eqn{\beta = \exp(\zeta)}{beta = exp(zeta)} so that \eqn{\beta > 0}.
#' Coefficients are then reported on the link (transformed) scale; use
#' \code{summary()} or the \code{alpha_values} / \code{beta_values} elements
#' of the returned object for response-scale values.
#' Constrained estimation is available for the Poisson and NB methods only.
#'
#' \strong{Population size estimation.}
#' The total unauthorized population is estimated as
#' \eqn{\hat{\xi} = \sum_i N_i^{\hat{\alpha}_i}}{xi_hat = sum(N_i^alpha_hat_i)}.
#' See \code{\link{popsize}} for bias correction and confidence intervals.
#'
#' @return An object of class \code{"uncounted"}, a list containing:
#'   \describe{
#'     \item{\code{coefficients}}{Named vector of all estimated coefficients.}
#'     \item{\code{alpha_coefs}, \code{beta_coefs}}{Coefficient sub-vectors
#'       for the alpha and beta equations.}
#'     \item{\code{alpha_values}, \code{beta_values}}{Per-observation alpha
#'       and beta on the response scale (after applying link inverse when
#'       constrained).}
#'     \item{\code{vcov}}{HC-robust variance-covariance matrix.}
#'     \item{\code{vcov_model}}{Model-based (homoscedastic) variance-covariance
#'       matrix, used for bias correction.}
#'     \item{\code{fitted.values}}{Fitted values \eqn{\hat{m}_i}{m_hat_i}.}
#'     \item{\code{residuals}}{Raw residuals \eqn{m_i - \hat{m}_i}{m_i - m_hat_i}.}
#'     \item{\code{gamma}}{Estimated or fixed gamma value (or \code{NULL}).}
#'     \item{\code{theta}}{NB dispersion parameter (only for \code{method = "nb"}).}
#'     \item{\code{loglik}}{Log-likelihood (for Poisson and NB methods).}
#'     \item{\code{method}}{The estimation method used.}
#'   }
#'
#' @references
#' Zhang, L.-C. (2008). Developing methods for determining the number of
#' unauthorized foreigners in Norway. \emph{Documents} 2008/11, Statistics
#' Norway. \url{https://www.ssb.no/a/english/publikasjoner/pdf/doc_200811_en/doc_200811_en.pdf}
#'
#' Beręsewicz, M. and Pawlukiewicz, K. (2020). Estimation of the number of
#' irregular foreigners in Poland using non-linear count regression models.
#' \emph{arXiv preprint} arXiv:2008.09407.
#'
#' Santos Silva, J. M. C. and Tenreyro, S. (2006). The log of gravity.
#' \emph{The Review of Economics and Statistics}, 88(4), 641--658.
#'
#' @examples
#' # Simulate data: 50 groups with known population structure
#' set.seed(42)
#' sim_data <- data.frame(
#'   N = sample(1000:50000, 50, replace = TRUE)
#' )
#' sim_data$n <- rpois(50, lambda = sim_data$N * 0.05)
#' alpha_true <- 0.6
#' beta_true <- 1.2
#' gamma_true <- 0.01
#' mu <- sim_data$N^alpha_true * (gamma_true + sim_data$n / sim_data$N)^beta_true
#' sim_data$m <- rpois(50, lambda = mu)
#'
#' # Poisson with estimated gamma (default)
#' fit_pois <- estimate_hidden_pop(
#'   data = sim_data, observed = ~ m,
#'   auxiliary = ~ n, reference_pop = ~ N
#' )
#' summary(fit_pois)
#'
#' # OLS with fixed gamma
#' fit_ols <- estimate_hidden_pop(
#'   data = sim_data, observed = ~ m,
#'   auxiliary = ~ n, reference_pop = ~ N,
#'   method = "ols", gamma = 0.01
#' )
#' summary(fit_ols)
#'
#' # Constrained Poisson (alpha in (0,1), beta > 0)
#' fit_constr <- estimate_hidden_pop(
#'   data = sim_data, observed = ~ m,
#'   auxiliary = ~ n, reference_pop = ~ N,
#'   constrained = TRUE
#' )
#' summary(fit_constr)
#'
#' @export
estimate_hidden_pop <- function(data,
                                observed,
                                auxiliary,
                                reference_pop,
                                method = c("poisson", "nb", "ols", "nls", "iols"),
                                cov_alpha = NULL,
                                cov_beta = NULL,
                                gamma = "estimate",
                                gamma_bounds = c(1e-10, 0.5),
                                theta_start = 1,
                                vcov = "HC3",
                                weights = NULL,
                                constrained = FALSE,
                                countries = NULL,
                                cluster = NULL) {

  method <- match.arg(method)
  call <- match.call()
  ## Force-evaluate formula arguments so the stored call contains actual
  ## formula objects, not symbols. This ensures loo() can re-evaluate the
  ## call without needing the original calling environment (e.g., Shiny).
  call$observed <- observed
  call$auxiliary <- auxiliary
  call$reference_pop <- reference_pop
  if (!is.null(cov_alpha)) call$cov_alpha <- cov_alpha
  if (!is.null(cov_beta)) call$cov_beta <- cov_beta
  if (!is.null(countries)) call$countries <- countries
  if (!is.null(cluster)) call$cluster <- cluster

  ## Determine vcov type
  vcov_is_function <- is.function(vcov)
  vcov_type <- if (vcov_is_function) "HC0" else vcov

  # ---- Extract variables ----
  m <- eval(observed[[2]], data)
  n_aux <- eval(auxiliary[[2]], data)
  N <- eval(reference_pop[[2]], data)

  stopifnot(length(m) == length(N), length(m) == length(n_aux))
  stopifnot(all(N > 0))

  ratio <- n_aux / N
  log_N <- log(N)

  # ---- Determine gamma configuration ----
  estimate_gamma <- identical(gamma, "estimate")
  gamma_fixed <- is.numeric(gamma)
  no_gamma <- is.null(gamma)

  if (no_gamma) {
    if (any(n_aux == 0)) {
      warning("auxiliary contains zeros; log(n/N) will be -Inf. ",
              "Consider using gamma = 'estimate' or a fixed gamma value.")
    }
    log_rate <- log(ratio)
    gamma_value <- NULL
    gamma_start <- NULL
  } else if (gamma_fixed) {
    log_rate <- log(gamma + ratio)
    gamma_value <- gamma
    gamma_start <- NULL
  } else if (estimate_gamma) {
    gamma_start_val <- mean(gamma_bounds)
    log_rate <- log(gamma_start_val + ratio)
    gamma_value <- NULL
    gamma_start <- gamma_start_val
  }

  # ---- Build design matrices for alpha and beta ----
  X_alpha <- .build_model_matrix(cov_alpha, data)
  X_beta <- .build_model_matrix(cov_beta, data)

  # Name columns
  if (!is.null(cov_alpha)) {
    colnames(X_alpha) <- paste0("alpha:", colnames(X_alpha))
  } else {
    colnames(X_alpha) <- "alpha"
  }
  if (!is.null(cov_beta)) {
    colnames(X_beta) <- paste0("beta:", colnames(X_beta))
  } else {
    colnames(X_beta) <- "beta"
  }

  # ---- Extract cluster variable (if provided) ----
  cluster_vec <- if (!is.null(cluster)) {
    eval(cluster[[2]], data)
  } else NULL

  # ---- Dispatch to estimator ----
  # Estimators compute internal HC vcov (used as fallback);

  # final vcov is recomputed post-hoc using sandwich package.
  result <- switch(method,
    "ols" = {
      if (estimate_gamma) {
        .fit_ols_gamma(m, N, ratio, log_N, X_alpha, X_beta,
                       gamma_start, gamma_bounds,
                       weights = weights, vcov_type = vcov_type)
      } else {
        .fit_ols(m, N, ratio, log_N, log_rate, X_alpha, X_beta,
                 weights = weights, vcov_type = vcov_type)
      }
    },
    "nls" = {
      .fit_nls(m, N, ratio, log_N, log_rate, X_alpha, X_beta,
               gamma_value = gamma_value,
               estimate_gamma = estimate_gamma,
               gamma_start = gamma_start,
               gamma_bounds = gamma_bounds,
               weights = weights, vcov_type = vcov_type)
    },
    "poisson" = {
      .fit_poisson(m, N, ratio, log_N, log_rate, X_alpha, X_beta,
                   gamma_value = gamma_value,
                   estimate_gamma = estimate_gamma,
                   gamma_start = gamma_start,
                   gamma_bounds = gamma_bounds,
                   weights = weights, vcov_type = vcov_type,
                   constrained = constrained)
    },
    "nb" = {
      .fit_nb(m, N, ratio, log_N, log_rate, X_alpha, X_beta,
              gamma_value = gamma_value,
              estimate_gamma = estimate_gamma,
              gamma_start = gamma_start,
              gamma_bounds = gamma_bounds,
              theta_start = theta_start,
              weights = weights, vcov_type = vcov_type,
              constrained = constrained)
    },
    "iols" = {
      if (estimate_gamma) {
        .fit_iols_gamma(m, N, ratio, log_N, X_alpha, X_beta,
                        gamma_start, gamma_bounds,
                        weights = weights, vcov_type = vcov_type)
      } else {
        .fit_iols(m, N, ratio, log_N, log_rate, X_alpha, X_beta,
                  gamma_value = gamma_value,
                  weights = weights, vcov_type = vcov_type)
      }
    }
  )

  # ---- Build output object ----
  p_alpha <- ncol(X_alpha)
  p_beta <- ncol(X_beta)

  names(result$alpha_coefs) <- colnames(X_alpha)
  names(result$beta_coefs) <- colnames(X_beta)

  all_coefs <- c(result$alpha_coefs, result$beta_coefs)
  rownames(result$vcov) <- colnames(result$vcov) <- names(all_coefs)

  # Compute alpha/beta per observation (response scale for constrained)
  alpha_values <- if (!is.null(result$alpha_values)) {
    result$alpha_values
  } else {
    as.numeric(X_alpha %*% result$alpha_coefs)
  }
  beta_values <- if (!is.null(result$beta_values)) {
    result$beta_values
  } else {
    as.numeric(X_beta %*% result$beta_coefs)
  }

  # Warn if alpha out of sensible range
  if (any(alpha_values > 1, na.rm = TRUE)) {
    warning("Some alpha values > 1 (max = ", round(max(alpha_values), 3),
            "). Population size estimates may be unreliable. ",
            "Consider using constrained = TRUE or simplifying cov_alpha.")
  }
  if (any(alpha_values < 0, na.rm = TRUE)) {
    warning("Some alpha values < 0 (min = ", round(min(alpha_values), 3),
            "). Consider using constrained = TRUE.")
  }

  out <- list(
    coefficients = all_coefs,
    alpha_coefs = result$alpha_coefs,
    beta_coefs = result$beta_coefs,
    alpha_values = alpha_values,
    beta_values = beta_values,
    vcov = NULL,  # computed below via sandwich
    vcov_full = NULL,
    vcov_model = result$vcov_model,
    fitted.values = result$fitted,
    residuals = result$residuals,
    log_mu = result$log_mu,
    n_obs = result$n_obs,
    df.residual = result$df.residual,
    loglik = result$loglik,
    theta = result$theta,
    gamma = if (estimate_gamma) result$gamma_estimated
            else if (gamma_fixed) gamma
            else NULL,
    gamma_estimated = !is.null(result$gamma_estimated),
    constrained = constrained,
    method = method,
    vcov_type = vcov_type,
    vcov_spec = vcov,  # store original vcov specification
    call = call,
    data = data,
    m = m,
    N = N,
    n_aux = n_aux,
    X_alpha = X_alpha,
    X_beta = X_beta,
    p_alpha = p_alpha,
    p_beta = p_beta,
    ## Sandwich ingredients (stored for bread/estfun methods)
    model_matrix_full = result$model_matrix_full,
    bread_weights = result$bread_weights,
    score_residuals = result$score_residuals,
    cov_alpha_vars = if (!is.null(cov_alpha)) {
      data[, all.vars(cov_alpha), drop = FALSE]
    } else {
      NULL
    },
    countries_var = if (!is.null(countries)) eval(countries[[2]], data) else NULL,
    cluster_var = cluster_vec,
    obs_weights = weights,
    convergence = if (!is.null(result$convergence)) result$convergence else 0L,
    call_args = list(
      observed = observed, auxiliary = auxiliary, reference_pop = reference_pop
    )
  )

  # Store sigma2 for OLS/NLS
  if (method %in% c("ols", "nls") && !is.null(result$sigma2)) {
    out$sigma2 <- result$sigma2
  }

  class(out) <- "uncounted"

  # ---- Compute vcov ----
  p_ab <- p_alpha + p_beta

  if (method == "nb" && !is.null(result$hessian_nll) && is.character(vcov)) {
    # NB with character vcov: dedicated sandwich path including theta
    V_nb_full <- .compute_nb_sandwich(
      hessian_nll = result$hessian_nll,
      score_full  = result$score_full,
      n_obs       = nrow(data),
      vcov_type   = vcov,
      cluster     = cluster_vec
    )
    out$vcov <- V_nb_full[seq_len(p_ab), seq_len(p_ab), drop = FALSE]
    out$vcov_full <- V_nb_full
    out$theta_se <- sqrt(max(0, V_nb_full[p_ab + 1, p_ab + 1]))
    out$score_full <- result$score_full
    out$hessian_nll <- result$hessian_nll
    # Store both requested and actual vcov types.
    # NB theta-aware sandwich supports HC0 and HC1 only; HC2+ downgraded to HC1.
    actual_type <- if (vcov %in% c("HC2","HC3","HC4","HC4m","HC5")) "HC1" else vcov
    out$vcov_requested <- vcov
    out$vcov_type <- actual_type
  } else {
    # Non-NB, or NB with user-supplied vcov function: use sandwich package
    V_full <- .compute_vcov_sandwich(out, vcov, cluster_vec)
    has_gamma_col <- ncol(result$model_matrix_full) > p_ab
    if (has_gamma_col) {
      out$vcov_full <- V_full
      out$vcov <- V_full[seq_len(p_ab), seq_len(p_ab), drop = FALSE]
    } else {
      out$vcov_full <- V_full
      out$vcov <- V_full
    }
    # Store NB ingredients even when user vcov function is used
    if (method == "nb" && !is.null(result$hessian_nll)) {
      out$score_full <- result$score_full
      out$hessian_nll <- result$hessian_nll
    }
  }

  out
}


# ---- S3 methods ----

#' @export
print.uncounted <- function(x, ...) {
  cat("Unauthorized population estimation\n")
  vcov_label <- if (is.function(x$vcov_spec)) {
    "user-supplied function"
  } else if (!is.null(x$cluster_var)) {
    paste0("CL(", x$vcov_type, "), cluster: ", deparse(x$call$cluster))
  } else {
    x$vcov_type
  }
  cat("Method:", toupper(x$method), "| vcov:", vcov_label, "\n")
  cat("N obs:", x$n_obs, "\n")
  if (!is.null(x$gamma)) {
    cat("Gamma:", round(x$gamma, 6),
        if (x$gamma_estimated) "(estimated)" else "(fixed)", "\n")
  }
  if (!is.null(x$theta)) {
    cat("Theta (NB dispersion):", round(x$theta, 4), "\n")
  }
  if (!is.null(x$loglik)) {
    cat("Log-likelihood:", round(x$loglik, 2), "\n")
  }
  if (isTRUE(x$constrained)) {
    cat("Constrained: alpha in (0,1), beta > 0\n")
    cat("\nCoefficients (link scale: logit-alpha, log-beta):\n")
    print(round(x$coefficients, 6))
    cat("\nPopulation parameters (response scale):\n")
    .print_response_params(x)
  } else {
    cat("\nCoefficients:\n")
    print(round(x$coefficients, 6))
  }

  # Population size table
  cat("\n-----------------------\n")
  cat("Population size estimation results:\n")
  .print_popsize_table(popsize(x))

  invisible(x)
}

#' Summary of an Uncounted Population Model
#'
#' Prints a detailed summary of the fitted model, including coefficient
#' estimates with robust standard errors, z- or t-values, p-values, and
#' goodness-of-fit statistics (log-likelihood, AIC, BIC, deviance when
#' available). The output concludes with estimated population sizes
#' \eqn{\hat{\xi} = \sum_i N_i^{\hat{\alpha}_i}}{xi_hat = sum(N_i^alpha_hat_i)}
#' for each covariate group defined by \code{cov_alpha}, along with
#' bias-corrected estimates and confidence intervals computed by
#' \code{\link{popsize}}.
#'
#' For OLS and NLS fits, inference uses t-statistics with
#' \eqn{n - p} degrees of freedom. For Poisson and NB fits, z-statistics
#' (normal approximation) are used. Standard errors are always
#' HC-robust (type controlled by \code{vcov_type} in the original call).
#'
#' @param object An \code{"uncounted"} object returned by
#'   \code{\link{estimate_hidden_pop}}.
#' @param total Logical; if \code{TRUE}, include a total row in the
#'   population-size table. Default \code{FALSE}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the coefficient summary matrix (estimates,
#'   standard errors, test statistics, and p-values).
#'
#' @seealso \code{\link{estimate_hidden_pop}}, \code{\link{popsize}}
#'
#' @export
summary.uncounted <- function(object, total = FALSE, ...) {
  coefs <- object$coefficients
  se <- sqrt(diag(object$vcov))
  z <- coefs / se
  n <- object$n_obs
  p <- length(coefs)

  if (object$method %in% c("ols", "nls")) {
    df <- n - p
    pval <- 2 * pt(-abs(z), df = df)
  } else {
    pval <- 2 * pnorm(-abs(z))
  }

  tab <- cbind(
    Estimate = coefs,
    `Std. Error` = se,
    `z value` = z,
    `Pr(>|z|)` = pval
  )

  cat("Unauthorized population estimation\n")
  vcov_label <- if (is.function(object$vcov_spec)) {
    "user-supplied function"
  } else if (!is.null(object$cluster_var)) {
    paste0("CL(", object$vcov_type, "), cluster: ", deparse(object$call$cluster))
  } else {
    object$vcov_type
  }
  cat("Method:", toupper(object$method), "| vcov:", vcov_label, "\n")
  cat("N obs:", object$n_obs, "\n")
  if (!is.null(object$gamma)) {
    cat("Gamma:", round(object$gamma, 6),
        if (object$gamma_estimated) "(estimated)" else "(fixed)", "\n")
  }
  if (!is.null(object$theta)) {
    cat("Theta (NB dispersion):", round(object$theta, 4), "\n")
  }
  if (!is.null(object$loglik)) {
    cat("Log-likelihood:", round(object$loglik, 2), "\n")
  }
  if (!is.null(object$loglik)) {
    n_par <- .count_params(object)
    aic <- -2 * object$loglik + 2 * n_par
    bic <- -2 * object$loglik + log(object$n_obs) * n_par
    cat("AIC:", round(aic, 2), " BIC:", round(bic, 2), "\n")
  }
  dev <- tryCatch(deviance(object), error = function(e) NULL)
  if (!is.null(dev)) {
    cat("Deviance:", round(dev, 2), "\n")
  }
  if (isTRUE(object$constrained)) {
    cat("\nCoefficients (link scale: logit for alpha, log for beta):\n")
  } else {
    cat("\nCoefficients:\n")
  }
  printCoefmat(tab, P.values = TRUE, has.Pvalue = TRUE)

  # Response-scale summary for constrained models
  if (isTRUE(object$constrained)) {
    cat("\nResponse-scale parameters (alpha in (0,1), beta > 0):\n")
    .print_response_summary(object)
  }

  # Population size table
  cat("\n-----------------------\n")
  cat("Population size estimation results:\n")
  .print_popsize_table(popsize(object, total = total), total = total)

  invisible(tab)
}

#' @export
coef.uncounted <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.uncounted <- function(object, ...) {
  object$vcov
}

#' @export
fitted.uncounted <- function(object, ...) {
  object$fitted.values
}

#' @export
weights.uncounted <- function(object, ...) {
  object$obs_weights
}

#' Update and re-fit an uncounted model
#'
#' @param object An \code{uncounted} object.
#' @param ... Arguments to override in the original call (e.g.,
#'   \code{weights = w}, \code{vcov = "HC0"}).
#' @param evaluate Logical. If \code{TRUE} (default), evaluate the updated
#'   call and return the new fit. If \code{FALSE}, return the unevaluated call.
#' @return An \code{uncounted} object (if \code{evaluate = TRUE}) or
#'   an unevaluated call (if \code{evaluate = FALSE}).
#' @export
update.uncounted <- function(object, ..., evaluate = TRUE) {
  call <- object$call
  extras <- match.call(expand.dots = FALSE)$...
  for (nm in names(extras)) {
    call[[nm]] <- extras[[nm]]
  }
  if (!evaluate) return(call)
  eval(call, parent.frame())
}

#' Print response-scale alpha/beta per group for constrained models
#' @noRd
.print_response_params <- function(x) {
  .print_response_summary(x)
}

#' Print response-scale alpha/beta with delta-method SE for constrained models
#' @noRd
.print_response_summary <- function(x) {
  X_alpha <- x$X_alpha
  X_beta <- x$X_beta
  V <- x$vcov
  p_alpha <- x$p_alpha
  p_beta <- x$p_beta

  V_alpha <- V[seq_len(p_alpha), seq_len(p_alpha), drop = FALSE]
  V_beta <- V[p_alpha + seq_len(p_beta), p_alpha + seq_len(p_beta), drop = FALSE]

  # Unique alpha groups
  X_alpha_key <- apply(X_alpha, 1, paste, collapse = "|")
  unique_alpha_keys <- unique(X_alpha_key)

  # Unique beta groups
  X_beta_key <- apply(X_beta, 1, paste, collapse = "|")
  unique_beta_keys <- unique(X_beta_key)

  # Alpha rows
  alpha_rows <- lapply(unique_alpha_keys, function(key) {
    idx <- which(X_alpha_key == key)[1]
    x_row <- X_alpha[idx, , drop = FALSE]
    eta <- as.numeric(x_row %*% x$alpha_coefs)
    alpha_resp <- .inv_logit(eta)
    # Delta method: SE(alpha) = logit'(eta) * SE(eta)
    se_eta <- sqrt(as.numeric(x_row %*% V_alpha %*% t(x_row)))
    se_alpha <- .inv_logit_deriv(eta) * se_eta
    label <- .make_group_label(idx, x$cov_alpha_vars)
    data.frame(
      group = label,
      alpha = alpha_resp,
      `SE(alpha)` = se_alpha,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
  alpha_tab <- do.call(rbind, alpha_rows)
  rownames(alpha_tab) <- alpha_tab$group
  alpha_tab$group <- NULL

  cat("  Alpha (response scale):\n")
  print(round(alpha_tab, 4))

  # Beta rows (unique patterns)
  beta_rows <- lapply(unique_beta_keys, function(key) {
    idx <- which(X_beta_key == key)[1]
    z_row <- X_beta[idx, , drop = FALSE]
    zeta <- as.numeric(z_row %*% x$beta_coefs)
    beta_resp <- exp(zeta)
    # Delta method: SE(beta) = exp(zeta) * SE(zeta)
    se_zeta <- sqrt(as.numeric(z_row %*% V_beta %*% t(z_row)))
    se_beta <- beta_resp * se_zeta
    # Label from beta covariates
    label <- paste0("beta=", round(beta_resp, 4))
    if (!is.null(x$cov_alpha_vars)) {
      # Try to find a meaningful label — use the first obs with this pattern
      label <- paste0("(", idx, ")")
    }
    data.frame(
      group = key,
      beta = beta_resp,
      `SE(beta)` = se_beta,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
  beta_tab <- do.call(rbind, beta_rows)
  rownames(beta_tab) <- NULL
  cat("  Beta (response scale):\n")
  print(round(beta_tab[, c("beta", "SE(beta)"), drop = FALSE], 4))
}

#' Format population size table for printing
#' @noRd
.print_popsize_table <- function(ps, total = FALSE) {
  fmt <- function(x) {
    ifelse(is.finite(x) & abs(x) < 1e12,
           format(round(x), big.mark = ",", trim = TRUE),
           ifelse(is.finite(x),
                  formatC(x, format = "e", digits = 2),
                  "Inf"))
  }

  has_bc <- "estimate_bc" %in% names(ps) &&
            !all(ps$estimate_bc == ps$estimate)
  has_obs <- "observed" %in% names(ps)

  if (has_bc) {
    tab <- data.frame(
      Observed = if (has_obs) fmt(ps$observed) else NULL,
      Estimate = fmt(ps$estimate),
      `Estimate (BC)` = fmt(ps$estimate_bc),
      `CI lower` = fmt(ps$lower),
      `CI upper` = fmt(ps$upper),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    tab <- data.frame(
      Observed = if (has_obs) fmt(ps$observed) else NULL,
      Estimate = fmt(ps$estimate),
      `CI lower` = fmt(ps$lower),
      `CI upper` = fmt(ps$upper),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  rownames(tab) <- ps$group

  if (total && nrow(ps) > 1) {
    tot_info <- attr(ps, "total")
    tot_est <- if (!is.null(tot_info)) tot_info$estimate else sum(ps$estimate)
    tot_est_bc <- if (!is.null(tot_info)) tot_info$estimate_bc else sum(ps$estimate_bc)
    tot_lower <- if (!is.null(tot_info)) tot_info$lower else sum(ps$lower)
    tot_upper <- if (!is.null(tot_info)) tot_info$upper else sum(ps$upper)
    tot_row <- if (has_bc) {
      data.frame(
        Observed = if (has_obs) fmt(sum(ps$observed)) else NULL,
        Estimate = fmt(tot_est),
        `Estimate (BC)` = fmt(tot_est_bc),
        `CI lower` = fmt(tot_lower),
        `CI upper` = fmt(tot_upper),
        check.names = FALSE, stringsAsFactors = FALSE, row.names = "Total"
      )
    } else {
      data.frame(
        Observed = if (has_obs) fmt(sum(ps$observed)) else NULL,
        Estimate = fmt(tot_est),
        `CI lower` = fmt(tot_lower),
        `CI upper` = fmt(tot_upper),
        check.names = FALSE, stringsAsFactors = FALSE, row.names = "Total"
      )
    }
    tab <- rbind(tab, tot_row)
  }

  if (has_bc) {
    cat("  (BC = bias-corrected using model-based variance)\n")
  }

  # Warn if estimates look unreliable
  if (any(ps$upper / pmax(ps$estimate, 1) > 1e6, na.rm = TRUE)) {
    cat("  (!) Warning: some estimates have extremely wide CIs.\n")
    cat("      This suggests the model may be poorly identified.\n\n")
  }

  print(tab)
}
