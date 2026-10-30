#' Estimate the Size of an Unauthorized Migrant Population
#'
#' Fits a power-law model relating observed counts of unauthorized migrants
#' to reference population sizes and auxiliary registration counts. The model
#' is estimated using OLS, NLS, Poisson MLE, or Negative Binomial MLE, and
#' supports covariate-varying parameters through interaction formulas.
#'
#' @param data A data frame containing all variables.
#' @param observed One-sided formula for the observed count (e.g., \code{~ m}).
#'   Must be whole numbers for \code{method = "nb"}. \code{method = "poisson"}
#'   also accepts non-integer values (pseudo-maximum likelihood); \code{$loglik}
#'   is then a pseudo-log-likelihood, flagged by \code{$pseudo_loglik}.
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
#' @param cov_gamma Optional formula for covariates in gamma (e.g.,
#'   \code{~ sex} or \code{~ year}). When specified, gamma varies across
#'   observations via \eqn{\gamma_i = \exp(\mathbf{X}_{\gamma,i}' \boldsymbol{\delta})}{gamma_i = exp(X_gamma_i' delta)},
#'   where \eqn{\boldsymbol{\delta}}{delta} are estimated coefficients on the
#'   log scale. Requires \code{gamma = "estimate"} and
#'   \code{method} \code{"poisson"} or \code{"nb"}.
#'   Default \code{NULL} means a single scalar gamma (intercept only).
#' @param gamma_bounds Numeric vector of length 2: lower and upper bounds for
#'   the estimated gamma parameter. Default \code{c(1e-10, 0.5)}.
#' @param gamma_start Optional starting value for the estimated gamma (the
#'   gamma intercept when \code{cov_gamma} is used); must lie inside
#'   \code{gamma_bounds}. Default \code{NULL} screens 25 log-spaced values in
#'   \code{gamma_bounds} and starts from the one with the smallest Poisson
#'   profile deviance. Ignored unless \code{gamma = "estimate"}, and by
#'   \code{method = "ols"}, which profiles gamma over its own grid. Refits
#'   inside \code{bootstrap_popsize()}, \code{loo()} and the profile functions
#'   start from the original estimate. The screen and the optimiser's trace
#'   are stored in \code{$optim_info}.
#' @param theta_start Starting value for the NB dispersion parameter (used
#'   only when \code{method = "nb"}). Default 1.
#' @param link_rho Link for the detection component
#'   \eqn{\rho_i = h(\beta_i \log(\gamma_i + n_i / N_i))}. Supported values are
#'   \code{"power"} (default, \eqn{\rho_i = \exp(\eta_i)}),
#'   \code{"cloglog"} (\eqn{\rho_i = 1 - \exp(-\exp(\eta_i))}), and
#'   \code{"logit"} (\eqn{\rho_i = 1 / (1 + \exp(-\eta_i))}), and
#'   \code{"probit"} (\eqn{\rho_i = \Phi(\eta_i)}).
#'   The legacy alias \code{"logistic"} is accepted and normalized to
#'   \code{"logit"}.
#'   Bounded links are currently available for \code{method = "poisson"},
#'   \code{"nb"}, and \code{"nls"}.
#' @param estimator Estimation criterion for count models:
#'   \code{"mle"} (default), \code{"gmm"}, or \code{"el"}.
#'   Non-MLE estimators are currently available only for
#'   \code{method = "poisson"} and \code{method = "nb"}.
#' @param vcov Controls the variance-covariance estimator. Can be:
#'   \itemize{
#'     \item A character string specifying the HC type: \code{"HC0"} through
#'       \code{"HC5"}, or \code{"HC4m"}. Default \code{"HC3"}. When \code{cluster}
#'       is also provided, clustered HC2+ requests are downgraded to
#'       \code{"HC1"} before calling \code{sandwich::vcovCL()}.
#'       For count models fitted with \code{estimator = "gmm"} or
#'       \code{estimator = "el"}, the default robust covariance is
#'       \code{"HC1"} and HC2+ requests are downgraded to \code{"HC1"}.
#'     \item A function that takes the fitted \code{uncounted} object and returns
#'       a variance-covariance matrix. For example,
#'       \code{sandwich::vcovHC} or
#'       \code{function(x) sandwich::vcovCL(x, cluster = data$country_code)}.
#'   }
#' @param cluster Optional one-sided formula identifying a cluster variable
#'   for cluster-robust variance estimation (e.g., \code{~ country_code}).
#'   When provided and \code{vcov} is a character string, the variance is
#'   computed using \code{sandwich::vcovCL()} with the requested HC type, except
#'   that HC2+ requests are reported in \code{vcov_requested} and computed as
#'   \code{"HC1"}. Ignored when \code{vcov} is a function.
#' @param weights Optional numeric vector of observation weights.
#' @param constrained Logical. If \code{TRUE} and
#'   \code{method \%in\% c("poisson", "nb")}, the fitted per-observation
#'   parameters are constrained on the response scale via
#'   \eqn{\alpha_i = \mathrm{logit}^{-1}(\mathbf{X}_{\alpha,i}' \mathbf{a})}
#'   and \eqn{\beta_i = \exp(\mathbf{X}_{\beta,i}' \mathbf{b})}.
#'   This keeps \eqn{\alpha_i \in (0, 1)} and \eqn{\beta_i > 0}. The reported
#'   coefficient vectors \code{alpha_coefs} and \code{beta_coefs} remain on
#'   the transformed scales (logit for alpha, log for beta); use
#'   \code{summary()} or \code{alpha_values} / \code{beta_values} for
#'   response-scale values. Default \code{FALSE}. See Details.
#' @param countries One-sided formula identifying a country or group variable
#'   used when reporting population size estimates (e.g., \code{~ country}).
#'   This does \strong{not} enable cluster-robust variance; use \code{cluster}
#'   for that.
#'
#' @details
#' \strong{Model specification.}
#' For observation \eqn{i}, the expected observed count \eqn{m_i} is modelled as:
#'
#' \deqn{E(m_i) = \xi_i \rho_i}{E(m_i) = xi_i rho_i}
#'
#' where \eqn{\xi_i = E(M_i \mid N_i)}{xi_i = E(M_i | N_i)} is the theoretical
#' unauthorized population size and
#' \eqn{\rho_i = E(p_i \mid N_i, n_i)}{rho_i = E(p_i | N_i, n_i)} is the
#' theoretical detection rate. In the baseline empirical specification,
#' \eqn{\xi_i = N_i^{\alpha_i}}{xi_i = N_i^alpha_i} and, under the default
#' \code{link_rho = "power"},
#' \eqn{\rho_i = (\gamma_i + n_i / N_i)^{\beta_i}}{rho_i = (gamma_i + n_i / N_i)^beta_i}.
#' More generally, the package writes the detection component as
#' \deqn{\eta_i = \beta_i \log(\gamma_i + n_i / N_i), \qquad \rho_i = h(\eta_i)}{eta_i = beta_i log(gamma_i + n_i / N_i), rho_i = h(eta_i)}
#' where \eqn{h}{h} is chosen by \code{link_rho}. Supported links are:
#' \describe{
#'   \item{\code{"power"}}{\eqn{\rho_i = \exp(\eta_i)}{rho_i = exp(eta_i)}}
#'   \item{\code{"cloglog"}}{\eqn{\rho_i = 1 - \exp(-\exp(\eta_i))}{rho_i = 1 - exp(-exp(eta_i))}}
#'   \item{\code{"logit"}}{\eqn{\rho_i = 1 / (1 + \exp(-\eta_i))}{rho_i = 1 / (1 + exp(-eta_i))}}
#'   \item{\code{"probit"}}{\eqn{\rho_i = \Phi(\eta_i)}{rho_i = Phi(eta_i)}}
#' }
#' The \code{"power"} link reproduces the paper's baseline power-law detection
#' model exactly; the bounded links are package extensions. On the log scale,
#'
#' \deqn{\log E(m_i) = \alpha_i \log N_i + \log(\rho_i)}{log E(m_i) = alpha_i * log(N_i) + log(rho_i)}
#'
#' so for every supported \code{link_rho}, \eqn{\xi_i = E(m_i) / \rho_i = N_i^{\alpha_i}}{xi_i = E(m_i) / rho_i = N_i^alpha_i}.
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
#'     uncertainty in \eqn{\theta} estimation. For \code{method = "poisson"}
#'     and \code{method = "nb"}, the same moment conditions can also be fit
#'     with \code{estimator = "gmm"} or \code{estimator = "el"} using
#'     \pkg{momentfit}.}
#' }
#'
#' \strong{Constrained vs. unconstrained estimation.}
#' When \code{constrained = FALSE} (default), \eqn{\alpha} and \eqn{\beta}
#' are estimated on the real line without restrictions. When
#' \code{constrained = TRUE}, the linear predictors are
#' \eqn{\eta_{\alpha,i} = \mathbf{X}_{\alpha,i}' \mathbf{a}}{eta_alpha_i = X_alpha_i' a}
#' and \eqn{\eta_{\beta,i} = \mathbf{X}_{\beta,i}' \mathbf{b}}{eta_beta_i = X_beta_i' b},
#' with response-scale parameters
#' \eqn{\alpha_i = \mathrm{logit}^{-1}(\eta_{\alpha,i})}{alpha_i = logit^-1(eta_alpha_i)}
#' and \eqn{\beta_i = \exp(\eta_{\beta,i})}{beta_i = exp(eta_beta_i)}.
#' This constrains the fitted \eqn{\alpha_i}{alpha_i} values to \eqn{(0, 1)}
#' and the fitted \eqn{\beta_i}{beta_i} values to be strictly positive.
#' The stored coefficient vectors \code{alpha_coefs} and \code{beta_coefs}
#' remain on the transformed scales (logit for alpha, log for beta); use
#' \code{summary()} or the \code{alpha_values} / \code{beta_values}
#' elements of the returned object for response-scale values.
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
#'       matrix, used for bias correction. For Poisson and NB count models,
#'       this is the Fisher-style inverse information for the mean-model
#'       parameters, evaluated at the fitted coefficients even when
#'       \code{estimator = "gmm"} or \code{estimator = "el"}.}
#'     \item{\code{rho_values}}{Per-observation fitted detection component
#'       \eqn{\hat\rho_i}{rho_hat_i}.}
#'     \item{\code{fitted.values}}{Fitted values \eqn{\hat{m}_i}{m_hat_i}.}
#'     \item{\code{residuals}}{Raw residuals \eqn{m_i - \hat{m}_i}{m_i - m_hat_i}.}
#'     \item{\code{gamma}}{Estimated or fixed gamma value (or \code{NULL}).}
#'     \item{\code{theta}}{NB dispersion parameter (only for \code{method = "nb"}).}
#'     \item{\code{loglik}}{Log-likelihood (for MLE count fits only).}
#'     \item{\code{method}}{The model family used.}
#'     \item{\code{estimator}}{The estimation criterion used.}
#'     \item{\code{link_rho}}{The detection-link transformation used.}
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
#' \emph{The Review of Economics and Statistics}, 88(4), 641--658. \doi{10.1162/rest.88.4.641}
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
                                cov_gamma = NULL,
                                gamma_bounds = c(1e-10, 0.5),
                                gamma_start = NULL,
                                theta_start = 1,
                                link_rho = c("power", "cloglog", "logit", "probit"),
                                estimator = c("mle", "gmm", "el"),
                                vcov = "HC3",
                                weights = NULL,
                                constrained = FALSE,
                                countries = NULL,
                                cluster = NULL) {

  method <- match.arg(method)
  link_rho <- .normalize_link_rho(link_rho)
  estimator <- match.arg(estimator)
  vcov_missing <- missing(vcov)
  call <- match.call()
  ## Force-evaluate formula arguments so the stored call contains actual
  ## formula objects, not symbols. This ensures loo() can re-evaluate the
  ## call without needing the original calling environment (e.g., Shiny).
  call$data <- data
  call$observed <- observed
  call$auxiliary <- auxiliary
  call$reference_pop <- reference_pop
  call$method <- method
  call$gamma <- gamma
  call$gamma_bounds <- gamma_bounds
  if (!is.null(gamma_start)) call$gamma_start <- gamma_start
  call$theta_start <- theta_start
  call$constrained <- constrained
  if (!is.null(cov_alpha)) call$cov_alpha <- cov_alpha
  if (!is.null(cov_beta)) call$cov_beta <- cov_beta
  if (!is.null(cov_gamma)) call$cov_gamma <- cov_gamma
  if (!is.null(weights)) call$weights <- weights
  if (!is.null(countries)) call$countries <- countries
  if (!is.null(cluster)) call$cluster <- cluster
  call$link_rho <- link_rho
  call$estimator <- estimator

  ## Determine vcov type
  vcov_is_function <- is.function(vcov)
  vcov_type <- if (vcov_is_function) "HC0" else vcov
  if (estimator != "mle" && !method %in% c("poisson", "nb")) {
    stop("'estimator' is only supported for method = 'poisson' and method = 'nb'.",
         call. = FALSE)
  }
  if (constrained && !method %in% c("poisson", "nb")) {
    stop("'constrained = TRUE' is only available for method = 'poisson' and ",
         "method = 'nb'.", call. = FALSE)
  }
  if (link_rho != "power" && method %in% c("ols", "iols")) {
    stop("'link_rho' values other than 'power' are currently available only ",
         "for method = 'poisson', method = 'nb', and method = 'nls'.",
         call. = FALSE)
  }
  actual_vcov_type <- if (estimator != "mle" && !vcov_is_function) {
    if (vcov_missing && identical(vcov_type, "HC3")) {
      "HC1"
    } else {
      .normalize_moment_vcov_type(vcov_type)
    }
  } else {
    vcov_type
  }

  # ---- Extract variables ----
  m <- eval(observed[[2]], data)
  n_aux <- eval(auxiliary[[2]], data)
  N <- eval(reference_pop[[2]], data)

  # ---- Validate the response and the two count inputs ----
  .check_count <- function(x, label, formula) {
    if (!is.numeric(x)) {
      stop("'", label, "' (", deparse(formula), ") must be numeric, not ",
           class(x)[1], ".", call. = FALSE)
    }
    if (anyNA(x)) {
      stop("'", label, "' (", deparse(formula), ") contains ", sum(is.na(x)),
           " missing value(s). Remove or impute them before fitting.",
           call. = FALSE)
    }
    if (any(!is.finite(x))) {
      stop("'", label, "' (", deparse(formula), ") contains non-finite values.",
           call. = FALSE)
    }
    if (any(x < 0)) {
      stop("'", label, "' (", deparse(formula), ") contains negative values; ",
           "counts and population sizes must be non-negative.", call. = FALSE)
    }
    invisible(TRUE)
  }
  .check_count(m, "observed", observed)
  .check_count(n_aux, "auxiliary", auxiliary)
  .check_count(N, "reference_pop", reference_pop)

  if (length(m) != length(N) || length(m) != length(n_aux)) {
    stop("'observed', 'auxiliary' and 'reference_pop' must have the same ",
         "length (got ", length(m), ", ", length(n_aux), ", ", length(N),
         ").", call. = FALSE)
  }
  if (any(N <= 0)) {
    stop("'reference_pop' must be strictly positive; ", sum(N <= 0),
         " observation(s) are zero or negative.", call. = FALSE)
  }
  if (any(N < 1)) {
    warning("Some reference populations are below 1, so log(N) < 0 and the ",
            "corresponding N^alpha contributions shrink as alpha grows. ",
            "The population-size estimator assumes N_i >= 1.", call. = FALSE)
  }

  # ---- Validate gamma and its bounds ----
  if (!is.null(gamma)) {
    if (is.character(gamma)) {
      if (!identical(gamma, "estimate")) {
        stop("'gamma' must be \"estimate\", a single non-negative number, or ",
             "NULL; got \"", gamma[1], "\".", call. = FALSE)
      }
    } else if (is.numeric(gamma)) {
      if (length(gamma) != 1L || !is.finite(gamma) || gamma < 0) {
        stop("A numeric 'gamma' must be a single finite non-negative value.",
             call. = FALSE)
      }
    } else {
      stop("'gamma' must be \"estimate\", a single non-negative number, or ",
           "NULL; got an object of class ", class(gamma)[1], ".",
           call. = FALSE)
    }
  }
  if (identical(gamma, "estimate")) {
    if (!is.numeric(gamma_bounds) || length(gamma_bounds) != 2L ||
        any(!is.finite(gamma_bounds))) {
      stop("'gamma_bounds' must be a numeric vector of length 2 with finite ",
           "values.", call. = FALSE)
    }
    if (gamma_bounds[1] < 0 || gamma_bounds[2] <= gamma_bounds[1]) {
      stop("'gamma_bounds' must satisfy 0 <= lower < upper; got c(",
           gamma_bounds[1], ", ", gamma_bounds[2], ").", call. = FALSE)
    }
  }

  # ---- Validate weights ----
  if (!is.null(weights)) {
    if (!is.numeric(weights)) {
      stop("'weights' must be numeric, not ", class(weights)[1], ".",
           call. = FALSE)
    }
    if (length(weights) != length(m)) {
      stop("'weights' must have one value per observation (got ",
           length(weights), " for ", length(m), " observations).",
           call. = FALSE)
    }
    if (anyNA(weights) || any(!is.finite(weights))) {
      stop("'weights' must not contain missing or non-finite values.",
           call. = FALSE)
    }
    if (any(weights < 0)) {
      stop("'weights' must be non-negative.", call. = FALSE)
    }
    if (all(weights == 0)) {
      stop("'weights' are all zero, so no observation contributes to the fit.",
           call. = FALSE)
    }
  }

  non_integer_m <- any(abs(m - round(m)) > sqrt(.Machine$double.eps))
  if (identical(method, "nb") && non_integer_m) {
    stop("method = 'nb' requires whole-number counts in 'observed'; the negative ",
         "binomial likelihood is not defined for non-integer values. ",
         "method = 'poisson' accepts them (pseudo-maximum likelihood).",
         call. = FALSE)
  }

  if (identical(method, "nb")) {
    if (!is.numeric(theta_start) || length(theta_start) != 1L ||
        !is.finite(theta_start) || theta_start <= 0) {
      stop("'theta_start' must be a single positive finite value.",
           call. = FALSE)
    }
  }

  ratio <- n_aux / N
  log_N <- log(N)

  # ---- Determine gamma configuration ----
  estimate_gamma <- identical(gamma, "estimate")
  gamma_fixed <- is.numeric(gamma)
  no_gamma <- is.null(gamma)

  if (no_gamma) {
    if (any(n_aux == 0)) {
      stop("'gamma = NULL' uses log(n/N), which is undefined for the ",
           sum(n_aux == 0), " observation(s) where the auxiliary count is ",
           "zero. Use gamma = \"estimate\" or a fixed positive gamma.",
           call. = FALSE)
    }
    log_rate <- log(ratio)
    gamma_value <- NULL
    gamma_start <- NULL
  } else if (gamma_fixed) {
    if (gamma == 0 && any(n_aux == 0)) {
      stop("'gamma = 0' uses log(n/N), which is undefined for the ",
           sum(n_aux == 0), " observation(s) where the auxiliary count is ",
           "zero. Use gamma = \"estimate\" or a fixed positive gamma.",
           call. = FALSE)
    }
    log_rate <- log(gamma + ratio)
    gamma_value <- gamma
    gamma_start <- NULL
  } else if (estimate_gamma) {
    ## log_rate is set once the design matrices exist (start screen below)
    gamma_value <- NULL
    if (!is.null(gamma_start)) {
      if (!is.numeric(gamma_start) || length(gamma_start) != 1L ||
          !is.finite(gamma_start) || gamma_start < gamma_bounds[1] ||
          gamma_start > gamma_bounds[2]) {
        stop("'gamma_start' must be a single value inside gamma_bounds = c(",
             gamma_bounds[1], ", ", gamma_bounds[2], ").", call. = FALSE)
      }
    }
  }

  # ---- Validate cov_gamma ----
  if (!is.null(cov_gamma)) {
    if (!estimate_gamma) {
      stop("'cov_gamma' requires gamma = 'estimate'. Cannot use cov_gamma with ",
           "fixed or NULL gamma.", call. = FALSE)
    }
    if (!method %in% c("poisson", "nb")) {
      stop("'cov_gamma' is only supported for method = 'poisson' and method = 'nb'.",
           call. = FALSE)
    }
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

  # Build design matrix for gamma covariates
  X_gamma <- if (!is.null(cov_gamma)) {
    Xg <- .build_model_matrix(cov_gamma, data)
    colnames(Xg) <- paste0("gamma:", colnames(Xg))
    Xg
  } else if (estimate_gamma) {
    matrix(1, nrow = nrow(data), ncol = 1,
           dimnames = list(NULL, "gamma"))
  } else {
    NULL
  }
  p_gamma <- if (!is.null(X_gamma)) ncol(X_gamma) else 0L

  # ---- Starting value for an estimated gamma ----
  gamma_screen <- NULL
  if (estimate_gamma && method != "ols") {
    if (is.null(gamma_start)) {
      gamma_screen <- .screen_gamma_start(m, log_N, ratio, X_alpha, X_beta,
                                          gamma_bounds, weights)
      gamma_start <- gamma_screen$start
    }
    log_rate <- log(gamma_start + ratio)
  }

  # ---- Extract cluster variable (if provided) ----
  cluster_vec <- if (!is.null(cluster)) {
    eval(cluster[[2]], data)
  } else NULL
  if (!is.null(cluster_vec)) {
    if (anyNA(cluster_vec)) {
      stop("'cluster' (", deparse(cluster), ") contains ", sum(is.na(cluster_vec)),
           " missing value(s). Remove or impute them before fitting.",
           call. = FALSE)
    }
    n_clusters <- length(unique(cluster_vec[!is.na(cluster_vec)]))
    if (n_clusters < 2L) {
      stop("Cluster-robust variance requires at least 2 clusters.",
           call. = FALSE)
    }
    if (!vcov_is_function) {
      actual_vcov_type <- .normalize_cluster_vcov_type(actual_vcov_type)
    }
  }

  # ---- Dispatch to estimator ----
  # Estimators compute internal HC vcov (used as fallback);

  # final vcov is recomputed post-hoc using sandwich package.
  result <- switch(method,
    "ols" = {
      if (estimate_gamma) {
        .fit_ols_gamma(m, N, ratio, log_N, X_alpha, X_beta,
                       gamma_start, gamma_bounds,
                       weights = weights)
      } else {
        .fit_ols(m, N, ratio, log_N, log_rate, X_alpha, X_beta,
                 weights = weights)
      }
    },
    "nls" = {
      .fit_nls(m, N, ratio, log_N, log_rate, X_alpha, X_beta,
               gamma_value = gamma_value,
               estimate_gamma = estimate_gamma,
               gamma_start = gamma_start,
               gamma_bounds = gamma_bounds,
               weights = weights,
               link_rho = link_rho)
    },
    "poisson" = {
      if (estimator == "mle") {
        .fit_poisson(m, N, ratio, log_N, log_rate, X_alpha, X_beta,
                     gamma_value = gamma_value,
                     estimate_gamma = estimate_gamma,
                     gamma_start = gamma_start,
                     gamma_bounds = gamma_bounds,
                     X_gamma = X_gamma,
                     weights = weights,
                     constrained = constrained,
                     link_rho = link_rho)
      } else {
        .fit_poisson_moment(m, N, ratio, log_N, log_rate, X_alpha, X_beta,
                            gamma_value = gamma_value,
                            estimate_gamma = estimate_gamma,
                            gamma_start = gamma_start,
                            gamma_bounds = gamma_bounds,
                            X_gamma = X_gamma,
                            weights = weights,
                            constrained = constrained,
                            link_rho = link_rho,
                            estimator = estimator)
      }
    },
    "nb" = {
      if (estimator == "mle") {
        .fit_nb(m, N, ratio, log_N, log_rate, X_alpha, X_beta,
                gamma_value = gamma_value,
                estimate_gamma = estimate_gamma,
                gamma_start = gamma_start,
                gamma_bounds = gamma_bounds,
                theta_start = theta_start,
                X_gamma = X_gamma,
                weights = weights,
                constrained = constrained,
                link_rho = link_rho)
      } else {
        .fit_nb_moment(m, N, ratio, log_N, log_rate, X_alpha, X_beta,
                       gamma_value = gamma_value,
                       estimate_gamma = estimate_gamma,
                       gamma_start = gamma_start,
                       gamma_bounds = gamma_bounds,
                       theta_start = theta_start,
                       X_gamma = X_gamma,
                       weights = weights,
                       constrained = constrained,
                       link_rho = link_rho,
                       estimator = estimator)
      }
    },
    "iols" = {
      if (estimate_gamma) {
        stop("gamma = 'estimate' is not yet supported for method = 'iols'. ",
             "Use a fixed gamma value instead.", call. = FALSE)
      }
      .fit_iols(m, N, ratio, log_N, log_rate, X_alpha, X_beta,
                gamma_value = gamma_value,
                weights = weights)
    }
  )

  # ---- Build output object ----
  p_alpha <- ncol(X_alpha)
  p_beta <- ncol(X_beta)

  names(result$alpha_coefs) <- colnames(X_alpha)
  names(result$beta_coefs) <- colnames(X_beta)

  # Track whether cov_gamma was user-specified (vs. the implicit scalar design)
  has_cov_gamma <- !is.null(cov_gamma)

  # When cov_gamma specified, include gamma coefs in all_coefs
  all_coefs <- if (has_cov_gamma && !is.null(result$gamma_coefs)) {
    c(result$alpha_coefs, result$beta_coefs, result$gamma_coefs)
  } else {
    c(result$alpha_coefs, result$beta_coefs)
  }
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
  gamma_obs <- if (!is.null(result$gamma_values)) {
    result$gamma_values
  } else if (estimate_gamma && !is.null(result$gamma_estimated)) {
    rep(result$gamma_estimated, length(ratio))
  } else if (gamma_fixed) {
    rep(gamma, length(ratio))
  } else {
    NULL
  }
  rate_values <- .rate_from_gamma(ratio, gamma_obs)
  rho_values <- .rho_from_eta(.eta_from_rate(beta_values, rate_values),
                              link_rho = link_rho)

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

  # ---- Optimiser trace and gamma-bound check ----
  optim_info <- NULL
  if (estimate_gamma || !is.null(result$optim_info)) {
    optim_info <- c(
      list(gamma_start = if (estimate_gamma) gamma_start else NULL,
           gamma_start_source = if (!estimate_gamma) NULL
                                else if (method == "ols") "ols_grid"
                                else if (is.null(gamma_screen)) "user"
                                else "screen",
           gamma_screen = gamma_screen),
      result$optim_info
    )
  }
  if (!is.null(result$optim_info) && isTRUE(result$optim_info$max_grad >
                                             1e-4 * (1 + abs(result$loglik)))) {
    warning("The optimiser stopped with a non-negligible gradient (max |grad| = ",
            signif(result$optim_info$max_grad, 3), "); the estimates may not be ",
            "a maximum. Try another 'gamma_start' or narrower 'gamma_bounds'. ",
            "See $optim_info.", call. = FALSE)
  }
  if (estimate_gamma) {
    gamma_level <- if (has_cov_gamma && !is.null(result$gamma_coefs)) {
      exp(unname(result$gamma_coefs[1]))
    } else {
      result$gamma_estimated
    }
    if (!is.null(gamma_level) && is.finite(gamma_level)) {
      at_bound <- c(lower = abs(log(gamma_level) - log(gamma_bounds[1])) < 1e-4,
                    upper = abs(log(gamma_level) - log(gamma_bounds[2])) < 1e-4)
      if (any(at_bound)) {
        message("Estimated gamma ", if (has_cov_gamma) "(intercept) ",
                "is at its ", names(at_bound)[at_bound], " bound (",
                signif(gamma_level, 3), "); the optimum may lie outside ",
                "gamma_bounds = c(", gamma_bounds[1], ", ", gamma_bounds[2], ").")
      }
      optim_info$gamma_at_bound <- if (any(at_bound)) names(at_bound)[at_bound] else "none"
    }
  }

  out <- list(
    coefficients = all_coefs,
    alpha_coefs = result$alpha_coefs,
    beta_coefs = result$beta_coefs,
    alpha_values = alpha_values,
    beta_values = beta_values,
    rho_values = rho_values,
    vcov = NULL,  # computed below via sandwich
    vcov_full = NULL,
    vcov_model = result$vcov_model,
    fitted.values = result$fitted,
    residuals = result$residuals,
    log_mu = result$log_mu,
    ## Offset in force when this fit was produced (see profile_dependence()).
    ## Stored so that anything recomputing log(mu) later, such as predict(),
    ## reproduces the fitted values instead of reading whatever the global
    ## state happens to be at call time.
    dependence_offset = .current_dependence_offset(),
    n_obs = result$n_obs,
    df.residual = result$df.residual,
    loglik = result$loglik,
    theta = result$theta,
    gamma = if (has_cov_gamma) NULL  # varying gamma: no scalar slot
            else if (estimate_gamma) result$gamma_estimated
            else if (gamma_fixed) gamma
            else NULL,
    gamma_estimated = !is.null(result$gamma_estimated),
    gamma_coefs = result$gamma_coefs,
    gamma_values = result$gamma_values,
    X_gamma = X_gamma,
    p_gamma = p_gamma,
    has_cov_gamma = has_cov_gamma,
    constrained = constrained,
    method = method,
    estimator = estimator,
    link_rho = link_rho,
    vcov_type = actual_vcov_type,
    vcov_requested = if (vcov_is_function) NULL else vcov_type,
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
    theta_se = result$theta_se,
    score_full = result$score_full,
    hessian_nll = result$hessian_nll,
    backend_vcov_full = result$backend_vcov_full,
    backend_vcov_model_full = result$backend_vcov_model_full,
    moment_backend = result$moment_backend,
    optim_info = optim_info,
    pseudo_loglik = non_integer_m && method == "poisson",
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
  # Number of first-class coefs (alpha+beta, or alpha+beta+gamma when cov_gamma)
  p_coefs <- length(all_coefs)

  # Helper: extract the appropriate vcov submatrix.
  # When cov_gamma: vcov includes alpha+beta+gamma = p_coefs
  # When scalar gamma (p_gamma <= 1): vcov includes only alpha+beta = p_ab
  .extract_vcov <- function(V_full) {
    n_total <- ncol(V_full)
    if (n_total == p_coefs) {
      list(vcov_full = V_full, vcov = V_full)
    } else if (n_total > p_coefs) {
      list(vcov_full = V_full,
           vcov = V_full[seq_len(p_coefs), seq_len(p_coefs), drop = FALSE])
    } else {
      list(vcov_full = V_full, vcov = V_full)
    }
  }

  if (method == "nb" && estimator == "mle" &&
      !is.null(result$hessian_nll) && is.character(vcov)) {
    # NB with character vcov: dedicated sandwich path including theta
    V_nb_full <- .compute_nb_sandwich(
      hessian_nll = result$hessian_nll,
      score_full  = result$score_full,
      n_obs       = nrow(data),
      vcov_type   = actual_vcov_type,
      cluster     = cluster_vec
    )
    ## V_nb_full is ordered (alpha, beta, log-theta, gamma...): the reported
    ## coefficients skip theta, unlike the other methods' vcov_full.
    coef_idx <- c(seq_len(p_ab), if (has_cov_gamma) p_ab + 1L + seq_len(p_gamma))
    out$vcov <- V_nb_full[coef_idx, coef_idx, drop = FALSE]
    dimnames(out$vcov) <- list(names(all_coefs), names(all_coefs))
    out$vcov_full <- V_nb_full
    ## The optimiser works on log(theta): delta method to the theta scale
    out$theta_se <- as.numeric(out$theta) * sqrt(max(0, V_nb_full[p_ab + 1, p_ab + 1]))
    out$score_full <- result$score_full
    out$hessian_nll <- result$hessian_nll
    # Store both requested and actual vcov types.
    # NB theta-aware sandwich supports HC0 and HC1 only; HC2+ downgraded to HC1.
    actual_type <- if (vcov %in% c("HC2","HC3","HC4","HC4m","HC5")) "HC1" else actual_vcov_type
    out$vcov_requested <- vcov_type
    out$vcov_type <- actual_type
  } else if (is.function(vcov)) {
    # User-supplied vcov function: first build object with default HC vcov
    # (so update() inside vcovFWB etc. doesn't recurse), then apply function.
    default_type <- if (estimator != "mle" && method %in% c("poisson", "nb")) {
      "HC1"
    } else {
      "HC3"
    }
    V_default <- .compute_vcov_sandwich(out, default_type, cluster_vec)
    vv <- .extract_vcov(V_default)
    out$vcov_full <- vv$vcov_full
    out$vcov <- vv$vcov
    # Replace vcov in the stored call so update() doesn't re-invoke the function
    out$vcov_spec <- default_type
    out$call[["vcov"]] <- default_type
    # Now apply the user function on the fully constructed object
    V_user <- vcov(out)
    vv <- .extract_vcov(V_user)
    out$vcov_full <- vv$vcov_full
    out$vcov <- vv$vcov
    # Restore the original vcov function for printing and update()
    out$vcov_spec <- vcov
    out$call[["vcov"]] <- vcov
    # Store NB ingredients
    if (method == "nb" && !is.null(result$hessian_nll)) {
      out$score_full <- result$score_full
      out$hessian_nll <- result$hessian_nll
    }
  } else {
    # Character vcov type: use sandwich package directly
    V_full <- .compute_vcov_sandwich(out, actual_vcov_type, cluster_vec)
    vv <- .extract_vcov(V_full)
    out$vcov_full <- vv$vcov_full
    out$vcov <- vv$vcov
    # Store NB ingredients
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
  cat("Method:", toupper(x$method), "| estimator:", toupper(x$estimator),
      "| link_rho:", x$link_rho, "| vcov:", vcov_label, "\n")
  cat("N obs:", x$n_obs, "\n")
  if (!is.null(x$gamma_values) && length(x$gamma_values) > 0 && isTRUE(x$has_cov_gamma)) {
    cat("Gamma: covariate-varying (response scale), range [",
        round(min(x$gamma_values), 6), ", ",
        round(max(x$gamma_values), 6), "]\n", sep = "")
  } else if (!is.null(x$gamma)) {
    cat("Gamma:", round(x$gamma, 6),
        if (x$gamma_estimated) "(estimated)" else "(fixed)", "\n")
  }
  if (!is.null(x$theta)) {
    cat("Theta (NB dispersion):", round(x$theta, 4), "\n")
  }
  if (!is.null(x$loglik) && length(x$loglik) == 1 && is.finite(x$loglik)) {
    cat("Log-likelihood:", round(x$loglik, 2), "\n")
  }
  if (isTRUE(x$constrained)) {
    cat("Constrained: alpha_values = logit^-1(X_alpha %*% alpha_coefs), ",
        "beta_values = exp(X_beta %*% beta_coefs)\n", sep = "")
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
#' (normal approximation) are used. Standard errors come from the covariance
#' requested by the \code{vcov} argument of the original
#' \code{\link{estimate_hidden_pop}} call: an HC type (\code{"HC3"} by
#' default, \code{"HC1"} for \code{estimator = "gmm"} and
#' \code{estimator = "el"}), or whatever a user-supplied \code{vcov} function
#' returns. The type in force is shown in the printed header.
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

  use_t <- object$method %in% c("ols", "nls")
  if (use_t) {
    df <- n - p
    pval <- 2 * pt(-abs(z), df = df)
  } else {
    pval <- 2 * pnorm(-abs(z))
  }

  stat_name <- if (use_t) "t value" else "z value"
  p_name <- if (use_t) "Pr(>|t|)" else "Pr(>|z|)"
  tab <- cbind(Estimate = coefs, `Std. Error` = se, z, pval)
  colnames(tab)[3:4] <- c(stat_name, p_name)

  cat("Unauthorized population estimation\n")
  vcov_label <- if (is.function(object$vcov_spec)) {
    "user-supplied function"
  } else if (!is.null(object$cluster_var)) {
    paste0("CL(", object$vcov_type, "), cluster: ", deparse(object$call$cluster))
  } else {
    object$vcov_type
  }
  cat("Method:", toupper(object$method), "| estimator:", toupper(object$estimator),
      "| link_rho:", object$link_rho, "| vcov:", vcov_label, "\n")
  cat("N obs:", object$n_obs, "\n")
  if (!is.null(object$gamma_values) && length(object$gamma_values) > 0 &&
      !is.null(object$p_gamma) && isTRUE(object$has_cov_gamma)) {
    cat("Gamma: covariate-varying (response scale), range [",
        round(min(object$gamma_values), 6), ", ",
        round(max(object$gamma_values), 6), "]\n", sep = "")
  } else if (!is.null(object$gamma)) {
    cat("Gamma:", round(object$gamma, 6),
        if (object$gamma_estimated) "(estimated)" else "(fixed)", "\n")
  }
  if (!is.null(object$theta)) {
    cat("Theta (NB dispersion):", round(object$theta, 4), "\n")
  }
  if (!is.null(object$loglik) && length(object$loglik) == 1 &&
      is.finite(object$loglik)) {
    cat("Log-likelihood:", round(object$loglik, 2), "\n")
  }
  if (!is.null(object$loglik) && length(object$loglik) == 1 &&
      is.finite(object$loglik)) {
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
    cat("\nResponse-scale fitted parameters ",
        "(alpha_values in (0,1), beta_values > 0):\n", sep = "")
    .print_response_summary(object)
  }

  # Response-scale gamma summary for cov_gamma models
  if (!is.null(object$p_gamma) && isTRUE(object$has_cov_gamma) &&
      !is.null(object$gamma_coefs) && !is.null(object$X_gamma)) {
    cat("\nGamma (response scale):\n")
    .print_gamma_response(object)
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
  # If vcov was a function and update doesn't explicitly override it,
  # replace with a simple built-in type to prevent infinite recursion from
  # vcovFWB etc.
  if (is.function(object$vcov_spec) && !"vcov" %in% names(extras)) {
    call[["vcov"]] <- if (identical(object$estimator %in% c("gmm", "el"), TRUE)) {
      "HC1"
    } else {
      "HC3"
    }
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

#' Print response-scale gamma values per unique covariate group
#' @noRd
.print_gamma_response <- function(x) {
  X_gamma <- x$X_gamma
  gamma_coefs <- x$gamma_coefs
  p_gamma <- x$p_gamma

  # Get the gamma block of vcov_full (gamma coefs are on log scale)
  # For Poisson: gamma cols are at positions p_ab+1 : p_ab+p_gamma in vcov_full
  # For NB: gamma cols are at p_ab+1+1 : p_ab+1+p_gamma (theta at p_ab+1)
  V_gamma <- NULL
  if (!is.null(x$vcov_full)) {
    nms <- colnames(x$vcov_full)
    gamma_idx <- grep("^gamma", nms)
    if (length(gamma_idx) == p_gamma) {
      V_gamma <- x$vcov_full[gamma_idx, gamma_idx, drop = FALSE]
    }
  }

  # Unique gamma groups
  X_gamma_key <- apply(X_gamma, 1, paste, collapse = "|")
  unique_keys <- unique(X_gamma_key)

  gamma_rows <- lapply(unique_keys, function(key) {
    idx <- which(X_gamma_key == key)[1]
    z_row <- X_gamma[idx, , drop = FALSE]
    log_gamma <- as.numeric(z_row %*% gamma_coefs)
    gamma_resp <- exp(log_gamma)
    # Delta method: SE(gamma) = gamma * SE(log_gamma)
    se_log_gamma <- if (!is.null(V_gamma)) {
      sqrt(max(0, as.numeric(z_row %*% V_gamma %*% t(z_row))))
    } else {
      NA_real_
    }
    se_gamma <- gamma_resp * se_log_gamma
    data.frame(
      gamma = gamma_resp,
      `SE(gamma)` = se_gamma,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
  gamma_tab <- do.call(rbind, gamma_rows)
  rownames(gamma_tab) <- unique_keys
  print(round(gamma_tab, 6))
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
            !all(is.na(ps$estimate_bc))
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
