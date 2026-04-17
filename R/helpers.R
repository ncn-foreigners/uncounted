# ---- Constrained parameter transforms ----

#' Inverse logit (logistic function): R -> (0, 1)
#' @noRd
.inv_logit <- function(x) stats::plogis(x)

#' Derivative of inverse logit
#' @noRd
.inv_logit_deriv <- function(x) {
  s <- .inv_logit(x)
  s * (1 - s)
}

#' Logit: (0, 1) -> R
#' @noRd
.logit <- function(p) stats::qlogis(p)

# ---- Detection-rate link helpers ----

#' Stable log(1 - exp(-x)) for x > 0
#' @noRd
.log1mexp_neg <- function(x) {
  ifelse(x <= log(2), log(-expm1(-x)), log1p(-exp(-x)))
}

#' Normalize detection-link name
#' @noRd
.normalize_link_rho <- function(link_rho) {
  if (is.null(link_rho)) return("power")
  if (length(link_rho) > 1) link_rho <- link_rho[[1]]
  if (identical(link_rho, "logistic")) link_rho <- "logit"
  match.arg(link_rho, c("power", "cloglog", "logit"))
}

#' Compute rate term gamma + n/N
#' @noRd
.rate_from_gamma <- function(ratio, gamma_values = NULL) {
  if (is.null(gamma_values)) ratio else gamma_values + ratio
}

#' Compute log-rate term log(gamma + n/N)
#' @noRd
.log_rate_from_gamma <- function(ratio, gamma_values = NULL) {
  log(.rate_from_gamma(ratio, gamma_values))
}

#' Compute eta = beta * log(gamma + n/N)
#' @noRd
.eta_from_rate <- function(beta_values, rate_values) {
  beta_values * log(rate_values)
}

#' Compute rho from eta
#' @noRd
.rho_from_eta <- function(eta, link_rho = "power") {
  link_rho <- .normalize_link_rho(link_rho)
  switch(link_rho,
    power = exp(eta),
    cloglog = -expm1(-exp(eta)),
    logit = .inv_logit(eta)
  )
}

#' Compute log(rho) from eta
#' @noRd
.log_rho_from_eta <- function(eta, link_rho = "power") {
  link_rho <- .normalize_link_rho(link_rho)
  switch(link_rho,
    power = eta,
    cloglog = .log1mexp_neg(exp(eta)),
    logit = stats::plogis(eta, log.p = TRUE)
  )
}

#' Compute d log(rho) / d eta
#' @noRd
.dlog_rho_deta <- function(eta, link_rho = "power") {
  link_rho <- .normalize_link_rho(link_rho)
  switch(link_rho,
    power = rep(1, length(eta)),
    cloglog = {
      exp_eta <- exp(eta)
      out <- exp_eta / expm1(exp_eta)
      out[!is.finite(out) & eta > 0] <- 0
      out
    },
    logit = .inv_logit(-eta)
  )
}

#' Compute log(mu) from alpha, beta, rate, and link choice
#' @noRd
.compute_log_mu <- function(alpha_values, log_N, beta_values, rate_values,
                            link_rho = "power") {
  eta_values <- .eta_from_rate(beta_values, rate_values)
  alpha_values * log_N + .log_rho_from_eta(eta_values, link_rho = link_rho)
}

#' Normalize requested covariance type for moment estimators
#' @noRd
.normalize_moment_vcov_type <- function(vcov_type) {
  if (vcov_type %in% c("HC0", "HC1")) return(vcov_type)
  if (vcov_type %in% c("HC2", "HC3", "HC4", "HC4m", "HC5")) {
    warning("For estimator = 'gmm' and estimator = 'el', HC2+ types are not ",
            "supported. Falling back to HC1.", call. = FALSE)
    return("HC1")
  }
  vcov_type
}

#' Normalize requested covariance type for a fitted object
#' @noRd
.normalize_object_vcov_type <- function(object, vcov_type) {
  if (inherits(object, "uncounted") && identical(object$estimator %in% c("gmm", "el"), TRUE)) {
    return(.normalize_moment_vcov_type(vcov_type))
  }
  vcov_type
}

# ---- Design matrix helpers ----

#' Build design matrix from a covariate formula
#'
#' @param cov_formula One-sided formula (e.g., ~ sex + year) or NULL
#' @param data Data frame
#' @return Model matrix (n x p). If NULL, returns intercept-only (column of 1s).
#' @noRd
.build_model_matrix <- function(cov_formula, data) {
  if (is.null(cov_formula)) {
    return(matrix(1, nrow = nrow(data), ncol = 1,
                  dimnames = list(NULL, "(Intercept)")))
  }
  X <- model.matrix(cov_formula, data = data)
  if (ncol(X) == 0) {
    # ~0 or ~-1 produces empty matrix; treat as intercept-only
    return(matrix(1, nrow = nrow(data), ncol = 1,
                  dimnames = list(NULL, "(Intercept)")))
  }
  X
}

#' Safe matrix solve with fallback to generalized inverse
#' @noRd
.solve_safe <- function(A, b = NULL) {
  tryCatch({
    if (is.null(b)) solve(A) else solve(A, b)
  }, error = function(e) {
    tryCatch({
      if (is.null(b)) qr.solve(A) else qr.solve(A, b)
    }, error = function(e2) {
      Ainv <- MASS::ginv(A)
      if (is.null(b)) Ainv else Ainv %*% b
    })
  })
}

#' Compute sandwich HC variance-covariance matrix
#'
#' @param X Design matrix (n x p)
#' @param residuals Residuals vector (n)
#' @param weights Working weights for bread (n). For OLS: rep(1,n). For Poisson: mu.
#' @param hat_values Hat matrix diagonal (n). Needed for HC2+.
#' @param vcov_type Character: "HC0", "HC1", "HC2", "HC3", "HC4", "HC4m", "HC5"
#' @param cluster Optional integer/factor vector (n) for cluster-robust variance.
#'   When provided, score contributions are summed within clusters (Liang-Zeger 1986)
#'   and a small-sample correction G/(G-1) * (n-1)/(n-p) is applied.
#' @return p x p variance-covariance matrix
#' @noRd
.compute_sandwich_vcov <- function(X, residuals, weights = NULL,
                                   hat_values = NULL, vcov_type = "HC3",
                                   cluster = NULL) {
  n <- nrow(X)
  p <- ncol(X)

  if (is.null(weights)) weights <- rep(1, n)

  # Bread: (X' W X)^{-1}
  wX <- X * sqrt(weights)
  bread <- .solve_safe(crossprod(wX))

  if (!is.null(cluster)) {
    ## ── Cluster-robust sandwich ──
    ## Map HC type to cluster-robust (CR) type:
    ##   HC0 -> CR0 (no small-sample correction)
    ##   HC1 -> CR1 (G/(G-1) * (n-1)/(n-p) correction; Liang-Zeger 1986)
    ##   HC2+ not supported with clustering; fall back to CR1 with warning

    cr_type <- switch(vcov_type,
      "HC0" = "CR0",
      "HC1" = "CR1",
      {
        warning(sprintf(
          paste0("vcov_type = '%s' is not supported with cluster-robust variance. ",
                 "Falling back to CR1 (equivalent to HC1 with Liang-Zeger ",
                 "small-sample correction). Use 'HC0' for CR0 (no correction) ",
                 "or 'HC1' for CR1."),
          vcov_type
        ), call. = FALSE)
        "CR1"
      }
    )

    ## Score per observation: w_i * e_i * x_i = wX_i * (e_i * sqrt(w_i))
    score_i <- wX * as.numeric(residuals * sqrt(weights))

    ## Sum scores within clusters
    cl <- as.factor(cluster)
    G <- nlevels(cl)
    if (G < 2) {
      stop("Cluster-robust variance requires at least 2 clusters.",
           call. = FALSE)
    }
    score_g <- rowsum(score_i, cl, reorder = FALSE)

    ## Meat = sum_g (S_g' S_g)
    meat <- crossprod(score_g)

    ## Small-sample correction
    correction <- if (cr_type == "CR0") {
      1
    } else {
      ## CR1: G/(G-1) * (n-1)/(n-p)  (Liang-Zeger 1986)
      (G / (G - 1)) * ((n - 1) / (n - p))
    }
    return(bread %*% (correction * meat) %*% bread)
  }

  ## ── Observation-level HC variance ──
  # Hat values (needed for HC2+)
  h <- if (!is.null(hat_values)) {
    hat_values
  } else if (vcov_type %in% c("HC2", "HC3", "HC4", "HC4m", "HC5")) {
    .hat_values_wls(X, weights)
  } else {
    NULL
  }

  # Meat adjustment factor per observation
  e2 <- residuals^2
  adj <- switch(vcov_type,
    "HC0" = rep(1, n),
    "HC1" = rep(n / (n - p), n),
    "HC2" = 1 / (1 - h),
    "HC3" = 1 / (1 - h)^2,
    "HC4" = {
      delta <- pmin(4, n * h / p)
      1 / (1 - h)^delta
    },
    "HC4m" = {
      delta <- pmin(1, n * h / p) + pmin(1.5, n * h / p)
      1 / (1 - h)^delta
    },
    "HC5" = {
      k <- pmin(n * h / p, pmax(4, 0.7 * n * h / p))
      1 / sqrt((1 - h)^k)
    },
    rep(1, n)  # fallback to HC0
  )

  # Meat: (sqrt(W)*X)' diag(adj * e^2) (sqrt(W)*X)
  meat <- crossprod(wX * sqrt(adj * e2))

  bread %*% meat %*% bread
}

#' Compute model-based (non-sandwich) variance-covariance matrix
#'
#' For bias correction of xi, the paper uses the standard (homoscedastic)
#' variance estimate rather than HC-robust, because HC3 can be inflated
#' by high-leverage observations, leading to overcorrection.
#'
#' @param X Design matrix (n x p)
#' @param weights Working weights (n). OLS: rep(1,n). Poisson: mu. NB: mu*theta/(theta+mu).
#' @param sigma2 Residual variance for OLS/NLS (scalar). NULL for Poisson/NB.
#' @return p x p model-based variance-covariance matrix
#' @noRd
.compute_model_vcov <- function(X, weights, sigma2 = NULL) {
  wX <- X * sqrt(weights)
  V <- .solve_safe(crossprod(wX))
  if (!is.null(sigma2)) V <- sigma2 * V  # OLS/NLS: sigma^2 * (X'X)^{-1}
  V
}

#' NB-specific sandwich vcov using the full score and Hessian
#'
#' For the NB model, theta enters the likelihood in a way that cannot be
#' factored through the mean-model Jacobian (unlike alpha/beta/gamma).
#' This function computes the sandwich vcov using the full per-observation
#' score vector and the optim Hessian, bypassing the generic sandwich
#' package machinery (which assumes \code{estfun / model.matrix} factoring).
#'
#' @param hessian_nll Hessian of the negative log-likelihood from \code{optim()},
#'   on the optimizer scale (log-theta, log-gamma). Dimension k x k where
#'   k = p_alpha + p_beta + 1 for theta, optionally + 1 for gamma.
#' @param score_full Per-observation score matrix (n x k), on the same scale.
#' @param n_obs Number of observations.
#' @param vcov_type HC type: \code{"HC0"} or \code{"HC1"}. HC2+ not supported
#'   for NB with theta; falls back to HC1 with a message.
#' @param cluster Optional factor vector for cluster-robust variance.
#' @return k x k sandwich variance-covariance matrix.
#' @noRd
.compute_nb_sandwich <- function(hessian_nll, score_full, n_obs,
                                  vcov_type = "HC0", cluster = NULL) {
  H_inv <- .solve_safe(hessian_nll)
  p <- ncol(score_full)

  if (!is.null(cluster)) {
    cl <- as.factor(cluster)
    G <- nlevels(cl)
    if (G < 2) {
      stop("Cluster-robust variance requires at least 2 clusters.",
           call. = FALSE)
    }
    score_g <- rowsum(score_full, cl, reorder = FALSE)
    meat <- crossprod(score_g)
    # HC0 (CR0): no correction; HC1 (CR1): small-sample correction
    correction <- if (vcov_type == "HC0") {
      1
    } else {
      (G / (G - 1)) * ((n_obs - 1) / (n_obs - p))
    }
    return(H_inv %*% (correction * meat) %*% H_inv)
  }

  meat <- crossprod(score_full)
  if (vcov_type %in% c("HC2", "HC3", "HC4", "HC4m", "HC5")) {
    vcov_type <- "HC1"
  }
  if (vcov_type == "HC1") {
    meat <- meat * n_obs / (n_obs - p)
  }
  H_inv %*% meat %*% H_inv
}

#' Compute hat values for weighted least squares
#' @noRd
.hat_values_wls <- function(X, weights) {
  wX <- X * sqrt(weights)
  Q <- qr(wX)
  Q_mat <- qr.Q(Q)
  rowSums(Q_mat^2)
}

#' Compute log-rate term: log(gamma + n/N) handling gamma configurations
#'
#' @param ratio n/N vector
#' @param gamma_values Gamma value(s). Scalar or vector matching observations.
#' @return log(gamma + ratio) vector
#' @noRd
.log_rate <- function(ratio, gamma_values = NULL) {
  .log_rate_from_gamma(ratio, gamma_values = gamma_values)
}

#' Expand gamma from group-level to observation-level
#'
#' @param gamma_params Named vector of gamma values per group
#' @param gamma_groups Factor or character vector assigning each obs to a group
#' @return Numeric vector of gamma values (one per observation)
#' @noRd
.expand_gamma <- function(gamma_params, gamma_groups) {
  gamma_params[as.character(gamma_groups)]
}

#' Compute vcov using the sandwich package
#'
#' Dispatches to \code{sandwich::vcovHC()}, \code{sandwich::vcovCL()}, or a
#' user-supplied function, depending on the \code{vcov} argument.
#'
#' @param object An \code{uncounted} object (must have bread/estfun methods).
#' @param vcov A character string (HC type) or a function.
#' @param cluster Optional cluster vector (already extracted from formula).
#' @return p x p variance-covariance matrix (full, including gamma if estimated).
#' @noRd
.compute_vcov_sandwich <- function(object, vcov, cluster = NULL) {
  if (is.function(vcov)) {
    ## User-supplied function: call it on the object
    V <- vcov(object)
    return(V)
  }

  vcov <- .normalize_object_vcov_type(object, vcov)

  ## Character string: use sandwich package
  if (!requireNamespace("sandwich", quietly = TRUE)) {
    stop("Package 'sandwich' is required for vcov computation. ",
         "Install with install.packages('sandwich')")
  }

  if (!is.null(cluster)) {
    ## Cluster-robust: sandwich::vcovCL
    V <- sandwich::vcovCL(object, cluster = cluster, type = vcov)
  } else {
    ## Observation-level HC: sandwich::vcovHC
    V <- sandwich::vcovHC(object, type = vcov)
  }
  V
}
