# ---- Shared state for count-model moment estimators ----

#' @noRd
.count_param_state <- function(par, log_N, ratio, X_alpha, X_beta,
                               estimate_gamma = FALSE, gamma_value = NULL,
                               X_gamma = NULL, constrained = FALSE,
                               link_rho = "power", has_theta = FALSE) {
  n_obs <- length(log_N)
  p_alpha <- ncol(X_alpha)
  p_beta <- ncol(X_beta)
  p_gamma <- if (!is.null(X_gamma)) ncol(X_gamma) else 0L

  idx <- 0L
  a_raw <- par[idx + seq_len(p_alpha)]
  idx <- idx + p_alpha
  b_raw <- par[idx + seq_len(p_beta)]
  idx <- idx + p_beta

  theta <- NULL
  if (has_theta) {
    theta <- exp(par[idx + 1L])
    idx <- idx + 1L
  }

  if (estimate_gamma) {
    gamma_coefs <- par[idx + seq_len(p_gamma)]
    gamma_values <- exp(pmin(as.numeric(X_gamma %*% gamma_coefs), 10))
  } else {
    gamma_coefs <- NULL
    gamma_values <- if (is.null(gamma_value)) NULL else rep(gamma_value, n_obs)
  }

  eta_alpha <- as.numeric(X_alpha %*% a_raw)
  zeta_beta <- as.numeric(X_beta %*% b_raw)
  alpha_values <- if (constrained) .inv_logit(eta_alpha) else eta_alpha
  beta_values <- if (constrained) exp(zeta_beta) else zeta_beta
  dalpha <- if (constrained) alpha_values * (1 - alpha_values) else rep(1, n_obs)
  dbeta <- if (constrained) beta_values else rep(1, n_obs)

  rate_values <- .rate_from_gamma(ratio, gamma_values)
  eta_values <- .eta_from_rate(beta_values, rate_values)
  dlogrho_values <- .dlog_rho_deta(eta_values, link_rho = link_rho)
  log_mu <- pmin(.compute_log_mu(alpha_values, log_N, beta_values, rate_values,
                                 link_rho = link_rho), 20)
  mu <- pmax(exp(log_mu), 1e-300)

  list(
    alpha_coefs = a_raw,
    beta_coefs = b_raw,
    gamma_coefs = gamma_coefs,
    gamma_values = gamma_values,
    theta = theta,
    alpha_values = alpha_values,
    beta_values = beta_values,
    dalpha = dalpha,
    dbeta = dbeta,
    rate_values = rate_values,
    eta_values = eta_values,
    dlogrho_values = dlogrho_values,
    log_mu = log_mu,
    mu = mu
  )
}

#' @noRd
.count_coef_jacobian <- function(state, log_N, X_alpha, X_beta, X_gamma = NULL) {
  alpha_cols <- (log_N * state$dalpha) * X_alpha
  beta_cols <- (state$dlogrho_values * log(state$rate_values) * state$dbeta) * X_beta
  gamma_cols <- NULL
  if (!is.null(X_gamma) && !is.null(state$gamma_values)) {
    gamma_cols <- (state$dlogrho_values * state$beta_values *
                     state$gamma_values / state$rate_values) * X_gamma
    colnames(gamma_cols) <- colnames(X_gamma)
  }
  list(alpha = alpha_cols, beta = beta_cols, gamma = gamma_cols)
}

#' @noRd
.make_count_moment_data <- function(m, log_N, ratio, weights) {
  cbind(m = m, log_N = log_N, ratio = ratio, obs_weight = weights)
}

#' @noRd
.make_poisson_moment_function <- function(X_alpha, X_beta, estimate_gamma = FALSE,
                                          gamma_value = NULL, X_gamma = NULL,
                                          constrained = FALSE,
                                          link_rho = "power") {
  force(X_alpha)
  force(X_beta)
  force(estimate_gamma)
  force(gamma_value)
  force(X_gamma)
  force(constrained)
  force(link_rho)

  function(par, x) {
    state <- .count_param_state(
      par = par,
      log_N = x[, "log_N"],
      ratio = x[, "ratio"],
      X_alpha = X_alpha,
      X_beta = X_beta,
      estimate_gamma = estimate_gamma,
      gamma_value = gamma_value,
      X_gamma = X_gamma,
      constrained = constrained,
      link_rho = link_rho,
      has_theta = FALSE
    )
    jac <- .count_coef_jacobian(state, x[, "log_N"], X_alpha, X_beta, X_gamma)
    scalar_score <- x[, "obs_weight"] * (x[, "m"] - state$mu)
    out <- cbind(jac$alpha, jac$beta)
    if (!is.null(jac$gamma)) out <- cbind(out, jac$gamma)
    out * scalar_score
  }
}

#' @noRd
.make_nb_moment_function <- function(X_alpha, X_beta, estimate_gamma = FALSE,
                                     gamma_value = NULL, X_gamma = NULL,
                                     constrained = FALSE,
                                     link_rho = "power") {
  force(X_alpha)
  force(X_beta)
  force(estimate_gamma)
  force(gamma_value)
  force(X_gamma)
  force(constrained)
  force(link_rho)

  function(par, x) {
    state <- .count_param_state(
      par = par,
      log_N = x[, "log_N"],
      ratio = x[, "ratio"],
      X_alpha = X_alpha,
      X_beta = X_beta,
      estimate_gamma = estimate_gamma,
      gamma_value = gamma_value,
      X_gamma = X_gamma,
      constrained = constrained,
      link_rho = link_rho,
      has_theta = TRUE
    )
    jac <- .count_coef_jacobian(state, x[, "log_N"], X_alpha, X_beta, X_gamma)
    scalar_score <- x[, "obs_weight"] * (x[, "m"] - state$mu) *
      state$theta / (state$theta + state$mu)
    score_theta <- x[, "obs_weight"] * (
      digamma(x[, "m"] + state$theta) - digamma(state$theta) +
        log(state$theta) + 1 - log(state$theta + state$mu) -
        (x[, "m"] + state$theta) / (state$theta + state$mu)
    ) * state$theta

    out <- cbind(jac$alpha, jac$beta)
    out <- cbind(out * scalar_score, theta = score_theta)
    if (!is.null(jac$gamma)) {
      out <- cbind(out, jac$gamma * scalar_score)
    }
    out
  }
}

#' @noRd
.moment_vcov_par <- function(fit, estimator) {
  vcov_method <- methods::selectMethod("vcov", class(fit))
  if (identical(estimator, "gmm")) {
    return(vcov_method(fit))
  }
  vcov_method(fit)$vcov_par
}

#' @noRd
.moment_model_vcov_par <- function(fit, estimator) {
  vcov_method <- methods::selectMethod("vcov", class(fit))
  if (identical(estimator, "gmm")) {
    return(tryCatch(vcov_method(fit, sandwich = FALSE),
                    error = function(e) vcov_method(fit)))
  }
  .moment_vcov_par(fit, estimator)
}

#' @noRd
.moment_convergence_code <- function(fit, estimator) {
  if (identical(estimator, "gmm")) {
    conv <- fit@convergence
    if (is.list(conv) && !is.null(conv$code)) return(as.integer(conv$code))
    return(as.integer(conv))
  }
  if (isTRUE(fit@convergence == 0) && isTRUE(fit@lconvergence == 0)) 0L else 1L
}

#' @noRd
.fit_poisson_moment <- function(m, N, ratio, log_N, log_rate,
                                X_alpha, X_beta, gamma_value = NULL,
                                estimate_gamma = FALSE, gamma_start = 0.01,
                                gamma_bounds = c(1e-10, 1),
                                X_gamma = NULL,
                                weights = NULL, constrained = FALSE,
                                link_rho = "power",
                                estimator = c("gmm", "el")) {
  estimator <- match.arg(estimator)
  n_obs <- length(m)
  if (is.null(weights)) weights <- rep(1, n_obs)
  p_alpha <- ncol(X_alpha)
  p_beta <- ncol(X_beta)
  p_gamma <- if (!is.null(X_gamma)) ncol(X_gamma) else 0L

  start_fit <- .fit_poisson(
    m = m, N = N, ratio = ratio, log_N = log_N, log_rate = log_rate,
    X_alpha = X_alpha, X_beta = X_beta,
    gamma_value = gamma_value,
    estimate_gamma = estimate_gamma,
    gamma_start = gamma_start,
    gamma_bounds = gamma_bounds,
    X_gamma = X_gamma,
    weights = weights,
    vcov_type = "HC0",
    constrained = constrained,
    link_rho = link_rho
  )

  theta0 <- c(start_fit$alpha_coefs, start_fit$beta_coefs)
  if (estimate_gamma) theta0 <- c(theta0, start_fit$gamma_coefs)
  names(theta0) <- c(colnames(X_alpha), colnames(X_beta),
                     if (estimate_gamma) colnames(X_gamma))

  moment_data <- .make_count_moment_data(m, log_N, ratio, weights)
  moment_fn <- .make_poisson_moment_function(
    X_alpha = X_alpha,
    X_beta = X_beta,
    estimate_gamma = estimate_gamma,
    gamma_value = gamma_value,
    X_gamma = X_gamma,
    constrained = constrained,
    link_rho = link_rho
  )
  model <- momentfit::momentModel(moment_fn, x = moment_data, theta0 = theta0,
                                  vcov = "MDS")
  fit <- if (identical(estimator, "gmm")) {
    momentfit::gmmFit(model, type = "twostep", weights = "ident")
  } else {
    momentfit::gelFit(model, gelType = "EL", vcov = TRUE)
  }

  par_hat <- as.numeric(fit@theta)
  names(par_hat) <- names(theta0)
  state <- .count_param_state(
    par = par_hat,
    log_N = log_N,
    ratio = ratio,
    X_alpha = X_alpha,
    X_beta = X_beta,
    estimate_gamma = estimate_gamma,
    gamma_value = gamma_value,
    X_gamma = X_gamma,
    constrained = constrained,
    link_rho = link_rho,
    has_theta = FALSE
  )
  jac <- .count_coef_jacobian(state, log_N, X_alpha, X_beta, X_gamma)
  Z <- cbind(jac$alpha, jac$beta)
  if (!is.null(jac$gamma)) Z <- cbind(Z, jac$gamma)
  colnames(Z) <- c(colnames(X_alpha), colnames(X_beta),
                   if (!is.null(jac$gamma)) colnames(X_gamma))

  vcov_full_backend <- .moment_vcov_par(fit, estimator)
  vcov_model_backend <- .moment_model_vcov_par(fit, estimator)
  coef_idx <- seq_len(ncol(Z))

  list(
    alpha_coefs = state$alpha_coefs,
    beta_coefs = state$beta_coefs,
    alpha_values = state$alpha_values,
    beta_values = state$beta_values,
    gamma_coefs = state$gamma_coefs,
    gamma_values = state$gamma_values,
    rho_values = .rho_from_eta(state$eta_values, link_rho = link_rho),
    fitted = state$mu,
    residuals = m - state$mu,
    log_mu = state$log_mu,
    vcov = vcov_full_backend[coef_idx, coef_idx, drop = FALSE],
    vcov_model = vcov_model_backend[coef_idx, coef_idx, drop = FALSE],
    vcov_full = vcov_full_backend[coef_idx, coef_idx, drop = FALSE],
    backend_vcov_full = vcov_full_backend,
    backend_vcov_model_full = vcov_model_backend,
    model_matrix_full = Z,
    bread_weights = weights * state$mu,
    score_residuals = weights * (m - state$mu),
    n_obs = n_obs,
    df.residual = n_obs - length(theta0),
    loglik = NA_real_,
    theta = NULL,
    gamma_estimated = if (estimate_gamma) {
      if (p_gamma == 1) state$gamma_values[1] else NULL
    } else NULL,
    convergence = .moment_convergence_code(fit, estimator),
    moment_backend = fit
  )
}

#' @noRd
.fit_nb_moment <- function(m, N, ratio, log_N, log_rate,
                           X_alpha, X_beta, gamma_value = NULL,
                           estimate_gamma = FALSE, gamma_start = 0.01,
                           gamma_bounds = c(1e-10, 1), theta_start = 1,
                           X_gamma = NULL,
                           weights = NULL, constrained = FALSE,
                           link_rho = "power",
                           estimator = c("gmm", "el")) {
  estimator <- match.arg(estimator)
  n_obs <- length(m)
  if (is.null(weights)) weights <- rep(1, n_obs)
  p_alpha <- ncol(X_alpha)
  p_beta <- ncol(X_beta)
  p_gamma <- if (!is.null(X_gamma)) ncol(X_gamma) else 0L

  start_fit <- .fit_nb(
    m = m, N = N, ratio = ratio, log_N = log_N, log_rate = log_rate,
    X_alpha = X_alpha, X_beta = X_beta,
    gamma_value = gamma_value,
    estimate_gamma = estimate_gamma,
    gamma_start = gamma_start,
    gamma_bounds = gamma_bounds,
    theta_start = theta_start,
    X_gamma = X_gamma,
    weights = weights,
    vcov_type = "HC0",
    constrained = constrained,
    link_rho = link_rho
  )

  theta0 <- c(start_fit$alpha_coefs, start_fit$beta_coefs, theta = log(start_fit$theta))
  if (estimate_gamma) theta0 <- c(theta0, start_fit$gamma_coefs)
  names(theta0) <- c(colnames(X_alpha), colnames(X_beta), "theta",
                     if (estimate_gamma) colnames(X_gamma))

  moment_data <- .make_count_moment_data(m, log_N, ratio, weights)
  moment_fn <- .make_nb_moment_function(
    X_alpha = X_alpha,
    X_beta = X_beta,
    estimate_gamma = estimate_gamma,
    gamma_value = gamma_value,
    X_gamma = X_gamma,
    constrained = constrained,
    link_rho = link_rho
  )
  model <- momentfit::momentModel(moment_fn, x = moment_data, theta0 = theta0,
                                  vcov = "MDS")
  fit <- if (identical(estimator, "gmm")) {
    momentfit::gmmFit(model, type = "twostep", weights = "ident")
  } else {
    momentfit::gelFit(model, gelType = "EL", vcov = TRUE)
  }

  par_hat <- as.numeric(fit@theta)
  names(par_hat) <- names(theta0)
  state <- .count_param_state(
    par = par_hat,
    log_N = log_N,
    ratio = ratio,
    X_alpha = X_alpha,
    X_beta = X_beta,
    estimate_gamma = estimate_gamma,
    gamma_value = gamma_value,
    X_gamma = X_gamma,
    constrained = constrained,
    link_rho = link_rho,
    has_theta = TRUE
  )
  jac <- .count_coef_jacobian(state, log_N, X_alpha, X_beta, X_gamma)
  Z <- cbind(jac$alpha, jac$beta)
  if (!is.null(jac$gamma)) Z <- cbind(Z, jac$gamma)
  colnames(Z) <- c(colnames(X_alpha), colnames(X_beta),
                   if (!is.null(jac$gamma)) colnames(X_gamma))

  nb_score_factor <- weights * (m - state$mu) * state$theta / (state$theta + state$mu)
  score_theta <- weights * (
    digamma(m + state$theta) - digamma(state$theta) +
      log(state$theta) + 1 - log(state$theta + state$mu) -
      (m + state$theta) / (state$theta + state$mu)
  ) * state$theta
  score_full <- cbind(Z * nb_score_factor, theta = score_theta)
  if (!is.null(jac$gamma)) {
    gamma_score <- jac$gamma * nb_score_factor
    score_full <- cbind(Z[, seq_len(p_alpha + p_beta), drop = FALSE] * nb_score_factor,
                        theta = score_theta, gamma_score)
    colnames(score_full) <- c(colnames(X_alpha), colnames(X_beta), "theta",
                              colnames(X_gamma))
  } else {
    colnames(score_full) <- c(colnames(X_alpha), colnames(X_beta), "theta")
  }

  vcov_full_backend <- .moment_vcov_par(fit, estimator)
  vcov_model_backend <- .moment_model_vcov_par(fit, estimator)
  p_ab <- p_alpha + p_beta
  coef_idx <- c(seq_len(p_ab), if (estimate_gamma) p_ab + 1L + seq_len(p_gamma))

  theta_idx <- p_ab + 1L
  theta_se <- state$theta * sqrt(max(0, vcov_full_backend[theta_idx, theta_idx]))

  list(
    alpha_coefs = state$alpha_coefs,
    beta_coefs = state$beta_coefs,
    alpha_values = state$alpha_values,
    beta_values = state$beta_values,
    gamma_coefs = state$gamma_coefs,
    gamma_values = state$gamma_values,
    rho_values = .rho_from_eta(state$eta_values, link_rho = link_rho),
    fitted = state$mu,
    residuals = m - state$mu,
    log_mu = state$log_mu,
    vcov = vcov_full_backend[coef_idx, coef_idx, drop = FALSE],
    vcov_model = vcov_model_backend[coef_idx, coef_idx, drop = FALSE],
    vcov_full = vcov_full_backend[coef_idx, coef_idx, drop = FALSE],
    backend_vcov_full = vcov_full_backend,
    backend_vcov_model_full = vcov_model_backend,
    model_matrix_full = Z,
    bread_weights = weights * state$mu * state$theta / (state$theta + state$mu),
    score_residuals = nb_score_factor,
    score_full = score_full,
    theta_se = theta_se,
    n_obs = n_obs,
    df.residual = n_obs - length(theta0),
    loglik = NA_real_,
    theta = state$theta,
    gamma_estimated = if (estimate_gamma) {
      if (p_gamma == 1) state$gamma_values[1] else NULL
    } else NULL,
    convergence = .moment_convergence_code(fit, estimator),
    moment_backend = fit
  )
}
