# ---- OLS estimator ----

#' @noRd
.fit_ols <- function(m, N, ratio, log_N, log_rate,
                     X_alpha, X_beta, weights = NULL, vcov_type = "HC3") {
  n_obs <- length(m)
  if (is.null(weights)) weights <- rep(1, n_obs)

  Z <- cbind(X_alpha * log_N, X_beta * log_rate)
  p_alpha <- ncol(X_alpha)
  p_beta <- ncol(X_beta)

  has_zeros <- any(m == 0)
  y <- if (has_zeros) log(m + 1) else log(m)

  fit <- if (all(weights == 1)) lm.fit(Z, y) else lm.wfit(Z, y, weights)
  coefs <- fit$coefficients
  names(coefs) <- colnames(Z)

  alpha_coefs <- coefs[seq_len(p_alpha)]
  beta_coefs <- coefs[p_alpha + seq_len(p_beta)]

  log_mu <- as.numeric(Z %*% coefs)
  mu <- exp(log_mu)
  resid_raw <- m - mu
  resid_log <- y - log_mu

  h <- .hat_values_wls(Z, weights)
  V <- .compute_sandwich_vcov(Z, resid_log, weights = weights,
                              hat_values = h, vcov_type = vcov_type,
)

  sigma2 <- sum(weights * resid_log^2) / (n_obs - length(coefs))
  V_model <- .compute_model_vcov(Z, weights, sigma2 = sigma2)

  list(alpha_coefs = alpha_coefs, beta_coefs = beta_coefs,
       fitted = mu, residuals = resid_raw, log_mu = log_mu,
       vcov = V, vcov_model = V_model,
       model_matrix_full = Z, bread_weights = weights,
       score_residuals = weights * resid_log,
       n_obs = n_obs, df.residual = n_obs - length(coefs),
       loglik = NULL, theta = NULL, gamma_estimated = NULL,
       sigma2 = sigma2)
}


# ---- OLS with profiled gamma ----

#' @noRd
.fit_ols_gamma <- function(m, N, ratio, log_N, X_alpha, X_beta,
                           gamma_start, gamma_bounds,
                           weights = NULL, vcov_type = "HC3") {
  n_obs <- length(m)
  if (is.null(weights)) weights <- rep(1, n_obs)
  p_alpha <- ncol(X_alpha)
  p_beta <- ncol(X_beta)

  has_zeros <- any(m == 0)
  y <- if (has_zeros) log(m + 1) else log(m)

  obj_fn <- function(g) {
    lr <- log(g + ratio)
    Z <- cbind(X_alpha * log_N, X_beta * lr)
    fit <- if (all(weights == 1)) lm.fit(Z, y) else lm.wfit(Z, y, weights)
    sum(weights * fit$residuals^2)
  }

  # Multi-start grid search to avoid local minima
  gamma_grid <- c(seq(1e-6, 0.01, length.out = 10),
                  seq(0.01, 0.1, length.out = 5),
                  seq(0.1, gamma_bounds[2], length.out = 5))
  grid_rss <- sapply(gamma_grid, obj_fn)
  best_idx <- which.min(grid_rss)
  best_grid <- gamma_grid[best_idx]

  # Narrow interval around best grid point
  lower_g <- if (best_idx > 1) gamma_grid[best_idx - 1] else gamma_bounds[1]
  upper_g <- if (best_idx < length(gamma_grid)) gamma_grid[best_idx + 1] else gamma_bounds[2]
  opt <- optimize(obj_fn, interval = c(lower_g, upper_g))
  gamma_hat <- opt$minimum

  if (gamma_hat > 0.5) {
    warning("Estimated gamma is large (> 0.5), suggesting possible identification issues.")
  }

  log_rate <- log(gamma_hat + ratio)
  Z <- cbind(X_alpha * log_N, X_beta * log_rate)

  fit <- if (all(weights == 1)) lm.fit(Z, y) else lm.wfit(Z, y, weights)
  coefs <- fit$coefficients
  names(coefs) <- colnames(Z)

  alpha_coefs <- coefs[seq_len(p_alpha)]
  beta_coefs <- coefs[p_alpha + seq_len(p_beta)]

  log_mu <- as.numeric(Z %*% coefs)
  mu <- exp(log_mu)
  resid_raw <- m - mu
  resid_log <- y - log_mu

  beta_lin <- as.numeric(X_beta %*% beta_coefs)
  J_gamma <- beta_lin / (gamma_hat + ratio)
  J_full <- cbind(Z, gamma = J_gamma)

  h <- .hat_values_wls(J_full, weights)
  V_full <- .compute_sandwich_vcov(J_full, resid_log, weights = weights,
                                   hat_values = h, vcov_type = vcov_type,
     )

  p_ab <- p_alpha + p_beta
  V <- V_full[seq_len(p_ab), seq_len(p_ab), drop = FALSE]

  sigma2 <- sum(weights * resid_log^2) / (n_obs - length(coefs) - 1)
  V_model_full <- .compute_model_vcov(J_full, weights, sigma2 = sigma2)
  V_model <- V_model_full[seq_len(p_ab), seq_len(p_ab), drop = FALSE]

  list(alpha_coefs = alpha_coefs, beta_coefs = beta_coefs,
       fitted = mu, residuals = resid_raw, log_mu = log_mu,
       vcov = V, vcov_model = V_model, vcov_full = V_full,
       model_matrix_full = J_full, bread_weights = weights,
       score_residuals = weights * resid_log,
       n_obs = n_obs, df.residual = n_obs - length(coefs) - 1,
       loglik = NULL, theta = NULL, gamma_estimated = gamma_hat,
       sigma2 = sigma2)
}


# ---- NLS estimator ----

#' @noRd
.fit_nls <- function(m, N, ratio, log_N, log_rate,
                     X_alpha, X_beta, gamma_value = NULL,
                     estimate_gamma = FALSE, gamma_start = 0.01,
                     gamma_bounds = c(1e-10, 1),
                     weights = NULL, vcov_type = "HC3") {
  n_obs <- length(m)
  if (is.null(weights)) weights <- rep(1, n_obs)
  p_alpha <- ncol(X_alpha)
  p_beta <- ncol(X_beta)

  obj_fn <- function(par) {
    a <- par[seq_len(p_alpha)]
    b <- par[p_alpha + seq_len(p_beta)]
    alpha_lin <- as.numeric(X_alpha %*% a)
    if (estimate_gamma) {
      g <- exp(par[p_alpha + p_beta + 1])
      lr <- log(g + ratio)
    } else {
      lr <- log_rate
    }
    beta_lin <- as.numeric(X_beta %*% b)
    mu <- exp(alpha_lin * log_N + beta_lin * lr)
    sum(weights * (m - mu)^2)
  }

  grad_fn <- function(par) {
    a <- par[seq_len(p_alpha)]
    b <- par[p_alpha + seq_len(p_beta)]
    alpha_lin <- as.numeric(X_alpha %*% a)
    if (estimate_gamma) {
      g <- exp(par[p_alpha + p_beta + 1])
      rate <- g + ratio
      lr <- log(rate)
    } else {
      lr <- log_rate
      rate <- exp(lr)
    }
    beta_lin <- as.numeric(X_beta %*% b)
    mu <- exp(alpha_lin * log_N + beta_lin * lr)
    resid <- m - mu
    g_alpha <- -2 * colSums(weights * resid * mu * log_N * X_alpha)
    g_beta <- -2 * colSums(weights * resid * mu * lr * X_beta)
    if (estimate_gamma) {
      g_gamma <- -2 * sum(weights * resid * mu * beta_lin / rate) * g
      c(g_alpha, g_beta, g_gamma)
    } else {
      c(g_alpha, g_beta)
    }
  }

  has_zeros <- any(m == 0)
  y_start <- if (has_zeros) log(m + 1) else log(m)
  Z_start <- cbind(X_alpha * log_N, X_beta * log_rate)
  start_par <- lm.fit(Z_start, y_start)$coefficients

  if (estimate_gamma) {
    start_par <- c(start_par, log(gamma_start))
    lower <- c(rep(-Inf, p_alpha + p_beta), log(gamma_bounds[1]))
    upper <- c(rep(Inf, p_alpha + p_beta), log(gamma_bounds[2]))
  } else {
    lower <- rep(-Inf, p_alpha + p_beta)
    upper <- rep(Inf, p_alpha + p_beta)
  }

  opt <- optim(start_par, obj_fn, gr = grad_fn, method = "L-BFGS-B",
               lower = lower, upper = upper,
               control = list(maxit = 1000))

  if (opt$convergence != 0) {
    warning("NLS optimization did not converge (code ", opt$convergence,
            "): ", opt$message, call. = FALSE)
  }

  a <- opt$par[seq_len(p_alpha)]
  b <- opt$par[p_alpha + seq_len(p_beta)]
  gamma_hat <- if (estimate_gamma) exp(opt$par[p_alpha + p_beta + 1]) else gamma_value

  alpha_lin <- as.numeric(X_alpha %*% a)
  lr_final <- if (estimate_gamma) log(gamma_hat + ratio) else log_rate
  beta_lin <- as.numeric(X_beta %*% b)
  log_mu <- alpha_lin * log_N + beta_lin * lr_final
  mu <- exp(log_mu)
  resid_raw <- m - mu

  J <- cbind(mu * log_N * X_alpha, mu * lr_final * X_beta)
  if (estimate_gamma) {
    J <- cbind(J, gamma = mu * beta_lin / (gamma_hat + ratio))
  }

  h <- .hat_values_wls(J, weights)
  V_full <- .compute_sandwich_vcov(J, resid_raw, weights = weights,
                                   hat_values = h, vcov_type = vcov_type,
     )
  p_ab <- p_alpha + p_beta
  V <- V_full[seq_len(p_ab), seq_len(p_ab), drop = FALSE]

  sigma2 <- sum(weights * resid_raw^2) / (n_obs - length(opt$par))
  V_model_full <- .compute_model_vcov(J, weights, sigma2 = sigma2)
  V_model <- V_model_full[seq_len(p_ab), seq_len(p_ab), drop = FALSE]

  ## NLS Jacobian: J = dmu/dtheta * mu (already computed above)
  J_full <- if (estimate_gamma) J else J
  list(alpha_coefs = a, beta_coefs = b,
       alpha_values = alpha_lin, beta_values = beta_lin,
       fitted = mu, residuals = resid_raw, log_mu = log_mu,
       vcov = V, vcov_model = V_model,
       vcov_full = if (estimate_gamma) V_full else V,
       model_matrix_full = J_full, bread_weights = weights,
       score_residuals = weights * resid_raw,
       n_obs = n_obs, df.residual = n_obs - length(opt$par),
       loglik = NULL, theta = NULL,
       gamma_estimated = if (estimate_gamma) gamma_hat else NULL,
       sigma2 = sigma2,
       convergence = opt$convergence,
       X_alpha = X_alpha, X_beta = X_beta)
}


# ---- Poisson MLE ----

#' @noRd
.fit_poisson <- function(m, N, ratio, log_N, log_rate,
                         X_alpha, X_beta, gamma_value = NULL,
                         estimate_gamma = FALSE, gamma_start = 0.01,
                         gamma_bounds = c(1e-10, 1),
                         weights = NULL, vcov_type = "HC3",
                         constrained = FALSE) {
  n_obs <- length(m)
  if (is.null(weights)) weights <- rep(1, n_obs)
  p_alpha <- ncol(X_alpha)
  p_beta <- ncol(X_beta)

  # Transform functions for constrained optimization
  # alpha = inverse logit(eta) -> (0,1); beta = exp(zeta) -> (0, Inf)
  .alpha_fn <- if (constrained) function(x) .inv_logit(x) else identity
  .beta_fn  <- if (constrained) function(x) exp(x) else identity

  nll <- function(par) {
    a <- par[seq_len(p_alpha)]
    b <- par[p_alpha + seq_len(p_beta)]
    alpha_lin <- .alpha_fn(as.numeric(X_alpha %*% a))
    if (estimate_gamma) {
      g <- exp(par[p_alpha + p_beta + 1])
      lr <- log(g + ratio)
    } else {
      lr <- log_rate
    }
    beta_lin <- .beta_fn(as.numeric(X_beta %*% b))
    log_mu <- pmin(alpha_lin * log_N + beta_lin * lr, 20)
    mu <- pmax(exp(log_mu), 1e-300)
    val <- -sum(weights * dpois(m, mu, log = TRUE))
    if (!is.finite(val)) return(1e20)
    val
  }

  grad_nll <- function(par) {
    a <- par[seq_len(p_alpha)]
    b <- par[p_alpha + seq_len(p_beta)]
    eta_alpha <- as.numeric(X_alpha %*% a)
    zeta_beta <- as.numeric(X_beta %*% b)
    alpha_lin <- .alpha_fn(eta_alpha)
    beta_lin <- .beta_fn(zeta_beta)
    if (estimate_gamma) {
      g <- exp(par[p_alpha + p_beta + 1])
      rate <- g + ratio
      lr <- log(rate)
    } else {
      lr <- log_rate
      rate <- exp(lr)
    }
    log_mu <- pmin(alpha_lin * log_N + beta_lin * lr, 20)
    mu <- pmax(exp(log_mu), 1e-300)
    score <- weights * (mu - m)

    # Chain rule for constrained: d/d(eta) = d/d(alpha) * d(alpha)/d(eta)
    if (constrained) {
      dalpha <- alpha_lin * (1 - alpha_lin)  # inverse logit derivative
      dbeta <- beta_lin                       # exp derivative
    } else {
      dalpha <- rep(1, n_obs)
      dbeta <- rep(1, n_obs)
    }

    g_alpha <- colSums(score * log_N * dalpha * X_alpha)
    g_beta <- colSums(score * lr * dbeta * X_beta)
    if (estimate_gamma) {
      g_gamma <- sum(score * beta_lin / rate) * g
      c(g_alpha, g_beta, g_gamma)
    } else {
      c(g_alpha, g_beta)
    }
  }

  # Starting values
  has_zeros <- any(m == 0)
  y_start <- if (has_zeros) log(m + 1) else log(m)
  Z_start <- cbind(X_alpha * log_N, X_beta * log_rate)
  start_par <- lm.fit(Z_start, y_start)$coefficients
  # Transform starting values to link scale
  if (constrained) {
    start_par[seq_len(p_alpha)] <- .logit(pmin(pmax(start_par[seq_len(p_alpha)], 0.01), 0.99))
    start_par[p_alpha + seq_len(p_beta)] <- log(pmax(start_par[p_alpha + seq_len(p_beta)], 0.01))
  }

  if (estimate_gamma) {
    start_par <- c(start_par, log(gamma_start))
    lower <- c(rep(-Inf, p_alpha + p_beta), log(gamma_bounds[1]))
    upper <- c(rep(Inf, p_alpha + p_beta), log(gamma_bounds[2]))
  } else {
    lower <- rep(-Inf, p_alpha + p_beta)
    upper <- rep(Inf, p_alpha + p_beta)
  }

  opt <- optim(start_par, nll, gr = grad_nll, method = "L-BFGS-B",
               lower = lower, upper = upper,
               control = list(maxit = 1000))

  if (opt$convergence != 0) {
    warning("Poisson optimization did not converge (code ", opt$convergence,
            "): ", opt$message, call. = FALSE)
  }

  a_raw <- opt$par[seq_len(p_alpha)]
  b_raw <- opt$par[p_alpha + seq_len(p_beta)]
  gamma_hat <- if (estimate_gamma) exp(opt$par[p_alpha + p_beta + 1]) else gamma_value

  # Transform to response scale
  eta_alpha <- as.numeric(X_alpha %*% a_raw)
  zeta_beta <- as.numeric(X_beta %*% b_raw)
  alpha_lin <- .alpha_fn(eta_alpha)
  beta_lin <- .beta_fn(zeta_beta)

  if (estimate_gamma) {
    rate_final <- gamma_hat + ratio
    lr_final <- log(rate_final)
  } else {
    lr_final <- log_rate
    rate_final <- exp(lr_final)
  }
  log_mu <- alpha_lin * log_N + beta_lin * lr_final
  mu <- exp(log_mu)
  resid_raw <- m - mu

  # Design matrix on link scale (d(log_mu)/d(eta,zeta))
  if (constrained) {
    dalpha <- alpha_lin * (1 - alpha_lin)
    dbeta <- beta_lin
    Z <- cbind(log_N * dalpha * X_alpha, lr_final * dbeta * X_beta)
  } else {
    Z <- cbind(log_N * X_alpha, lr_final * X_beta)
  }
  if (estimate_gamma) {
    Z <- cbind(Z, gamma = beta_lin / rate_final)
  }

  # Sandwich vcov (on link scale)
  w_glm <- weights * mu
  h <- .hat_values_wls(Z, w_glm)
  V_full <- .compute_sandwich_vcov(Z, resid_raw / sqrt(mu), weights = w_glm,
                                   hat_values = h, vcov_type = vcov_type,
     )

  p_ab <- p_alpha + p_beta
  V <- if (estimate_gamma) {
    V_full[seq_len(p_ab), seq_len(p_ab), drop = FALSE]
  } else {
    V_full
  }

  V_model_full <- .compute_model_vcov(Z, w_glm)
  V_model <- if (estimate_gamma) {
    V_model_full[seq_len(p_ab), seq_len(p_ab), drop = FALSE]
  } else {
    V_model_full
  }

  ## Z is the full model matrix (including gamma col if estimated)
  list(alpha_coefs = a_raw, beta_coefs = b_raw,
       alpha_values = alpha_lin, beta_values = beta_lin,
       fitted = mu, residuals = resid_raw, log_mu = log_mu,
       vcov = V, vcov_model = V_model,
       vcov_full = if (estimate_gamma) V_full else V,
       model_matrix_full = Z, bread_weights = w_glm,
       score_residuals = weights * resid_raw,
       n_obs = n_obs, df.residual = n_obs - length(opt$par),
       loglik = -opt$value, theta = NULL,
       gamma_estimated = if (estimate_gamma) gamma_hat else NULL,
       constrained = constrained,
       convergence = opt$convergence,
       X_alpha = X_alpha, X_beta = X_beta)
}


# ---- Negative Binomial MLE ----

#' @noRd
.fit_nb <- function(m, N, ratio, log_N, log_rate,
                    X_alpha, X_beta, gamma_value = NULL,
                    estimate_gamma = FALSE, gamma_start = 0.01,
                    gamma_bounds = c(1e-10, 1), theta_start = 1,
                    weights = NULL, vcov_type = "HC3",
                    constrained = FALSE) {
  n_obs <- length(m)
  if (is.null(weights)) weights <- rep(1, n_obs)
  p_alpha <- ncol(X_alpha)
  p_beta <- ncol(X_beta)

  .alpha_fn <- if (constrained) function(x) .inv_logit(x) else identity
  .beta_fn  <- if (constrained) function(x) exp(x) else identity

  nll <- function(par) {
    a <- par[seq_len(p_alpha)]
    b <- par[p_alpha + seq_len(p_beta)]
    theta <- exp(par[p_alpha + p_beta + 1])
    alpha_lin <- .alpha_fn(as.numeric(X_alpha %*% a))
    idx_gamma <- p_alpha + p_beta + 2
    if (estimate_gamma) {
      g <- exp(par[idx_gamma])
      lr <- log(g + ratio)
    } else {
      lr <- log_rate
    }
    beta_lin <- .beta_fn(as.numeric(X_beta %*% b))
    log_mu <- pmin(alpha_lin * log_N + beta_lin * lr, 20)
    mu <- pmax(exp(log_mu), 1e-300)
    val <- -sum(weights * dnbinom(m, size = theta, mu = mu, log = TRUE))
    if (!is.finite(val)) return(1e20)
    val
  }

  grad_nll <- function(par) {
    a <- par[seq_len(p_alpha)]
    b <- par[p_alpha + seq_len(p_beta)]
    theta <- exp(par[p_alpha + p_beta + 1])
    eta_alpha <- as.numeric(X_alpha %*% a)
    zeta_beta <- as.numeric(X_beta %*% b)
    alpha_lin <- .alpha_fn(eta_alpha)
    beta_lin <- .beta_fn(zeta_beta)
    idx_gamma <- p_alpha + p_beta + 2
    if (estimate_gamma) {
      g <- exp(par[idx_gamma])
      rate <- g + ratio
      lr <- log(rate)
    } else {
      lr <- log_rate
      rate <- exp(lr)
    }
    log_mu <- pmin(alpha_lin * log_N + beta_lin * lr, 20)
    mu <- pmax(exp(log_mu), 1e-300)

    w_nb <- weights * (mu - m) * theta / (theta + mu)

    if (constrained) {
      dalpha <- alpha_lin * (1 - alpha_lin)
      dbeta <- beta_lin
    } else {
      dalpha <- rep(1, n_obs)
      dbeta <- rep(1, n_obs)
    }

    g_alpha <- colSums(w_nb * log_N * dalpha * X_alpha)
    g_beta <- colSums(w_nb * lr * dbeta * X_beta)

    g_theta <- -sum(weights * (
      digamma(m + theta) - digamma(theta) +
        log(theta) + 1 - log(theta + mu) -
        (m + theta) / (theta + mu)
    )) * theta

    grad <- c(g_alpha, g_beta, g_theta)
    if (estimate_gamma) {
      g_gamma <- sum(w_nb * beta_lin / rate) * g
      grad <- c(grad, g_gamma)
    }
    grad
  }

  # Starting values: first get unconstrained OLS, then transform
  has_zeros <- any(m == 0)
  y_start <- if (has_zeros) log(m + 1) else log(m)
  Z_start <- cbind(X_alpha * log_N, X_beta * log_rate)
  start_ab <- lm.fit(Z_start, y_start)$coefficients

  if (constrained) {
    # For constrained: first run unconstrained NB to get good starting values
    # then transform to link scale
    start_ab_uc <- tryCatch({
      fit_uc <- .fit_nb(m, N, ratio, log_N, log_rate, X_alpha, X_beta,
                        gamma_value = gamma_value, estimate_gamma = estimate_gamma,
                        gamma_start = gamma_start, gamma_bounds = gamma_bounds,
                        theta_start = theta_start, weights = weights,
                        vcov_type = "HC0", constrained = FALSE)
      c(fit_uc$alpha_coefs, fit_uc$beta_coefs)
    }, error = function(e) start_ab)

    start_ab[seq_len(p_alpha)] <- .logit(pmin(pmax(start_ab_uc[seq_len(p_alpha)], 0.05), 0.95))
    start_ab[p_alpha + seq_len(p_beta)] <- log(pmax(start_ab_uc[p_alpha + seq_len(p_beta)], 0.01))
  }

  start_par <- c(start_ab, log(theta_start))

  if (estimate_gamma) {
    start_par <- c(start_par, log(gamma_start))
    lower <- c(rep(-Inf, p_alpha + p_beta), log(0.01), log(gamma_bounds[1]))
    upper <- c(rep(Inf, p_alpha + p_beta), log(1e6), log(gamma_bounds[2]))
  } else {
    lower <- c(rep(-Inf, p_alpha + p_beta), log(0.01))
    upper <- c(rep(Inf, p_alpha + p_beta), log(1e6))
  }

  opt <- optim(start_par, nll, gr = grad_nll, method = "L-BFGS-B",
               lower = lower, upper = upper,
               control = list(maxit = 1000), hessian = TRUE)

  if (opt$convergence != 0) {
    warning("NB optimization did not converge (code ", opt$convergence,
            "): ", opt$message, call. = FALSE)
  }

  a_raw <- opt$par[seq_len(p_alpha)]
  b_raw <- opt$par[p_alpha + seq_len(p_beta)]
  theta <- exp(opt$par[p_alpha + p_beta + 1])
  gamma_hat <- if (estimate_gamma) exp(opt$par[p_alpha + p_beta + 2]) else gamma_value

  eta_alpha <- as.numeric(X_alpha %*% a_raw)
  zeta_beta <- as.numeric(X_beta %*% b_raw)
  alpha_lin <- .alpha_fn(eta_alpha)
  beta_lin <- .beta_fn(zeta_beta)

  if (estimate_gamma) {
    rate_final <- gamma_hat + ratio
    lr_final <- log(rate_final)
  } else {
    lr_final <- log_rate
    rate_final <- exp(lr_final)
  }
  log_mu <- alpha_lin * log_N + beta_lin * lr_final
  mu <- exp(log_mu)
  resid_raw <- m - mu

  if (constrained) {
    dalpha <- alpha_lin * (1 - alpha_lin)
    dbeta <- beta_lin
    Z <- cbind(log_N * dalpha * X_alpha, lr_final * dbeta * X_beta)
  } else {
    dalpha <- rep(1, n_obs)
    dbeta <- rep(1, n_obs)
    Z <- cbind(log_N * X_alpha, lr_final * X_beta)
  }
  if (estimate_gamma) {
    Z <- cbind(Z, gamma = beta_lin / rate_final)
  }

  w_nb <- weights * mu * theta / (theta + mu)
  h <- .hat_values_wls(Z, w_nb)

  resid_pearson <- resid_raw / sqrt(mu + mu^2 / theta)
  V_full <- .compute_sandwich_vcov(Z, resid_pearson, weights = w_nb,
                                   hat_values = h, vcov_type = vcov_type,
     )

  p_ab <- p_alpha + p_beta
  V <- if (estimate_gamma) {
    V_full[seq_len(p_ab), seq_len(p_ab), drop = FALSE]
  } else {
    V_full
  }

  V_model_full <- .compute_model_vcov(Z, w_nb)
  V_model <- if (estimate_gamma) {
    V_model_full[seq_len(p_ab), seq_len(p_ab), drop = FALSE]
  } else {
    V_model_full
  }

  ## NB score residual: w_i * (m_i - mu_i) * theta/(theta + mu_i)
  nb_score_factor <- weights * resid_raw * theta / (theta + mu)

  ## Per-observation score on optimizer (log) scale for full NB sandwich.
  ## Sign: score of log-likelihood = -score of nll, so negate grad_nll terms.
  w_nb_signed <- -weights * (mu - m) * theta / (theta + mu)
  score_alpha <- w_nb_signed * log_N * dalpha * X_alpha
  score_beta  <- w_nb_signed * lr_final * dbeta * X_beta
  score_theta <- weights * (
    digamma(m + theta) - digamma(theta) +
    log(theta) + 1 - log(theta + mu) -
    (m + theta) / (theta + mu)
  ) * theta  # on log-theta scale
  score_full <- cbind(score_alpha, score_beta, theta = score_theta)
  if (estimate_gamma) {
    score_gamma <- w_nb_signed * beta_lin / rate_final * gamma_hat
    score_full <- cbind(score_full, gamma = score_gamma)
  }

  list(alpha_coefs = a_raw, beta_coefs = b_raw,
       alpha_values = alpha_lin, beta_values = beta_lin,
       fitted = mu, residuals = resid_raw, log_mu = log_mu,
       vcov = V, vcov_model = V_model,
       vcov_full = if (estimate_gamma) V_full else V,
       model_matrix_full = Z, bread_weights = w_nb,
       score_residuals = nb_score_factor,
       score_full = score_full,
       hessian_nll = opt$hessian,
       n_obs = n_obs, df.residual = n_obs - length(opt$par),
       loglik = -opt$value, theta = theta,
       gamma_estimated = if (estimate_gamma) gamma_hat else NULL,
       constrained = constrained,
       convergence = opt$convergence,
       X_alpha = X_alpha, X_beta = X_beta)
}
