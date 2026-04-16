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
                     weights = NULL, vcov_type = "HC3",
                     link_rho = "power") {
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
      rate <- g + ratio
    } else {
      rate <- exp(log_rate)
    }
    beta_lin <- as.numeric(X_beta %*% b)
    log_mu <- pmin(.compute_log_mu(alpha_lin, log_N, beta_lin, rate,
                                   link_rho = link_rho), 20)
    mu <- exp(log_mu)
    sum(weights * (m - mu)^2)
  }

  grad_fn <- function(par) {
    a <- par[seq_len(p_alpha)]
    b <- par[p_alpha + seq_len(p_beta)]
    alpha_lin <- as.numeric(X_alpha %*% a)
    if (estimate_gamma) {
      g <- exp(par[p_alpha + p_beta + 1])
      rate <- g + ratio
    } else {
      rate <- exp(log_rate)
    }
    beta_lin <- as.numeric(X_beta %*% b)
    eta <- .eta_from_rate(beta_lin, rate)
    dlogrho <- .dlog_rho_deta(eta, link_rho = link_rho)
    log_mu <- pmin(.compute_log_mu(alpha_lin, log_N, beta_lin, rate,
                                   link_rho = link_rho), 20)
    mu <- exp(log_mu)
    resid <- m - mu
    g_alpha <- -2 * colSums(weights * resid * mu * log_N * X_alpha)
    g_beta <- -2 * colSums(weights * resid * mu * (dlogrho * log(rate)) * X_beta)
    if (estimate_gamma) {
      g_gamma <- -2 * sum(weights * resid * mu * dlogrho * beta_lin * g / rate)
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
  rate_final <- if (estimate_gamma) gamma_hat + ratio else exp(log_rate)
  beta_lin <- as.numeric(X_beta %*% b)
  eta_final <- .eta_from_rate(beta_lin, rate_final)
  dlogrho_final <- .dlog_rho_deta(eta_final, link_rho = link_rho)
  log_mu <- .compute_log_mu(alpha_lin, log_N, beta_lin, rate_final,
                            link_rho = link_rho)
  mu <- exp(log_mu)
  resid_raw <- m - mu

  J <- cbind(mu * log_N * X_alpha,
             mu * (dlogrho_final * log(rate_final)) * X_beta)
  if (estimate_gamma) {
    J <- cbind(J, gamma = mu * dlogrho_final * beta_lin * gamma_hat / rate_final)
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
                         X_gamma = NULL,
                         weights = NULL, vcov_type = "HC3",
                         constrained = FALSE,
                         link_rho = "power") {
  n_obs <- length(m)
  if (is.null(weights)) weights <- rep(1, n_obs)
  p_alpha <- ncol(X_alpha)
  p_beta <- ncol(X_beta)
  p_gamma <- if (!is.null(X_gamma)) ncol(X_gamma) else 0L

  # Transform functions for constrained optimization
  # alpha = inverse logit(eta) -> (0,1); beta = exp(zeta) -> (0, Inf)
  .alpha_fn <- if (constrained) function(x) .inv_logit(x) else identity
  .beta_fn  <- if (constrained) function(x) exp(x) else identity

  nll <- function(par) {
    a <- par[seq_len(p_alpha)]
    b <- par[p_alpha + seq_len(p_beta)]
    alpha_lin <- .alpha_fn(as.numeric(X_alpha %*% a))
    if (estimate_gamma) {
      gamma_coefs_cur <- par[p_alpha + p_beta + seq_len(p_gamma)]
      gamma_lin <- exp(pmin(as.numeric(X_gamma %*% gamma_coefs_cur), 10))
      rate <- gamma_lin + ratio
    } else {
      rate <- exp(log_rate)
    }
    beta_lin <- .beta_fn(as.numeric(X_beta %*% b))
    log_mu <- pmin(.compute_log_mu(alpha_lin, log_N, beta_lin, rate,
                                   link_rho = link_rho), 20)
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
      gamma_coefs_cur <- par[p_alpha + p_beta + seq_len(p_gamma)]
      gamma_lin <- exp(pmin(as.numeric(X_gamma %*% gamma_coefs_cur), 10))
      rate <- gamma_lin + ratio
    } else {
      rate <- exp(log_rate)
    }
    eta <- .eta_from_rate(beta_lin, rate)
    dlogrho <- .dlog_rho_deta(eta, link_rho = link_rho)
    log_mu <- pmin(.compute_log_mu(alpha_lin, log_N, beta_lin, rate,
                                   link_rho = link_rho), 20)
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
    g_beta <- colSums(score * (dlogrho * log(rate)) * dbeta * X_beta)
    if (estimate_gamma) {
      dgamma_dcoef <- (gamma_lin / rate) * X_gamma  # n x p_gamma
      g_gamma <- colSums(score * dlogrho * beta_lin * dgamma_dcoef)
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
    gamma_start_coefs <- rep(0, p_gamma)
    gamma_start_coefs[1] <- log(gamma_start)
    start_par <- c(start_par, gamma_start_coefs)
    lower_gamma <- c(log(gamma_bounds[1]), rep(-Inf, max(0, p_gamma - 1)))
    upper_gamma <- c(log(gamma_bounds[2]), rep(Inf, max(0, p_gamma - 1)))
    lower <- c(rep(-Inf, p_alpha + p_beta), lower_gamma)
    upper <- c(rep(Inf, p_alpha + p_beta), upper_gamma)
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

  # Extract gamma on response scale
  if (estimate_gamma) {
    gamma_coefs_hat <- opt$par[p_alpha + p_beta + seq_len(p_gamma)]
    gamma_lin_hat <- exp(pmin(as.numeric(X_gamma %*% gamma_coefs_hat), 10))
    gamma_hat <- if (p_gamma == 1) gamma_lin_hat[1] else NULL
    names(gamma_coefs_hat) <- colnames(X_gamma)
  } else {
    gamma_coefs_hat <- NULL
    gamma_lin_hat <- NULL
    gamma_hat <- gamma_value
  }

  # Transform to response scale
  eta_alpha <- as.numeric(X_alpha %*% a_raw)
  zeta_beta <- as.numeric(X_beta %*% b_raw)
  alpha_lin <- .alpha_fn(eta_alpha)
  beta_lin <- .beta_fn(zeta_beta)

  if (estimate_gamma) {
    rate_final <- gamma_lin_hat + ratio
  } else {
    rate_final <- exp(log_rate)
  }
  eta_final <- .eta_from_rate(beta_lin, rate_final)
  dlogrho_final <- .dlog_rho_deta(eta_final, link_rho = link_rho)
  log_mu <- .compute_log_mu(alpha_lin, log_N, beta_lin, rate_final,
                            link_rho = link_rho)
  mu <- exp(log_mu)
  resid_raw <- m - mu

  # Design matrix on link scale (d(log_mu)/d(eta,zeta))
  if (constrained) {
    dalpha <- alpha_lin * (1 - alpha_lin)
    dbeta <- beta_lin
    Z <- cbind(log_N * dalpha * X_alpha,
               (dlogrho_final * log(rate_final)) * dbeta * X_beta)
  } else {
    Z <- cbind(log_N * X_alpha,
               (dlogrho_final * log(rate_final)) * X_beta)
  }
  if (estimate_gamma) {
    gamma_cols <- (dlogrho_final * beta_lin * gamma_lin_hat / rate_final) * X_gamma
    colnames(gamma_cols) <- colnames(X_gamma)
    Z <- cbind(Z, gamma_cols)
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

  ## Z is the full model matrix (including gamma cols if estimated)
  list(alpha_coefs = a_raw, beta_coefs = b_raw,
       alpha_values = alpha_lin, beta_values = beta_lin,
       gamma_coefs = gamma_coefs_hat, gamma_values = gamma_lin_hat,
       fitted = mu, residuals = resid_raw, log_mu = log_mu,
       vcov = V, vcov_model = V_model,
       vcov_full = if (estimate_gamma) V_full else V,
       model_matrix_full = Z, bread_weights = w_glm,
       score_residuals = weights * resid_raw,
       n_obs = n_obs, df.residual = n_obs - length(opt$par),
       loglik = -opt$value, theta = NULL,
       gamma_estimated = if (estimate_gamma) {
         gamma_lin_hat[1]
       } else NULL,
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
                    X_gamma = NULL,
                    weights = NULL, vcov_type = "HC3",
                    constrained = FALSE,
                    link_rho = "power") {
  n_obs <- length(m)
  if (is.null(weights)) weights <- rep(1, n_obs)
  p_alpha <- ncol(X_alpha)
  p_beta <- ncol(X_beta)
  p_gamma <- if (!is.null(X_gamma)) ncol(X_gamma) else 0L

  .alpha_fn <- if (constrained) function(x) .inv_logit(x) else identity
  .beta_fn  <- if (constrained) function(x) exp(x) else identity

  # NB parameter layout: c(alpha[p_a], beta[p_b], log_theta[1], gamma[p_g])
  idx_theta <- p_alpha + p_beta + 1
  idx_gamma_start <- p_alpha + p_beta + 2

  nll <- function(par) {
    a <- par[seq_len(p_alpha)]
    b <- par[p_alpha + seq_len(p_beta)]
    theta <- exp(par[idx_theta])
    alpha_lin <- .alpha_fn(as.numeric(X_alpha %*% a))
    if (estimate_gamma) {
      gamma_coefs_cur <- par[idx_gamma_start - 1 + seq_len(p_gamma)]
      gamma_lin <- exp(pmin(as.numeric(X_gamma %*% gamma_coefs_cur), 10))
      rate <- gamma_lin + ratio
    } else {
      rate <- exp(log_rate)
    }
    beta_lin <- .beta_fn(as.numeric(X_beta %*% b))
    log_mu <- pmin(.compute_log_mu(alpha_lin, log_N, beta_lin, rate,
                                   link_rho = link_rho), 20)
    mu <- pmax(exp(log_mu), 1e-300)
    val <- -sum(weights * dnbinom(m, size = theta, mu = mu, log = TRUE))
    if (!is.finite(val)) return(1e20)
    val
  }

  grad_nll <- function(par) {
    a <- par[seq_len(p_alpha)]
    b <- par[p_alpha + seq_len(p_beta)]
    theta <- exp(par[idx_theta])
    eta_alpha <- as.numeric(X_alpha %*% a)
    zeta_beta <- as.numeric(X_beta %*% b)
    alpha_lin <- .alpha_fn(eta_alpha)
    beta_lin <- .beta_fn(zeta_beta)
    if (estimate_gamma) {
      gamma_coefs_cur <- par[idx_gamma_start - 1 + seq_len(p_gamma)]
      gamma_lin <- exp(pmin(as.numeric(X_gamma %*% gamma_coefs_cur), 10))
      rate <- gamma_lin + ratio
    } else {
      rate <- exp(log_rate)
    }
    eta <- .eta_from_rate(beta_lin, rate)
    dlogrho <- .dlog_rho_deta(eta, link_rho = link_rho)
    log_mu <- pmin(.compute_log_mu(alpha_lin, log_N, beta_lin, rate,
                                   link_rho = link_rho), 20)
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
    g_beta <- colSums(w_nb * (dlogrho * log(rate)) * dbeta * X_beta)

    g_theta <- -sum(weights * (
      digamma(m + theta) - digamma(theta) +
        log(theta) + 1 - log(theta + mu) -
        (m + theta) / (theta + mu)
    )) * theta

    grad <- c(g_alpha, g_beta, g_theta)
    if (estimate_gamma) {
      dgamma_dcoef <- (gamma_lin / rate) * X_gamma  # n x p_gamma
      g_gamma <- colSums(w_nb * dlogrho * beta_lin * dgamma_dcoef)
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
                        theta_start = theta_start, X_gamma = X_gamma,
                        weights = weights,
                        vcov_type = "HC0", constrained = FALSE,
                        link_rho = link_rho)
      c(fit_uc$alpha_coefs, fit_uc$beta_coefs)
    }, error = function(e) start_ab)

    start_ab[seq_len(p_alpha)] <- .logit(pmin(pmax(start_ab_uc[seq_len(p_alpha)], 0.05), 0.95))
    start_ab[p_alpha + seq_len(p_beta)] <- log(pmax(start_ab_uc[p_alpha + seq_len(p_beta)], 0.01))
  }

  start_par <- c(start_ab, log(theta_start))

  if (estimate_gamma) {
    gamma_start_coefs <- rep(0, p_gamma)
    gamma_start_coefs[1] <- log(gamma_start)
    start_par <- c(start_par, gamma_start_coefs)
    lower_gamma <- c(log(gamma_bounds[1]), rep(-Inf, max(0, p_gamma - 1)))
    upper_gamma <- c(log(gamma_bounds[2]), rep(Inf, max(0, p_gamma - 1)))
    lower <- c(rep(-Inf, p_alpha + p_beta), log(0.01), lower_gamma)
    upper <- c(rep(Inf, p_alpha + p_beta), log(1e6), upper_gamma)
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
  theta <- exp(opt$par[idx_theta])

  # Extract gamma on response scale
  if (estimate_gamma) {
    gamma_coefs_hat <- opt$par[idx_gamma_start - 1 + seq_len(p_gamma)]
    gamma_lin_hat <- exp(pmin(as.numeric(X_gamma %*% gamma_coefs_hat), 10))
    gamma_hat <- if (p_gamma == 1) gamma_lin_hat[1] else NULL
    names(gamma_coefs_hat) <- colnames(X_gamma)
  } else {
    gamma_coefs_hat <- NULL
    gamma_lin_hat <- NULL
    gamma_hat <- gamma_value
  }

  eta_alpha <- as.numeric(X_alpha %*% a_raw)
  zeta_beta <- as.numeric(X_beta %*% b_raw)
  alpha_lin <- .alpha_fn(eta_alpha)
  beta_lin <- .beta_fn(zeta_beta)

  if (estimate_gamma) {
    rate_final <- gamma_lin_hat + ratio
  } else {
    rate_final <- exp(log_rate)
  }
  eta_final <- .eta_from_rate(beta_lin, rate_final)
  dlogrho_final <- .dlog_rho_deta(eta_final, link_rho = link_rho)
  log_mu <- .compute_log_mu(alpha_lin, log_N, beta_lin, rate_final,
                            link_rho = link_rho)
  mu <- exp(log_mu)
  resid_raw <- m - mu

  if (constrained) {
    dalpha <- alpha_lin * (1 - alpha_lin)
    dbeta <- beta_lin
    Z <- cbind(log_N * dalpha * X_alpha,
               (dlogrho_final * log(rate_final)) * dbeta * X_beta)
  } else {
    dalpha <- rep(1, n_obs)
    dbeta <- rep(1, n_obs)
    Z <- cbind(log_N * X_alpha,
               (dlogrho_final * log(rate_final)) * X_beta)
  }
  if (estimate_gamma) {
    gamma_cols <- (dlogrho_final * beta_lin * gamma_lin_hat / rate_final) * X_gamma
    colnames(gamma_cols) <- colnames(X_gamma)
    Z <- cbind(Z, gamma_cols)
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
  score_beta  <- w_nb_signed * (dlogrho_final * log(rate_final)) * dbeta * X_beta
  score_theta <- weights * (
    digamma(m + theta) - digamma(theta) +
    log(theta) + 1 - log(theta + mu) -
    (m + theta) / (theta + mu)
  ) * theta  # on log-theta scale
  score_full <- cbind(score_alpha, score_beta, theta = score_theta)
  if (estimate_gamma) {
    # Per-observation score for gamma coefficients
    dgamma_dcoef_final <- (gamma_lin_hat / rate_final) * X_gamma
    score_gamma_mat <- w_nb_signed * dlogrho_final * beta_lin * dgamma_dcoef_final
    colnames(score_gamma_mat) <- colnames(X_gamma)
    score_full <- cbind(score_full, score_gamma_mat)
  }

  list(alpha_coefs = a_raw, beta_coefs = b_raw,
       alpha_values = alpha_lin, beta_values = beta_lin,
       gamma_coefs = gamma_coefs_hat, gamma_values = gamma_lin_hat,
       fitted = mu, residuals = resid_raw, log_mu = log_mu,
       vcov = V, vcov_model = V_model,
       vcov_full = if (estimate_gamma) V_full else V,
       model_matrix_full = Z, bread_weights = w_nb,
       score_residuals = nb_score_factor,
       score_full = score_full,
       hessian_nll = opt$hessian,
       n_obs = n_obs, df.residual = n_obs - length(opt$par),
       loglik = -opt$value, theta = theta,
       gamma_estimated = if (estimate_gamma) {
         gamma_lin_hat[1]
       } else NULL,
       constrained = constrained,
       convergence = opt$convergence,
       X_alpha = X_alpha, X_beta = X_beta)
}


# ---- iOLS estimator (Benatia, Bellego & Pape, 2024) ----
#
# Two-phase iterated OLS targeting GPML (Gamma PML) score equations.
# Phase 1: warm-up with increasing delta and empirical centering.
# Phase 2: limiting GPML transformation to remove O(delta^-2) bias.
# At convergence, beta solves: Z'(m/mu - 1) = 0 (GPML first-order conditions).

#' @noRd
.fit_iols <- function(m, N, ratio, log_N, log_rate,
                      X_alpha, X_beta, gamma_value = NULL,
                      estimate_gamma = FALSE, gamma_start = 0.01,
                      gamma_bounds = c(1e-10, 0.5),
                      weights = NULL, vcov_type = "HC3",
                      delta_grid = c(1, 10, 100, 1000),
                      rho = 1,
                      tol = 1e-6, max_iter = 100) {
  n_obs <- length(m)
  if (is.null(weights)) weights <- rep(1, n_obs)
  p_alpha <- ncol(X_alpha)
  p_beta <- ncol(X_beta)

  Z <- cbind(X_alpha * log_N, X_beta * log_rate)
  p <- ncol(Z)

  # Initialize with log(m+1) OLS
  y_start <- log(m + 1)
  beta_hat <- if (all(weights == 1)) {
    lm.fit(Z, y_start)$coefficients
  } else {
    lm.wfit(Z, y_start, weights)$coefficients
  }

  .ols_step <- function(Z, y, w) {
    if (all(w == 1)) lm.fit(Z, y)$coefficients
    else lm.wfit(Z, y, w)$coefficients
  }

  # Phase 1: warm-up with increasing delta and empirical centering
  for (delta in delta_grid) {
    for (iter in 1:max_iter) {
      log_mu <- pmin(as.numeric(Z %*% beta_hat), 20)
      mu <- exp(log_mu)
      lhs <- log(m + delta * mu)
      c_delta <- weighted.mean(lhs, weights)
      y_tilde <- lhs - c_delta
      beta_new <- .ols_step(Z, y_tilde, weights)
      if (max(abs(beta_new - beta_hat)) < tol) {
        beta_hat <- beta_new
        break
      }
      beta_hat <- beta_new
    }
  }

  # Phase 2: exact GPML limiting transformation
  # y_tilde = log(mu) + (m/mu - 1) / (1 + rho)
  # At convergence: Z'(m/mu - 1) = 0 (GPML score equations)
  # Try increasing rho until contraction condition is met
  for (rho_try in c(rho, 2, 5, 10, 50, 100)) {
    beta_ph2 <- beta_hat
    converged_ph2 <- FALSE
    for (iter in 1:max_iter) {
      log_mu <- pmin(pmax(as.numeric(Z %*% beta_ph2), -20), 20)
      mu <- exp(log_mu)
      u <- m / pmax(mu, 1e-300)
      y_tilde <- log_mu + (u - 1) / (1 + rho_try)
      beta_new <- .ols_step(Z, y_tilde, weights)
      if (max(abs(beta_new - beta_ph2)) < tol) {
        beta_ph2 <- beta_new
        converged_ph2 <- TRUE
        break
      }
      if (max(abs(beta_new)) > 100) break
      beta_ph2 <- beta_new
    }
    if (converged_ph2) {
      beta_hat <- beta_ph2
      break
    }
  }

  # Final quantities
  coefs <- beta_hat
  names(coefs) <- colnames(Z)
  alpha_coefs <- coefs[seq_len(p_alpha)]
  beta_coefs <- coefs[p_alpha + seq_len(p_beta)]

  log_mu <- as.numeric(Z %*% coefs)
  mu <- exp(log_mu)
  resid_raw <- m - mu

  # Check GPML score condition: Z'(m/mu - 1) / n should be near zero
  gpml_resid <- as.numeric(m / pmax(mu, 1e-300) - 1)
  score_mean <- as.numeric(crossprod(Z * weights, gpml_resid)) / sum(weights)
  score_tol <- max(tol * 10, 1e-4)  # more generous than iteration tol
  convergence <- if (max(abs(score_mean)) < score_tol) 0L else 1L

  if (convergence != 0L) {
    warning("iOLS did not converge to GPML solution (max |score/n| = ",
            round(max(abs(score_mean)), 6), ").", call. = FALSE)
  }

  # GPML deviance as pseudo log-likelihood
  loglik_gpml <- -sum(weights * (m / mu - log(pmax(m, 1e-300) / mu) - 1))

  # Sandwich vcov using GPML score residual: (m/mu - 1)
  h <- .hat_values_wls(Z, weights)
  V <- .compute_sandwich_vcov(Z, gpml_resid, weights = weights,
                              hat_values = h, vcov_type = vcov_type)

  # GPML Fisher information: (Z'Z)^{-1} without sigma2 scaling.
  # Gamma variance function gives unit working weights, so the model-based
  # variance is (Z'Z)^{-1}, not sigma2*(Z'Z)^{-1}. See notes/bias-correction-iols.md.
  V_model <- .compute_model_vcov(Z, weights, sigma2 = NULL)
  sigma2 <- sum(weights * gpml_resid^2) / (n_obs - p)

  list(alpha_coefs = alpha_coefs, beta_coefs = beta_coefs,
       fitted = mu, residuals = resid_raw, log_mu = log_mu,
       vcov = V, vcov_model = V_model,
       model_matrix_full = Z, bread_weights = weights,
       score_residuals = weights * gpml_resid,
       n_obs = n_obs, df.residual = n_obs - p,
       loglik = loglik_gpml, theta = NULL,
       gamma_estimated = NULL,
       sigma2 = sigma2, convergence = convergence)
}


# ---- iOLS with profiled gamma ----

#' @noRd
.fit_iols_gamma <- function(m, N, ratio, log_N, X_alpha, X_beta,
                            gamma_start, gamma_bounds,
                            weights = NULL, vcov_type = "HC3",
                            delta_grid = c(1, 10, 100, 1000),
                            rho = 1,
                            tol = 1e-6, max_iter = 100) {
  n_obs <- length(m)
  if (is.null(weights)) weights <- rep(1, n_obs)
  p_alpha <- ncol(X_alpha)
  p_beta <- ncol(X_beta)

  .ols_step <- function(Z, y, w) {
    if (all(w == 1)) lm.fit(Z, y)$coefficients
    else lm.wfit(Z, y, w)$coefficients
  }

  # GPML deviance for gamma profiling
  .gpml_deviance <- function(m, mu, w) {
    sum(w * (m / pmax(mu, 1e-300) - log(pmax(m, 1e-300) / pmax(mu, 1e-300)) - 1))
  }

  # Initialize
  log_rate_init <- log(gamma_start + ratio)
  Z_init <- cbind(X_alpha * log_N, X_beta * log_rate_init)
  beta_hat <- .ols_step(Z_init, log(m + 1), weights)
  gamma_hat <- gamma_start

  # Phase 1: warm-up with increasing delta, profiling gamma via GPML deviance
  for (delta in delta_grid) {
    for (iter in 1:max_iter) {
      # Profile gamma
      obj_fn <- function(g) {
        lr <- log(g + ratio)
        Z_g <- cbind(X_alpha * log_N, X_beta * lr)
        mu <- exp(as.numeric(Z_g %*% beta_hat))
        .gpml_deviance(m, mu, weights)
      }
      gamma_hat <- optimize(obj_fn, c(gamma_bounds[1], gamma_bounds[2]))$minimum

      log_rate <- log(gamma_hat + ratio)
      Z <- cbind(X_alpha * log_N, X_beta * log_rate)
      mu <- exp(as.numeric(Z %*% beta_hat))
      lhs <- log(m + delta * mu)
      c_delta <- weighted.mean(lhs, weights)
      y_tilde <- lhs - c_delta
      beta_new <- .ols_step(Z, y_tilde, weights)
      if (max(abs(beta_new - beta_hat)) < tol) {
        beta_hat <- beta_new
        break
      }
      beta_hat <- beta_new
    }
  }

  # Phase 2: GPML limiting transformation with gamma profiling
  for (iter in 1:max_iter) {
    # Profile gamma
    obj_fn <- function(g) {
      lr <- log(g + ratio)
      Z_g <- cbind(X_alpha * log_N, X_beta * lr)
      mu <- exp(as.numeric(Z_g %*% beta_hat))
      .gpml_deviance(m, mu, weights)
    }
    gamma_hat <- optimize(obj_fn, c(gamma_bounds[1], gamma_bounds[2]))$minimum

    log_rate <- log(gamma_hat + ratio)
    Z <- cbind(X_alpha * log_N, X_beta * log_rate)
    mu <- exp(as.numeric(Z %*% beta_hat))
    u <- m / pmax(mu, 1e-300)
    c_inf <- log(rho) + (1 / (1 + rho)) * (u - 1)
    y_tilde <- log(m + rho * mu) - c_inf
    beta_new <- .ols_step(Z, y_tilde, weights)
    if (max(abs(beta_new - beta_hat)) < tol) {
      beta_hat <- beta_new
      break
    }
    beta_hat <- beta_new
  }

  # Final quantities
  log_rate <- log(gamma_hat + ratio)
  Z <- cbind(X_alpha * log_N, X_beta * log_rate)
  coefs <- beta_hat
  names(coefs) <- colnames(Z)
  alpha_coefs <- coefs[seq_len(p_alpha)]
  beta_coefs <- coefs[p_alpha + seq_len(p_beta)]

  log_mu <- as.numeric(Z %*% coefs)
  mu <- exp(log_mu)
  resid_raw <- m - mu

  gpml_resid <- as.numeric(m / pmax(mu, 1e-300) - 1)
  gpml_score <- as.numeric(crossprod(Z * weights, gpml_resid))
  convergence <- if (max(abs(gpml_score)) < 0.1) 0L else 1L

  loglik_gpml <- -sum(weights * (m / mu - log(pmax(m, 1e-300) / mu) - 1))

  # Sandwich with gamma column
  beta_lin <- as.numeric(X_beta %*% beta_coefs)
  J_gamma <- beta_lin / (gamma_hat + ratio)
  J_full <- cbind(Z, gamma = J_gamma)
  p_ab <- p_alpha + p_beta

  h <- .hat_values_wls(J_full, weights)
  V_full <- .compute_sandwich_vcov(J_full, gpml_resid, weights = weights,
                                   hat_values = h, vcov_type = vcov_type)
  V <- V_full[seq_len(p_ab), seq_len(p_ab), drop = FALSE]

  sigma2 <- sum(weights * gpml_resid^2) / (n_obs - p_ab - 1)
  V_model_full <- .compute_model_vcov(J_full, weights, sigma2 = sigma2)
  V_model <- V_model_full[seq_len(p_ab), seq_len(p_ab), drop = FALSE]

  list(alpha_coefs = alpha_coefs, beta_coefs = beta_coefs,
       fitted = mu, residuals = resid_raw, log_mu = log_mu,
       vcov = V, vcov_model = V_model, vcov_full = V_full,
       model_matrix_full = J_full, bread_weights = weights,
       score_residuals = weights * gpml_resid,
       n_obs = n_obs, df.residual = n_obs - p_ab - 1,
       loglik = loglik_gpml, theta = NULL,
       gamma_estimated = gamma_hat,
       sigma2 = sigma2, convergence = convergence)
}
