

# LOG-LIKELIHOOD
log_lik_zhang_model_cov <- function(alpha, beta, phi, m, n ,N, X, Z){
  
  mu <- N^(X %*% alpha) * (n/N)^(Z %*% beta)
  
  # check point
  if (any(!is.finite(mu)) || any(mu <= 0)) return(Inf)   # obvious why, log(mu) if mu <= 0  -> error
  if (!is.finite(phi) || phi <= 0) return(Inf)           # if phi = 0, then lgamma(phi)= Inf
  
  val <- m * log(mu) + phi * log(phi) -
    (m + phi) * log(mu + phi) +
    lgamma(m + phi) - lfactorial(m) - lgamma(phi)
  
  # check if there are any inf
  if (any(!is.finite(val))) return(Inf)               
  
  return(sum(val))
}


# GRADIENT
grad_log_lik_zhang_model_cov <- function(alpha, beta, phi, m, n ,N, X, Z){
  
  mu <- N^(X %*% alpha)*(n/N)^(Z %*% beta)
  logN <- log(N)
  logRatio <- log(n/N)
  
  # check point
  if (any(mu <= 0) || any(!is.finite(mu)) ||   # obvious
      phi <= 0 || !is.finite(phi)) {           # obvious
    return(rep(NA, ncol(X) + ncol(Z) + 1))     # to have correct gradient dimension and still try to go into optim()
  } 
  
  d_mu <- m/mu - (m + phi)/(phi + mu)
  d_alpha <- d_mu * mu * logN
  d_beta <- d_mu * mu * logRatio
  
  log_mu_phi <- log(mu + phi)
  log_m_phi  <- log(m + phi)
  
  # check point
  if (any(!is.finite(log_mu_phi)) ||        # to avoid problems in sum d_phi
      any(!is.finite(log_m_phi))) {
    return(rep(NA, ncol(X) + ncol(Z) + 1))  # to have correct gradient dimension and still try to go into optim()
  } 
  
  d_phi <- 1/(2*phi) - m/(mu + phi) - log(mu + phi) 
  - phi/(mu + phi) + log(m + phi) + (m + phi - 0.5)/(m + phi)
  
  c(sum(t(X) %*% d_alpha), sum(t(Z) %*% d_beta), sum(d_phi))  
  
}


# HESSIAN
hess_log_lik_zhang_model_cov <- function(alpha, beta, phi, m, n ,N, X, Z){
  
  mu <- N^(X %*% alpha)*(n/N)^(Z %*% beta)
  logN <- log(N)
  logRatio <- log(n/N)
  
  # check point
  if (any(mu <= 0) || phi <= 0 || any(!is.finite(mu))) {
    return(matrix(NA, nrow = ncol(X) + ncol(Z) + 1, ncol = ncol(X) + ncol(Z) + 1))
  }          # similar to gradient check point
  
  d_mu <- m/mu - (m + phi)/(phi + mu)
  d_mu2 <- -m/mu^2 + (m + phi)/(mu + phi)^2
  
  d_alpha2 <- d_mu2 * (mu * logN)^2 + d_mu * mu * logN^2
  d_beta2 <- d_mu2 * (mu * logRatio)^2 + d_mu * mu * logRatio^2
  d_phi2 <- 1/(2*phi^2) - (2*mu - m + phi)/(mu + phi)^2 + (mu + phi+ 0.5)/(m + phi)^2
  
  d_alpha_beta <- d_mu2 * mu^2 * logN * logRatio + d_mu * mu * logN * logRatio
  d_alpha_phi <- mu * logN * (m - mu)/(mu + phi)^2  
  d_beta_phi <- mu * logRatio * (m - mu)/(mu + phi)^2  
  
  # include covariates
  H11 <- t(X) %*% (d_alpha2 * X)      ##### SPRAWDZIĆ i PRZEMYŚLEĆ TO !
  H12 <- t(X) %*% (d_alpha_beta * Z)
  H13 <- t(X) %*% d_alpha_phi
  
  H22 <- t(Z) %*% (d_beta2 * X)
  H23 <- t(Z) %*% d_beta_phi
  
  H33 <- matrix(sum(d_phi2), ncol=1, nrow=1)
  
  # full hessian matrix
  trow <- cbind(H11, H12, H13)
  mrow <- cbind(H12, H22, H23)
  brow <- cbind(H13, H23, H33)
  rbind(trow, mrow, brow)
  
}

