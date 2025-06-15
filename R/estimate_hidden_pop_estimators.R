

### OLS MODEL ###

ols_model <- function(m, n, N, vcov = 'hessian'){
  
  df <- data.frame(
    y = m,
    x1 = log(N),
    x2 = log(n/N)
  )
  
  if (any(df$y <= 0)) {
    message('Variable observed included zeros or negative values. They were removed.')
    df <- subset(df, y > 0)
  } 
  ols_fit <- lm(log(y) ~ x1 + x2 - 1, data = df_ols)
  alpha_est <- unname(coef(ols_fit)[1])
  beta_est <- unname(coef(ols_fit)[2])
  xi_est <- sum(N^alpha_est)
  
  X <- model.matrix(ols_fit)
  
  # covariance matrix
  if (vcov == 'robust') {
    cov_matrix <- vcovHC(ols_fit, type = "HC1")
  } else {
    cov_matrix <- tryCatch(solve(t(X)%*%X), error = function(e) NULL)
  }
  
  # standard error for alpha_est
  var_alpha <- cov_matrix[1,1]
  st_alpha <- sqrt(var_alpha)

  # standard error for xi    -  method similar to confidence intervals
  st_xi <- sum(as.numeric(N)^as.vector(st_alpha))
  
  #confidence intervals for alpha
  confint_alpha <- confint(ols_fit)[1,]
  
  #confidence intervals for xi - M estimate
  confint_xi <- c(sum(N^confint_alpha[1]), sum(N^confint_alpha[2]))
  
  estimates <- c(xi = xi_est, alpha = alpha_est, beta = beta_est)
  
  return(
    list(estimates = estimates,
         cov_matrix = cov_matrix,
         confint_alpha = setNames(confint_alpha, c("lower", "upper")),
         confint_xi = setNames(confint_xi, c("lower", "upper")),
         st_error_alpha = st_alpha,
         st_error_xi = st_xi,
         ols_fit = ols_fit
    )
  )
}





### NLS MODEL with estimates form ols as starting points ###

nls_model <- function(m, n, N, vcov = 'hessian'){
  
  df <- data.frame(
    y = m,
    N = N,
    n = n
  )
  
  if (any(df$y <= 0)) {
    message('Variable observed included zeros or negative values. They were removed.')
  }
  df_ols <- subset(df, y > 0)
  df_ols$x1 <- log(df_ols$N)
  df_ols$x2 <- log(df_ols$n/df_ols$N)
  ols_fit <- lm(log(y) ~ x1 + x2 - 1, data = df_ols)
  estim_ols <- list(alpha = coef(ols_fit)[1], beta = coef(ols_fit)[2])
  if (coef(ols_fit)[1] < 0 || coef(ols_fit)[2] < 0) {
    warning("OLS produced negative starting values. They were changed to proceed with nls.")
  }
  estim_ols$alpha <- min(max(estim_ols$alpha, 0.01), 1)
  estim_ols$beta <- min(max(estim_ols$beta, 0.01), 1)
  
  df_nls <- subset(df, y>0)
  nls_fit <- nls(y ~ N^alpha * (n/N)^beta, 
                 data = df_nls, 
                 start = estim_ols)
  estim_nls <- coef(nls_fit)
  alpha_est <- estim_nls['alpha']
  beta_est <- estim_nls['beta']
  xi_est <- sum(N^alpha_est)
  
  X <- model.matrix(nls_fit)
  
  # covariance matrix
  if (vcov == 'robust') {
    cov_matrix <- vcovHC(nls_fit, type = "HC1")  ### ERROR
  } else {
    cov_matrix <- vcov(nls_fit)                  ### ERROR
  }
  
  # standard error for alpha_est
  var_alpha <- cov_matrix[1,1]
  st_alpha <- sqrt(var_alpha)
  
  # standard error for xi
  st_xi <- sum(as.numeric(N)^as.vector(st_alpha))
  
  #confidence intervals for alpha
  confint_alpha <- confint(nls_fit)[1,]
  #confidence intervals for xi - M estimate
  confint_xi <- c(sum(N^confint_alpha[1]), sum(N^confint_alpha[2]))
  
  estimates <- c(xi = xi_est, alpha_est, beta_est)
  
  return(
    list(estimates = estimates,
         cov_matrix = cov_matrix,
         confint_alpha = setNames(confint_alpha, c("lower", "upper")),
         confint_xi = setNames(confint_xi, c("lower", "upper")),
         st_error_alpha = st_alpha,
         st_error_xi = st_xi,
         nls_fit = nls_fit
    )
  )
  
}





### GLM - ZHANG MODEL ###

zhang_model_cov <- function(m, n, N, X = NULL, Z = NULL, vcov = 'hessian', family = 'poisson'){
  
  if (is.null(X)==TRUE)  {X <- matrix(1, nrow = length(n), ncol = 1)}
  if (is.null(Z)==TRUE)  {Z <- matrix(1, nrow = length(n), ncol = 1)}
  
  p1 <- ncol(X)
  p2 <- ncol(Z)
  
  df <- data.frame(
    y = m,
    x1 = log(N),
    x2 = log(n/N)
  )
  
  if (any(df$y <= 0)) {
    message('Variable observed included zeros or negative values. They were removed.')
    df <- subset(df, y > 0)
  }
  
  if (family == 'poisson'){
    glm_fit <- glm(y ~ x1 + x2 - 1, data = df, family = poisson(link = 'log'))
  }
  starting <- coef(glm_fit)
  if (any(!is.finite(starting))) stop('GLM - infinite values')
  
  alpha_start <- rep(starting[1], p1)  # if X, Z are given, alpha and beta must be vectors of correct length
  beta_start <- rep(starting[2], p2)
  phi_start <- 0.01
  start_par <- c(alpha_start, beta_start, phi_start)
  
  
  # log-likelihood, gradient and hessian
  
  log_lik <- function(param){     # we will take start_par as param in optim function
    alpha <- param[1:p1]
    beta <- param[(p1+1):(p1+p2)]
    phi <- param[p1+p2+1]
    -log_lik_zhang_model_cov(alpha, beta, phi, m, n, N, X, Z)
  }
  
  grad_log_lik <- function(param){
    alpha <- param[1:p1]
    beta <- param[(p1+1):(p1+p2)]
    phi <- param[p1+p2+1]
    -grad_log_lik_zhang_model_cov(alpha, beta, phi, m, n, N, X, Z)
  }
  
  hess_log_lik <- function(param){
    alpha <- param[1:p1]
    beta <- param[(p1+1):(p1+p2)]
    phi <- param[p1+p2+1]
    -hess_log_lik_zhang_model_cov(alpha, beta, phi, m, n, N, X, Z)
  }
  
  # optimization
  
  optimization <- optim(
    par = start_par,
    fn = log_lik,
    gr = grad_log_lik,
    method = 'BFGS',
    hessian = TRUE
  )
  
  alpha_est <- unname(optimization$par[1:p1])
  beta_est <- unname(optimization$par[(p1+1):(p1+p2)])
  phi_est <- unname(optimization$par[p1+p2+1])
  xi_est <- sum(as.numeric(N)^as.vector(X %*% alpha_est))     # target parameter estimator
  
  estimates <- list(xi = xi_est, alpha = alpha_est, beta = beta_est, phi = phi_est)
  
  hessian <- optimization$hessian
  
  # covariance matrix
  if (vcov == 'robust') {
    cov_matrix <- vcovHC(glm_fit, type = "HC1")
  } else {
    cov_matrix <- tryCatch(solve(hessian), error = function(e) NULL)
  }
  
  # standard error for alpha_est
  st_alpha <- rep(NA, length(alpha_est))
  for (i in 1:length(alpha_est)){
    var_alpha <- cov_matrix[i,i]
    st_alpha[i] <- sqrt(var_alpha)
  }
  
  # confidence intervals for alpha_est coordinates
  z <- qnorm(0.975)
  lower_alpha <- rep(NA, length(alpha_est))
  upper_alpha <- rep(NA, length(alpha_est))
  for (i in 1:length(alpha_est)){
    var_alpha <- cov_matrix[i,i]
    st_alpha <- sqrt(var_alpha)
    lower_alpha[i] <- alpha_est[i] - z*st_alpha
    upper_alpha[i] <- alpha_est[i] + z*st_alpha
  }
  
  # standard error for xi
  st_xi <- sum(as.numeric(N)^as.vector(X %*% st_alpha))
  
  # confidence intervals for xi
  ci_xi <- c(sum(as.numeric(N)^as.vector(X %*% lower_alpha)), sum(as.numeric(N)^as.vector(X %*% upper_alpha)))
  
  return(
    list(
      estimates = estimates,
      cov_matrix = cov_matrix,
      confint_alpha = data.frame(lower = lower_alpha, upper = upper_alpha),
      confint_xi = setNames(ci_xi, c('lower', 'upper')),
      st_error_alpha = st_alpha,
      st_error_xi = st_xi,
      optim_result = optimization)
  )
  
}


### LOG-LIKELIHOOD ###
log_lik_zhang_model_cov <- function(alpha, beta, phi, m, n ,N, X, Z){
  
  mu <- N^(as.vector(X %*% alpha)) * (n/N)^(as.vector(Z %*% beta))
  
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


### GRADIENT ###
grad_log_lik_zhang_model_cov <- function(alpha, beta, phi, m, n ,N, X, Z){
  
  mu <- N^(as.vector(X %*% alpha)) * (n/N)^(as.vector(Z %*% beta))
  
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
  
  d_phi <- 1/(2*phi) - m/(mu + phi) - log(mu + phi)- phi/(mu + phi) + log(m + phi) + (m + phi - 0.5)/(m + phi)
  
  grad_alpha <- as.vector(t(X) %*% d_alpha)
  grad_beta <- as.vector(t(Z) %*% d_beta)
  grad_phi <- sum(d_phi)
  
  return(c(grad_alpha, grad_beta, grad_phi))
  
}


### HESSIAN ###
hess_log_lik_zhang_model_cov <- function(alpha, beta, phi, m, n ,N, X, Z){
  
  mu <- N^(as.vector(X %*% alpha)) * (n/N)^(as.vector(Z %*% beta))
  
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
  H11 <- t(X) %*% (d_alpha2 * X)
  H12 <- t(X) %*% (d_alpha_beta * Z)
  H13 <- matrix(t(X) %*% d_alpha_phi, ncol=1)
  
  H22 <- t(Z) %*% (d_beta2 * Z)
  H23 <- matrix(t(Z) %*% d_beta_phi, ncol=1)
  
  H33 <- matrix(sum(d_phi2), ncol=1, nrow=1)
  
  # full hessian matrix
  top_row <- cbind(H11, H12, H13)
  middle_row <- cbind(H12, H22, H23)
  bottom_row <- cbind(H13, H23, H33)
  
  return(rbind(top_row, middle_row, bottom_row))
  
}

