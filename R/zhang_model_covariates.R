
zhang_model_cov <- function(m, n, N, start = 'glm', X = NULL, Z = NULL){

  if (is.null(X)==TRUE)  {X <- matrix(1, nrow = length(n), ncol = 1)}
  if (is.null(Z)==TRUE)  {Z <- matrix(1, nrow = length(n), ncol = 1)}

  p1 <- ncol(X)
  p2 <- ncol(Z)

  df <- data.frame(
    y = m,
    x1 = log(N),
    x2 = log(n/N)
  )

  glm_fit <- glm(y ~ x1 + x2 - 1, data = df, family = poisson(link = 'log'))
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
    -log_lik_zhang_model_cov(alpha, beta, phi, m, n, N, X, Z)    # negative because optim() minimizes
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

  alpha_est <- optimization$par[1:p1]
  beta_est <- optimization$par[(p1+1):(p1+p2)]
  phi_est <- optimization$par[p1+p2+1]
  xi_est <- sum(N^(X %*% alpha_est))     # target parameter estimator


  return(
    list(
      coefficients = c(alpha_est, beta_est, phi_est, xi_est),
      optim_result = optimization)
  )

}



