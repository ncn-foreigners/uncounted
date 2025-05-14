
zhang_model_cov <- function(m, n, N, start = 'glm', X = NULL, Z = NULL){
                            
  if (is.null(X)==TRUE)  {X <- matrix(1, nrow = length(n), ncol = 1)} 
  if (is.null(Z)==TRUE)  {Z <- matrix(1, nrow = length(n), ncol = 1)}
  
  df <- data.frame(
    y = m,
    x1 = log(N),
    x2 = log(n/N)
  )

  glm_fit <- glm(y ~ x1 + x2 - 1, data = df, family = poisson(link = 'log'))
  starting <- coef(glm_fit)
  if (any(!is.finite(starting))) stop('GLM - infinite values')
  
  # log-likelihood, gradient and hessian
  
  log_lik <- function(param){
    alpha <- param[1]
    beta <- param[2]
    phi <- param[3]               
    -log_lik_zhang_model_cov(alpha, beta, phi, m, n, N, X, Z)    # negative because optim() minimizes
  }
  
  grad_log_lik <- function(param){
    alpha <- param[1]
    beta <- param[2]
    phi <- param[3]                
    -grad_log_lik_zhang_model_cov(alpha, beta, phi, m, n, N, X, Z)
  }  
  
  hess_log_lik <- function(param){
    alpha <- param[1]
    beta <- param[2]
    phi <- param[3]                
    -hess_log_lik_zhang_model_cov(alpha, beta, phi, m, n, N, X, Z)
  }
  
  # optimization
  
  start_par <- c(starting[1], starting[2], 0.01)
  
  obliczenia <- optim(
    par = start_par,    
    fn = log_lik,     
    gr = grad_log_lik,   
    method = 'BFGS'
  )
  
  alpha_est <- obliczenia$par[1]
  beta_est <- obliczenia$par[2]
  phi_est <- obliczenia$par[3]
  xi_est <- sum(N^alpha_est)     # target parameter estimator
  
  
  return(
    list(
      coefficients = c(alpha_est, beta_est, phi_est, xi_est),
      optim_result = obliczenia)
  )
  
}



