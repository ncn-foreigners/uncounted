
zhang_model <- function(m, n, N){
  
  df <- data.frame(
    y = m,
    x1 = log(N),
    x2 = log(n/N)
  )
  
  glm_fit <- glm(y ~ x1 + x2 - 1, data = df, family = poisson(link = 'log'))
  starting <- coef(glm_fit)
  if (any(!is.finite(starting))) stop('GLM - infinite values')   # just in case


  # log-likelihood, gradient and hessian
  
  log_lik <- function(param){
    alpha <- param[1]
    beta <- param[2]
    phi <- param[3] 
    -log_lik_zhang_model(alpha, beta, phi, m, n, N)    # negative because optim() minimizes
  }
  
  grad_log_lik <- function(param){
    alpha <- param[1]
    beta <- param[2]
    phi <- param[3]         
    -grad_log_lik_zhang_model(alpha, beta, phi, m, n, N)
  }  
  
  hess_log_lik <- function(param){
    alpha <- param[1]
    beta <- param[2]
    phi <- param[3]         
    -hess_log_lik_zhang_model(alpha, beta, phi, m, n, N)
  }
  
  # optimization
  
  start_par <- c(starting[1], starting[2], 0.01)

  optimization <- optim(
    par = start_par,     # starting points
    fn = log_lik,        # a function to be minimized / maximized
    gr = grad_log_lik,   
    method = 'BFGS',  
    hessian = TRUE  
  )
  
  alpha_est <- optimization$par[1]
  beta_est <- optimization$par[2]
  phi_est <- optimization$par[3]
  xi_est <- sum(N^alpha_est)     # target parameter estimator
  

    return(
      list(
        coefficients = c(alpha_est, beta_est, phi_est, xi_est),
        optim_result = optimization)
    )
    
}



