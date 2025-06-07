
zhang_model_cov <- function(m, n, N, X = NULL, Z = NULL, start = "glm"){

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

  alpha_est <- unname(optimization$par[1:p1])
  beta_est <- unname(optimization$par[(p1+1):(p1+p2)])
  phi_est <- unname(optimization$par[p1+p2+1])
  xi_est <- sum(N^as.vector(X %*% alpha_est))     # target parameter estimator

  estimates <- list(xi = xi_est, alpha = alpha_est, beta = beta_est, phi = phi_est)

  # confidence intervals for alpha_est coordinates
  hessian <- optimization$hessian
  cov_matrix <- tryCatch(solve(hessian), error = function(e) NA)  # is it possible - check
  z <- qnorm(0.975)
  lower_alpha <- rep(NA, length(alpha_est))
  upper_alpha <- rep(NA, length(alpha_est))
  for (i in 1:length(alpha_est)){
    var_alpha <- cov_matrix[i,i]
    s_alpha <- sqrt(var_alpha)
    lower_alpha[i] <- alpha_est[i] - z*s_alpha
    upper_alpha[i] <- alpha_est[i] + z*s_alpha
  }

  # confidence intervals for xi
  ci_xi <- c(sum(N^as.vector(X %*% lower_alpha)), sum(N^as.vector(X %*% upper_alpha)))

  return(
    list(
      estimates = estimates,
      confint_alpha = data.frame(lower = lower_alpha, upper = upper_alpha),
      confint_xi = setNames(ci_xi, c('lower', 'upper')),
      optim_result = optimization)
  )

}

