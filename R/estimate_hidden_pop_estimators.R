#' @importFrom sandwich vcovHC
#' @importFrom stats confint qnorm residuals setNames var vcov lm nls glm coef model.matrix
#' @importFrom utils capture.output

## the OLS model
ols_model <- function(m,
                      n,
                      N,
                      vcov = 'hessian',
                      countries){

  df <- data.frame(
    y = m,
    x1 = log(N),
    x2 = log(n/N)
  )

  ols_fit <- lm(log(y) ~ x1 + x2 - 1, data = df)

  alpha_est <- unname(coef(ols_fit)[1])
  beta_est <- unname(coef(ols_fit)[2])
  xi_est <- sum(N^alpha_est)

  # estimates by nationality
  irreg_estimates <- N^alpha_est
  if(!is.null(countries)){
    by_nat_split <- data.frame(country = countries,
                               irreg_estimate = irreg_estimates)
    by_nationality <- aggregate(irreg_estimate ~ country, data = by_nat_split, sum)
  } else {
    by_nationality <- NULL
  }

  coef <- c(alpha = alpha_est, beta = beta_est)

  # covariance matrix
  if (vcov == 'robust') {
    cov_matrix <- sandwich::vcovHC(ols_fit, type = "HC1")
  } else {
    cov_matrix <- vcov(ols_fit)
  }

  # standard errors for coefficients
  se_coef <- data.frame(name = c('alpha', 'beta'),
                        Std.error = sqrt(diag(cov_matrix)))

  # standard error for xi - delta method
  se_alpha <- se_coef$Std.error[se_coef$name == 'alpha']
  se_xi <- abs(sum(N^alpha_est * log(N))) * se_alpha

  # confidence intervals for coefficients
  ci_alpha <- confint(ols_fit)[1,]
  ci_beta <- confint(ols_fit)[2,]
  conf_int_coef <- data.frame(name = c('alpha', 'beta'),
                               Lower = c(ci_alpha[1], ci_beta[1]),
                               Upper = c(ci_alpha[2], ci_beta[2]))

  # confidence intervals for xi - M estimate
  conf_int_xi <- data.frame(Lower = sum(N^ci_alpha[1]),
                            Upper = sum(N^ci_alpha[2]))

  # AIC, BIC values
  aic <- AIC(ols_fit)
  bic <- BIC(ols_fit)

  # summary lm
  summary_stat <- list(resid_se = summary(ols_fit)$sigma,
                       df_resid = summary(ols_fit)$df[2],
                       r_squared = summary(ols_fit)$r.squared,
                       adj_r_squared = summary(ols_fit)$adj.r.squared,
                       f_stat = summary(ols_fit)$fstatistic,
                       aic = aic,
                       bic = bic)

  results <- list(method = 'ols',
                  coefficients = coef,
                  xi_est = xi_est,
                  se_xi = se_xi,
                  se_coef = se_coef,
                  vcov_method = vcov,
                  vcov = cov_matrix,
                  conf_int_xi = conf_int_xi,
                  conf_int_coef = conf_int_coef,
                  summary_stat = summary_stat,
                  residuals = ols_fit$residuals,
                  fitted = fitted(ols_fit),
                  resid_stand = rstandard(ols_fit),
                  cooks = cooks.distance(ols_fit),
                  leverage = hatvalues(ols_fit),
                  m = m,
                  n = n,
                  N = N,
                  countries = countries,
                  by_nationality = by_nationality)

  return(results)
}



## nls model
nls_model <- function(m,
                      n,
                      N,
                      vcov = 'hessian',
                      countries){

  df_nls <- data.frame(
    y = m,
    N = N,
    n = n
  )

  df_nls$x1 <- log(df_nls$N)
  df_nls$x2 <- log(df_nls$n/df_nls$N)
  ols_fit <- lm(log(y) ~ x1 + x2 - 1, data = df_nls)
  estim_ols <- list(alpha = coef(ols_fit)[1], beta = coef(ols_fit)[2])

  if (coef(ols_fit)[1] < 0 || coef(ols_fit)[2] < 0) {
    warning("OLS produced negative starting values. They were changed to proceed with nls.")
  }
  estim_ols$alpha <- min(max(estim_ols$alpha, 0.01), 1)
  estim_ols$beta <- min(max(estim_ols$beta, 0.01), 1)

  nls_fit <- nls(y ~ N^alpha * (n/N)^beta,
                 data = df_nls,
                 start = estim_ols)

  estim_nls <- coef(nls_fit)
  alpha_est <- estim_nls['alpha']
  beta_est <- estim_nls['beta']
  xi_est <- sum(N^alpha_est)

  # estimates by nationality
  irreg_estimates <- N^alpha_est
  if(!is.null(countries)){
    by_nat_split <- data.frame(country = countries,
                               irreg_estimate = irreg_estimates)
    by_nationality <- aggregate(irreg_estimate ~ country, data = by_nat_split, sum)
  } else {
    by_nationality <- NULL
  }

  coef <- c(alpha_est, beta_est)

  # covariance matrix
  if (vcov == 'robust') {
    cov_matrix <- robust_vcov_nls_hc1(nls_fit)
  } else {
    cov_matrix <- vcov(nls_fit)
  }

  # standard errors for coefficients
  se_coef <- data.frame(name = c('alpha', 'beta'),
                        Std.error = sqrt(diag(cov_matrix)))

  # standard error for xi - delta method
  se_alpha <- se_coef$Std.error[se_coef$name == 'alpha']
  se_xi <- abs(sum(N^alpha_est * log(N))) * se_alpha

  # confidence intervals for coefficients
  ci_alpha <- confint(nls_fit)[1,]
  ci_beta <- confint(nls_fit)[2,]
  conf_int_coef <- data.frame(name = c('alpha', 'beta'),
                              Lower = c(ci_alpha[1], ci_beta[1]),
                              Upper = c(ci_alpha[2], ci_beta[2]))

  # confidence intervals for xi - M estimate
  conf_int_xi <- data.frame(Lower = sum(N^ci_alpha[1]),
                            Upper = sum(N^ci_alpha[2]))


  # AIC, BIC values
  aic <- AIC(nls_fit)
  bic <- BIC(nls_fit)

  # summary nls
  summary_stat <- list(resid_se = summary(nls_fit)$sigma,
                       df = summary(nls_fit)$df[2],
                       iter = nls_fit$convInfo$finIter,
                       convergence = nls_fit$convInfo$isConv,
                       aic = aic,
                       bic = bic)

  results <- list(method = 'nls',
                  coefficients = coef,
                  xi_est = xi_est,
                  se_coef = se_coef,
                  se_xi = se_xi,
                  vcov_method = vcov,
                  vcov = cov_matrix,
                  conf_int_xi = conf_int_xi,
                  conf_int_coef = conf_int_coef,
                  summary_stat = summary_stat,
                  residuals = resid(nls_fit),
                  fitted = fitted(nls_fit),
                  m = m,
                  n = n,
                  N = N,
                  countries = countries,
                  by_nationality = by_nationality)

  return(results)

}


 ## Poisson regression
glm_model <- function(m,
                      n,
                      N,
                      vcov = 'hessian',
                      countries){

  df <- data.frame(
    y = m,
    x1 = log(N),
    x2 = log(n/N)
  )

  glm_fit <- glm(y ~ x1 + x2 - 1, data = df, family = poisson(link = 'log'))

  alpha_est <- unname(coef(glm_fit)[1])
  beta_est <- unname(coef(glm_fit)[2])
  xi_est <- sum(N^alpha_est)

  # estimates by nationality
  irreg_estimates <- N^alpha_est
  if(!is.null(countries)){
    by_nat_split <- data.frame(country = countries,
                               irreg_estimate = irreg_estimates)
    by_nationality <- aggregate(irreg_estimate ~ country, data = by_nat_split, sum)
  } else {
    by_nationality <- NULL
  }

  coef <- c(alpha = alpha_est, beta = beta_est)

  # covariance matrix
  if (vcov == 'robust') {
    cov_matrix <- sandwich::vcovHC(glm_fit, type = "HC1")
  } else {
    cov_matrix <- vcov(glm_fit)
  }

  # standard errors for coefficients
  se_coef <- data.frame(name = c('alpha', 'beta'),
                        Std.error = sqrt(diag(cov_matrix)))

  # standard error for xi - delta method
  se_alpha <- se_coef$Std.error[se_coef$name == 'alpha']
  se_xi <- abs(sum(N^alpha_est * log(N))) * se_alpha

  # confidence intervals for coefficients
  ci_alpha <- confint(glm_fit)[1,]
  ci_beta <- confint(glm_fit)[2,]
  conf_int_coef <- data.frame(name = c('alpha', 'beta'),
                               Lower = c(ci_alpha[1], ci_beta[1]),
                               Upper = c(ci_alpha[2], ci_beta[2]))

  # confidence intervals for xi - M estimate
  conf_int_xi <- data.frame(Lower = sum(N^ci_alpha[1]),
                            Upper = sum(N^ci_alpha[2]))

  # AIC, BIC values
  aic <- AIC(glm_fit)
  bic <- BIC(glm_fit)

  # summary glm
  summary_stat <- list(resid_deviance = summary(glm_fit)$deviance,
                       df_resid = summary(glm_fit)$df.residual,
                       null_deviance = summary(glm_fit)$null.deviance,
                       df_null = summary(glm_fit)$df.null,
                       iter = summary(glm_fit)$iter,
                       aic = aic,
                       bic = bic)

  results <- list(method = 'glm - Poisson',
                  coefficients = coef,
                  xi_est = xi_est,
                  se_coef = se_coef,
                  se_xi = se_xi,
                  vcov_method = vcov,
                  vcov = cov_matrix,
                  conf_int_xi = conf_int_xi,
                  conf_int_coef = conf_int_coef,
                  summary_stat = summary_stat,
                  residuals = residuals(glm_fit, type = 'deviance'),
                  fitted = fitted(glm_fit),
                  resid_stand = rstandard(glm_fit, type = 'deviance'),
                  cooks = cooks.distance(glm_fit),
                  leverage = hatvalues(glm_fit),
                  m = m,
                  n = n,
                  N = N,
                  countries = countries,
                  by_nationality = by_nationality)

  return(results)

}



## mle model

mle_estim <- function(m,
                      n,
                      N,
                      X = NULL,
                      Z = NULL,
                      vcov = 'hessian',
                      countries,
                      df_cov,
                      family){

  log_lik_cov <- function(alpha, beta, phi=NULL, m, n, N, X, Z, family){

    mu <- N^(as.vector(X %*% alpha)) * (n/N)^(as.vector(Z %*% beta))

    # check point
    if (any(!is.finite(mu)) || any(mu <= 0)) return(Inf)   # obvious why, log(mu) if mu <= 0  -> error
    if (family == 'nb' & (!is.finite(phi) || phi <= 0)) return(Inf)           # if phi = 0, then lgamma(phi)= Inf

    val <- switch(family, 'poisson' = m * log(mu) - mu - lfactorial(m),
                  'nb' = m * log(mu) + phi * log(phi) -
                    (m + phi) * log(mu + phi) +
                    lgamma(m + phi) - lfactorial(m) - lgamma(phi))

    # check if there are any inf
    if (any(!is.finite(val))) return(Inf)

    return(sum(val))
  }

  grad_log_lik_cov <- function(alpha, beta, phi = NULL, m, n ,N, X, Z, family){

    mu <- N^(as.vector(X %*% alpha)) * (n/N)^(as.vector(Z %*% beta))

    logN <- log(N)
    logRatio <- log(n/N)

    # check point
    if (any(mu <= 0) || any(!is.finite(mu)) ||   # obvious
        (family == 'nb' & (phi <= 0 || !is.finite(phi)))) {           # obvious
      return(rep(NA, ncol(X) + ncol(Z) + 1))     # to have correct gradient dimension and still try to go into optim()
    }

    d_mu <- switch(family, 'poisson' = m/mu - 1,
                   'nb' = m/mu - (m + phi)/(phi + mu))
    d_alpha <- d_mu * mu * logN
    d_beta <- d_mu * mu * logRatio

    if (family == 'nb'){
      log_mu_phi <- log(mu + phi)
      log_m_phi  <- log(m + phi)
    }

    # check point
    if(family == 'nb'){
      if (any(!is.finite(log_mu_phi)) ||        # to avoid problems in sum d_phi
          any(!is.finite(log_m_phi))) {
        return(rep(NA, ncol(X) + ncol(Z) + 1))  # to have correct gradient dimension and still try to go into optim()
      }}

    if (family =='nb') {
      d_phi <- 1/(2*phi) - m/(mu + phi) - log(mu + phi)- phi/(mu + phi) + log(m + phi) + (m + phi - 0.5)/(m + phi)
    }

    grad_alpha <- as.vector(t(X) %*% d_alpha)
    grad_beta <- as.vector(t(Z) %*% d_beta)
    if (family == 'nb') {grad_phi <- sum(d_phi)}

    grads <- switch(family, 'poisson' = c(grad_alpha, grad_beta),
                    'nb' = c(grad_alpha, grad_beta, grad_phi))

    return(grads)

  }

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
  start_par <- switch(family, 'poisson' = c(alpha_start, beta_start),
                      'nb' = c(alpha_start, beta_start, 1/var(glm_fit$residuals)))

  # negative log-likelihood, gradient and hessian

  neg_log_lik <- function(param){     # we will take start_par as param in optim function
    alpha <- param[1:p1]
    beta <- param[(p1+1):(p1+p2)]
    phi <- switch(family, 'poisson' = NULL,
                  'nb' = param[p1+p2+1])
    -log_lik_cov(alpha, beta, phi, m, n, N, X, Z, family)
  }

  neg_grad_log_lik <- function(param){
    alpha <- param[1:p1]
    beta <- param[(p1+1):(p1+p2)]
    phi <- switch(family, 'poisson' = NULL,
                  'nb' = param[p1+p2+1])
    -grad_log_lik_cov(alpha, beta, phi, m, n, N, X, Z, family)
  }

  neg_hess_log_lik <- function(param){
    alpha <- param[1:p1]
    beta <- param[(p1+1):(p1+p2)]
    phi <- switch(family, 'poisson' = NULL,
                  'nb' = param[p1+p2+1])
    -hess_log_lik_cov(alpha, beta, phi, m, n, N, X, Z, family)
  }

  # optimization
  optimization <- optim(
    par = start_par,
    fn = neg_log_lik,
    gr = neg_grad_log_lik,
    method = 'BFGS',
    hessian = TRUE
  )

  alpha_est <- unname(optimization$par[1:p1])
  beta_est <- unname(optimization$par[(p1+1):(p1+p2)])
  if (family == 'nb'){phi_est <- unname(optimization$par[p1+p2+1])}
  xi_est <- sum(as.numeric(N)^as.vector(X %*% alpha_est))     # target parameter estimator

  # estimates by nationality
  irreg_estimates <- as.numeric(N)^as.vector(X %*% alpha_est)
  if(!is.null(countries)){
    by_nat_split <- data.frame(country = countries,
                               irreg_estimate = irreg_estimates)
    by_nationality <- aggregate(irreg_estimate ~ country, data = by_nat_split, sum)
  } else {
    by_nationality <- NULL
  }


  # estimates by covariates
  if(!is.null(df_cov)){
    covariate_vars <- colnames(df_cov)
    df_cov$irreg_estimate <- as.vector(irreg_estimates)
    formula <- as.formula(paste0('irreg_estimate ~', paste(covariate_vars, collapse = '+')))
    by_covariates <- aggregate(formula, data = df_cov, sum)
  } else {
    by_covariates <- NULL
  }

  names(alpha_est) <- paste0("alpha", seq_along(alpha_est))
  names(beta_est) <- paste0("beta", seq_along(beta_est))

  coef <- switch(family, 'poisson' = c(setNames(alpha_est, paste0('alpha', seq_along(alpha_est))),
                                       setNames(beta_est, paste0('beta', seq_along(beta_est)))),
                 'nb' = c(setNames(alpha_est, paste0('alpha', seq_along(alpha_est))),
                          setNames(beta_est, paste0('beta', seq_along(beta_est))),
                          phi = phi_est))

  hessian <- optimization$hessian

  # covariance matrix

  # mle robust covariance matrix
  robust_mle <- function(opt, m, n, N, X, Z){

    p1 <- ncol(X)
    p2 <- ncol(Z)

    grad_i <- function(param){

      phi_rob <- switch(family, 'poisson' = NULL,
                        'nb' = param[p1+p2+1])

      result <- lapply(1:length(m), function(i){
        grad_log_lik_cov(
          alpha = param[1 : p1],
          beta = param[(p1+1) : (p1+p2)],
          phi = phi_rob,
          m = m[i],
          n = n[i],
          N = N[i],
          matrix(X[i, ], nrow = 1),
          matrix(Z[i, ], nrow = 1),
          family
        )})

      return(result)
    }

    grads <- do.call(cbind, grad_i(opt$par))

    H <- opt$hessian
    inv_H <- MASS::ginv(H)

    I <- grads %*% t(grads)

    return(inv_H %*% I %*% inv_H)
  }

  if (vcov == 'robust') {
    #cov_matrix <- sandwich::vcovHC(glm_fit, type = 'HC1')
    cov_matrix <- robust_mle(optimization, m,n,N, X,Z)
  } else {
    cov_matrix <- tryCatch(solve(hessian), error = function(e) NULL)
  }

  # standard error and confidence intervals for alpha_est coordinates
  se_alpha <- rep(NA, length(alpha_est))
  z <- qnorm(0.975)
  lower_alpha <- rep(NA, length(alpha_est))
  upper_alpha <- rep(NA, length(alpha_est))
  for (i in 1:length(alpha_est)){
    var_alpha <- cov_matrix[i,i]
    se_alpha[i] <- sqrt(var_alpha)
    lower_alpha[i] <- alpha_est[i] - z*se_alpha[i]
    upper_alpha[i] <- alpha_est[i] + z*se_alpha[i]
  }

  ci_alpha <- data.frame(name = names(alpha_est),
                         Lower = lower_alpha,
                         Upper = upper_alpha)

  # for beta
  se_beta <- rep(NA, length(beta_est))
  z <- qnorm(0.975)
  lower_beta <- rep(NA, length(beta_est))
  upper_beta <- rep(NA, length(beta_est))
  for (j in 1:length(beta_est)){
    i <- p1 + j
    var_beta <- cov_matrix[i,i]
    se_beta[j] <- sqrt(var_beta)
    lower_beta[j] <- beta_est[j] - z*se_beta[j]
    upper_beta[j] <- beta_est[j] + z*se_beta[j]
  }

  ci_beta <- data.frame(name = names(beta_est),
                        Lower = lower_beta,
                        Upper = upper_beta)

  # for phi
  if(family == 'nb'){
    se_phi <- cov_matrix[p1+p2+1, p1+p2+1]
    lower_phi <- phi_est - z*se_phi
    upper_phi <- phi_est + z*se_phi
    ci_phi <- data.frame(name = 'phi',
                         Lower = lower_phi,
                         Upper = upper_phi)
  }

  # confidence intervals for coefficients
  conf_int_coef <- switch(family, 'poisson' = data.frame(name = c(names(alpha_est), names(beta_est)),
                                                         Lower = c(ci_alpha$Lower, ci_beta$Lower),
                                                         Upper = c(ci_alpha$Upper, ci_beta$Upper)),
                          'nb' = data.frame(name = c(names(alpha_est), names(beta_est),'phi'),
                                            Lower = c(ci_alpha$Lower, ci_beta$Lower, ci_phi$Lower),
                                            Upper = c(ci_alpha$Upper, ci_beta$Upper, ci_phi$Upper)))

  # standard errors for coefficients
  se_coef <- switch(family, 'poisson' = data.frame(name = c(paste0('alpha', seq_along(alpha_est)), paste0('beta', seq_along(beta_est))),
                                                   Std.error = c(se_alpha, se_beta)),
                    'nb' = data.frame(name = c(paste0('alpha', seq_along(alpha_est)), paste0('beta', seq_along(beta_est)), 'phi'),
                                      Std.error = c(se_alpha, se_beta, se_phi)))

  # standard error for xi
  N_Xalpha <- as.numeric(N)^as.vector(X %*% alpha_est)
  grad_g_alpha <- t(X) %*% (log(N) * N_Xalpha)
  se_xi <- se_alpha %*% abs(grad_g_alpha)

  # confidence intervals for xi
  conf_int_xi <- data.frame(Lower = sum(as.numeric(N)^as.vector(X %*% lower_alpha)),
                            Upper = sum(as.numeric(N)^as.vector(X %*% upper_alpha)))

  # AIC, BIC values
  LL <- switch(family, 'poisson' = log_lik_cov(alpha_est, beta_est, NULL, m, n, N, X,Z, family),
               'nb' = log_lik_cov(alpha_est, beta_est, phi_est, m, n, N, X,Z, family))
  k <- length(alpha_est) + length(beta_est) + 1
  C <- length(m)
  aic <- -2 * LL + 2 * k
  bic <- -2 * LL + k * log(C)

  # residuals
  residuals <- m - N^(X %*% alpha_est) * (n/N)^(Z %*% beta_est)  # should it be different type ???

  # fitted
  fitted <- N^(X %*% alpha_est) * (n/N)^(Z %*% beta_est)

  # summary mle
  summary_stat <- list(convergence = optimization$convergence,
                       iter = optimization$counts,
                       aic = aic,
                       bic = bic)

  results <- list(method = paste0('mle - ', family),
                  coefficients = coef,
                  xi_est = xi_est,
                  se_coef = se_coef,
                  se_xi = se_xi,
                  vcov_method = vcov,
                  vcov = cov_matrix,
                  conf_int_xi = conf_int_xi,
                  conf_int_coef = conf_int_coef,
                  summary_stat = summary_stat,
                  residuals = residuals,
                  fitted = fitted,
                  m = m,
                  n = n,
                  N = N,
                  countries = countries,
                  by_nationality = by_nationality,
                  by_covariates = by_covariates)

  return(results)

}

# nls robust
robust_vcov_nls_hc1 <- function(nls_model){

  # residuals
  res <- residuals(nls_model)

  # gradient
  X <- nls_model$m$gradient()
  n <- nrow(X)
  k <- ncol(X)

  meat <- t(X) %*% diag(res^2) %*% X

  bread <- MASS::ginv(t(X) %*% X)

  cov_matrix <- (n/(n-k)) * bread %*% meat %*% bread

  return(cov_matrix)
}


