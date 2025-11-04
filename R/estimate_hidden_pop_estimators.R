#' @importFrom sandwich vcovHC
#' @importFrom stats confint qnorm residuals setNames var vcov lm nls glm coef model.matrix
#' @importFrom utils capture.output
#' @importFrom boot boot
#' @importFrom fwb fwb
#' @importFrom fANCOVA wild.boot

## the OLS model
ols_model <- function(m,
                      n,
                      N,
                      vcov = 'hessian',
                      countries){

  df <- data.frame(y = m,
                   x1 = log(N),
                   x2 = log(n/N))

  ols_fit_fun <- function(data = df, weights = NULL, indices = NULL){

    if (!is.null(weights)) {
      weights <- weights
      ols_fit <- lm(log(y) ~ x1 + x2 - 1, data = data, weights = weights)
    } else if (!is.null(indices)) {
        data <- data[indices, ]
        ols_fit <- lm(log(y) ~ x1 + x2 - 1, data = data)
    } else {
      ols_fit <- lm(log(y) ~ x1 + x2 - 1, data = data)
    }

    alpha_est <- unname(coef(ols_fit)[1])
    beta_est <- unname(coef(ols_fit)[2])
    xi_est <- sum(N^alpha_est)

    return(c(alpha_est, beta_est, xi_est))
  }

  ols_fit_nonpar <- function(data, indices = indices) {
    ols_fit_fun(data = data, indices = indices)
  }

  ols_fit_fwb <- function(data, weights = weights) {
    ols_fit_fun(data = data, weights = weights)
  }

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
  } else if (vcov == 'hessian') {
    cov_matrix <- vcov(ols_fit)
  } else if (vcov == 'nonparametric') {
    nonpar_results <- boot(data = df, statistic = ols_fit_nonpar, R = 1000)
    cov_matrix <- cov(nonpar_results$t[,1:2]) # dla alpha i beta
    # se_xi <- sd(nonpar_results[,3])
  } else if (vcov == 'wild') {
    fitted_vals <- fitted(ols_fit)
    R <- 1000
    wild_residuals <- wild.boot(ols_fit$residuals, R)
    alpha_boot <- numeric(R)
    beta_boot <- numeric(R)
    xi_boot <- numeric(R)

    for (i in 1:R) {
      y_boot <- fitted_vals + wild_residuals[, i]
      boot_fit <- lm(y_boot ~ x1 + x2 - 1, data = df)
      coef_wild <- coef(boot_fit)
      alpha_boot[i] <- unname(coef_wild[1])
      beta_boot[i] <- unname(coef_wild[2])
      xi_boot[i] <- sum(N^alpha_boot[i])
    }

    boot_matrix <- cbind(alpha_boot, beta_boot)
    cov_matrix <- cov(boot_matrix)
    # se_xi <- sd(xi_boot)
  } else if (vcov == 'fwb') {
    fwb_results <- fwb(data = df, statistic = ols_fit_fwb, R = 1000, verbose = FALSE)
    cov_matrix <- cov(fwb_results[,1:2])
    # se_xi <- sd(fwb_results[,3])
  }

  # standard errors for coefficients
  se_coef <- data.frame(name = c('alpha', 'beta'),
                        Std.error = sqrt(diag(cov_matrix)))

  # standard error for xi - delta method
  se_alpha <- se_coef$Std.error[se_coef$name == 'alpha']
  se_xi <- abs(sum(N^alpha_est * log(N))) * se_alpha   #  TUTAJ POMYSLEC JAK Z BOOTSTRAPEM

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
  # } else if (vcov == 'nonparametric') {
  #
  # } else if (vcov == 'wild') {
  #
  # }

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
                      family,
                      bias_corr){

  log_lik_cov <- function(alpha, beta, phi=NULL, m, n, N, X, Z, family, weights = NULL){

    if(is.null(weights)) weights <- rep(1, length(m))

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

    return(sum(weights * val))
  }

  grad_log_lik_cov <- function(alpha, beta, phi = NULL, m, n ,N, X, Z, family, weights = NULL){

    if(is.null(weights)) weights <- rep(1, length(m))

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

    # check point
    if(family == 'nb'){
      if (any(!is.finite(log(mu + phi))) ||        # to avoid problems in sum d_phi
          any(!is.finite(log(m + phi)))) {
        return(rep(NA, ncol(X) + ncol(Z) + 1))  # to have correct gradient dimension and still try to go into optim()
      }}

    if (family =='nb') {
      d_phi <- 1/(2*phi) - m/(mu + phi) - log(mu + phi)- phi/(mu + phi) + log(m + phi) + (m + phi - 0.5)/(m + phi)
    }

    grad_alpha <- as.vector(t(X) %*% (weights * d_alpha))
    grad_beta <- as.vector(t(Z) %*% (weights * d_beta))
    if (family == 'nb') grad_phi <- sum(weights * d_phi)

    grads <- switch(family, 'poisson' = c(grad_alpha, grad_beta),
                    'nb' = c(grad_alpha, grad_beta, grad_phi))
    return(grads)

  }

  if (is.null(X)) X <- matrix(1, nrow = length(n), ncol = 1)
  if (is.null(Z)) Z <- matrix(1, nrow = length(n), ncol = 1)

  p1 <- ncol(X)
  p2 <- ncol(Z)

  df <- data.frame(
    y = m,
    x1 = log(N),
    x2 = log(n/N)
  )

  glm_fit <- glm(y ~ x1 + x2 - 1, data = df, family = poisson(link = 'log'))

  # starting points for optimization
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
  if (family == 'nb') phi_est <- unname(optimization$par[p1+p2+1])
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

  # covariance matrix calculation
  # hessian
  hessian <- optimization$hessian

  # function to use in bootstrap
  mle_fit_fun <- function(data, weights = NULL, indices = NULL){

    if(!is.null(indices)) data <- data[indices,]
    if(!is.null(weights)) weights <- weights

    m <- data$m
    n <- data$n
    N <- data$N
    X <- as.matrix(data[, grepl('^X', names(data)), drop = FALSE])
    Z <- as.matrix(data[, grepl('^Z', names(data)), drop = FALSE])
    if (ncol(X) == 0) X <- matrix(1, nrow = length(m), ncol = 1)
    if (ncol(Z) == 0) Z <- matrix(1, nrow = length(m), ncol = 1)
    p1 <- ncol(X)
    p2 <- ncol(Z)

    # starting values
    df_boot <- data.frame(y = m,
                          x1 = log(N),
                          x2 = log(n/N))

    glm_fit_boot <- glm(y ~ x1 + x2 - 1, data = df_boot, family = poisson(link = 'log'))

    starting_boot <- coef(glm_fit_boot)
    if (any(!is.finite(starting_boot))) stop('GLM - infinite values')

    alpha_start_boot <- rep(starting_boot[1], p1)
    beta_start_boot <- rep(starting_boot[2], p2)
    start_par_boot <- switch(family, 'poisson' = c(alpha_start_boot, beta_start_boot),
                        'nb' = c(alpha_start_boot, beta_start_boot, 1/var(glm_fit_boot$residuals)))

    optimization <- optim(par = start_par_boot,
                          fn = function(par){
                            alpha <- par[1:p1]
                            beta <- par[(p1 + 1):(p1 + p2)]
                            phi <- if (family == 'nb') {par[p1 + p2 + 1]} else NULL
                            -log_lik_cov(alpha, beta, phi, m, n, N, X, Z, family, weights)
                          },
                          gr = function(par) {
                            alpha <- par[1:p1]
                            beta <- par[(p1 + 1):(p1 + p2)]
                            phi <- if (family == 'nb') par[p1 + p2 + 1] else NULL
                            -grad_log_lik_cov(alpha, beta, phi, m, n, N, X, Z, family, weights)
                          },
                          method = 'BFGS',
                          hessian = FALSE)

    alpha_est <- unname(optimization$par[1:p1])
    beta_est <- unname(optimization$par[(p1+1):(p1+p2)])
    if (family == 'nb'){phi_est <- unname(optimization$par[p1+p2+1])}
    xi_est <- sum(as.numeric(N)^as.vector(X %*% alpha_est))

    coef_boot <- switch(family, 'poisson' = c(alpha_est, beta_est, xi_est),
                        'nb' = c(alpha_est, beta_est, phi_est, xi_est))

    return(coef_boot)
  }

  # support functions to have indices/weights as the second argument
  mle_fit_nonpar <- function(data, indices) {
    mle_fit_fun(data = data, indices = indices)
  }

  mle_fit_fwb <- function(data, weights) {
    mle_fit_fun(data = data, weights = weights)
  }

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
    cov_matrix <- robust_mle(optimization, m, n, N, X, Z)
  } else if (vcov == 'hessian'){
    cov_matrix <- tryCatch(solve(hessian), error = function(e) NULL)
  } else if (vcov == 'nonparametric') {

    df_boot <- cbind(data.frame(m = m, n = n, N = N), as.data.frame(X), as.data.frame(Z))
    names(df_boot)[4:(3 + p1)] <- paste0('X', 1:p1)
    names(df_boot)[(4 + p1):(3 + p1 + p2)] <- paste0('Z', 1:p2)

    nonpar_results <- boot(data = df_boot, statistic = mle_fit_nonpar, R = 1000)

    if (family == 'nb') {
      cov_matrix <- cov(nonpar_results$t[, 1:(p1 + p2 + 1)])
    } else if (family == 'poisson'){
      cov_matrix <- cov(nonpar_results$t[, 1:(p1 + p2)])
    }

  } else if (vcov == 'wild') {
    residuals <- m - N^(X %*% alpha_est) * (n/N)^(Z %*% beta_est)
    fitted <- N^(X %*% alpha_est) * (n/N)^(Z %*% beta_est)
    R <- 1000
    wild_residuals <- wild.boot(residuals, R)
    alpha_boot <- matrix(NA, nrow = R, ncol = p1)
    beta_boot  <- matrix(NA, nrow = R, ncol = p2)
    if(family == 'nb') {phi_boot <- numeric(R)}
    xi_boot <- numeric(R)

    for (i in 1:R) {
      y_boot <- round(pmax(fitted + wild_residuals[, i], 0))
      data_boot <- data.frame(m = y_boot, n = n, N = N)
      data_boot <- cbind(data_boot, as.data.frame(X)); names(data_boot)[4:(3 + p1)] <- paste0('X', 1:p1)
      data_boot <- cbind(data_boot, as.data.frame(Z)); names(data_boot)[(4 + p1):(3 + p1 + p2)] <- paste0('Z', 1:p2)

      # est <- mle_fit_fun(data_boot)
      # alpha_boot[i, ] <- est[1:p1]
      # beta_boot[i, ] <- est[(p1 + 1):(p1 + p2)]
      # if (family == 'nb') {
      #   phi_boot[i] <- est[p1 + p2 + 1]
      #   xi_boot[i] <- est[p1 + p2 + 2]
      # } else if (family =='poisson') {
      #   xi_boot[i] <- est[p1 + p2 + 1]
      # }
      tryCatch({est <- mle_fit_fun(data_boot)
      alpha_boot[i, ] <- est[1:p1]
      beta_boot[i, ] <- est[(p1 + 1):(p1 + p2)]
      if (family == 'nb') {
        phi_boot[i] <- est[p1 + p2 + 1]
        xi_boot[i] <- est[p1 + p2 + 2]
      } else if (family =='poisson') {
        xi_boot[i] <- est[p1 + p2 + 1]
      }
      }, error = function(e) {})
    }

    if (family == 'nb') {
      cov_matrix <- cov(cbind(alpha_boot, beta_boot, phi_boot))
    } else if (family == 'poisson') {
      cov_matrix <- cov(cbind(alpha_boot, beta_boot))
    }

  } else if (vcov == 'fwb') {

    df_boot <- cbind(data.frame(m = m, n = n, N = N), as.data.frame(X), as.data.frame(Z))
    names(df_boot)[4:(3 + p1)] <- paste0('X', 1:p1)
    names(df_boot)[(4 + p1):(3 + p1 + p2)] <- paste0('Z', 1:p2)
    fwb_results <- fwb(data = df_boot, statistic = mle_fit_fwb, R = 1000, verbose = FALSE)
    if (family == 'nb') {
      cov_matrix <- cov(fwb_results$t[, 1:(p1 + p2 + 1), drop = FALSE])
    } else if (family == 'poisson'){
      cov_matrix <- cov(fwb_results$t[, 1:(p1 + p2), drop = FALSE])
    }

  }

  # standard error and confidence intervals for alpha_est coordinates
  if (vcov %in% c('hessian', 'robust')){
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
  } else if (vcov == 'nonparametric') {
    se_alpha <- apply(nonpar_results$t[,1:p1, drop = FALSE], 2, sd, na.rm = TRUE)
    lower_alpha <- apply(nonpar_results$t[,1:p1, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE)
    upper_alpha <- apply(nonpar_results$t[,1:p1, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)
  } else if (vcov == 'wild') {
    se_alpha <- apply(alpha_boot, 2, sd, na.rm = TRUE)
    lower_alpha <- apply(alpha_boot, 2, quantile, probs = 0.025, na.rm = TRUE)
    upper_alpha <- apply(alpha_boot, 2, quantile, probs = 0.975, na.rm = TRUE)
  } else if (vcov == 'fwb') {
    se_alpha <- apply(fwb_results$t[,1:p1, drop = FALSE], 2, sd, na.rm = TRUE)
    lower_alpha <- apply(fwb_results$t[,1:p1, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE)
    upper_alpha <- apply(fwb_results$t[,1:p1, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)
  }

  ci_alpha <- data.frame(name = names(alpha_est),
                         Lower = lower_alpha,
                         Upper = upper_alpha)

  # for beta
  if (vcov %in% c('hessian', 'robust')){
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
  } else if (vcov == 'nonparametric') {
    se_beta <- apply(nonpar_results$t[,(p1 + 1):(p1+p2), drop = FALSE], 2, sd, na.rm = TRUE)
    lower_beta <- apply(nonpar_results$t[,(p1 + 1):(p1+p2), drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE)
    upper_beta <- apply(nonpar_results$t[,(p1 + 1):(p1+p2), drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)
  } else if (vcov == 'wild') {
    se_beta <- apply(beta_boot, 2, sd, na.rm = TRUE)
    lower_beta <- apply(beta_boot, 2, quantile, probs = 0.025, na.rm = TRUE)
    upper_beta <- apply(beta_boot, 2, quantile, probs = 0.975, na.rm = TRUE)
  } else if (vcov == 'fwb') {
    se_beta <- apply(fwb_results$t[,(p1 + 1):(p1+p2), drop = FALSE], 2, sd, na.rm = TRUE)
    lower_beta <- apply(fwb_results$t[,(p1 + 1):(p1+p2), drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE)
    upper_beta <- apply(fwb_results$t[,(p1 + 1):(p1+p2), drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)
  }

  ci_beta <- data.frame(name = names(beta_est),
                        Lower = lower_beta,
                        Upper = upper_beta)

  # for phi
  if(family == 'nb'){
    se_phi <- sqrt(cov_matrix[p1+p2+1, p1+p2+1])
    z <- qnorm(0.975)
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

  # bias corrected estimator
  if (!is.null(bias_corr)){

    if (vcov %in% c('hessian', 'robust') & bias_corr == 'with_alpha_bias'){
      df_boot <- cbind(data.frame(m = m, n = n, N = N), as.data.frame(X), as.data.frame(Z))
      fwb_results <- fwb(data = df_boot, statistic = mle_fit_fwb, R = 1000, verbose = FALSE)
      bias_alpha <- colMeans(fwb_results$t[,1:p1, drop = FALSE]) - alpha_est
    } else if (vcov %in% c('nonparametric', 'fwb', 'wild') & bias_corr == 'with_alpha_bias'){
      bias_alpha <- switch(vcov,
         nonparametric = colMeans(nonpar_results$t[,1:p1, drop = FALSE]) - alpha_est,
         fwb = colMeans(fwb_results$t[,1:p1, drop = FALSE]) - alpha_est,
         wild = colMeans(alpha_boot) - alpha_est
      )
    }

    cov_matrix_alpha <- cov_matrix[1:p1,1:p1, drop = FALSE]
    var_alpha_est <- diag(X %*% cov_matrix_alpha %*% t(X))

    bias_approx <- switch(bias_corr,
      with_alpha_bias = N^as.vector(X %*% alpha_est) * log(N) * (X %*% bias_alpha) + (N^as.vector(X %*% alpha_est) * log(N)^2 * (var_alpha_est + (X %*% bias_alpha)^2))/2,
      no_alpha_bias = (N^as.vector(X %*% alpha_est) * log(N)^2 * var_alpha_est)/2)  # when bias_alpha_est = 0

    # cat('Summary as.vector(X %*% alpha_est):\n')
    # print(summary(as.vector(X %*% alpha_est)))
    # cat('Summary bias_approx (as vector):\n')
    # print(summary(bias_approx))
    # cat('Sum of bias_approx:\n')

    xi_est_bc <- if (!is.null(bias_corr)) xi_est - sum(bias_approx) else NULL

    # # dla poprawionego:
    # # if vcov = hessian / robust:
    # #     nie ma znaczenia dla plug-in estymatora
    # if (vcov == 'nonparametric'){
    #   nonpar_results$t[,ncol(nonpar_results$t)] <- nonpar_results$t[,ncol(nonpar_results$t)] - bias_approx
    # } else if (vcov == 'fwb'){
    #   fwb_results$t[,ncol(fwb_results$t)] <- fwb_results$t[,ncol(fwb_results$t)] - bias_approx
    # } else if (vcov == 'wild'){
    #   xi_boot <- xi_boot - bias_approx
    # }
    # # próbka poprzednia bootstrapowa <-  próbka poprzednia bootstrapowa - bias_approx
  }

  # standard error and confidence intervals for xi
  if (vcov %in% c('hessian', 'robust')){
    N_Xalpha <- as.numeric(N)^as.vector(X %*% alpha_est)
    grad_g_alpha <- t(X) %*% (log(N) * N_Xalpha)
    se_xi <- se_alpha %*% abs(grad_g_alpha)

    conf_int_xi <- data.frame(Lower = sum(as.numeric(N)^as.vector(X %*% lower_alpha)),
                              Upper = sum(as.numeric(N)^as.vector(X %*% upper_alpha)))
  } else if (vcov == 'nonparametric') {
    se_xi <- sd(nonpar_results$t[, ncol(nonpar_results$t)])
    conf_int_xi<- data.frame(Lower = quantile(nonpar_results$t[, ncol(nonpar_results$t)], probs = 0.027),
                             Upper = quantile(nonpar_results$t[, ncol(nonpar_results$t)], probs = 0.975))
  } else if (vcov == 'wild') {
    se_xi <- sd(xi_boot)
    conf_int_xi <- data.frame(Lower = quantile(xi_boot, probs = 0.027),
                              Upper = quantile(xi_boot, probs = 0.975))
  } else if (vcov == 'fwb') {
    se_xi <- sd(fwb_results$t[, ncol(fwb_results$t)])
    conf_int_xi<- data.frame(Lower = quantile(fwb_results$t[, ncol(fwb_results$t)], probs = 0.027),
                             Upper = quantile(fwb_results$t[, ncol(fwb_results$t)], probs = 0.975))
  }

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
                  xi_est_bc = xi_est_bc,
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
                  X = X,
                  Z = Z,
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


