# nls model with estimates form ols as starting points

nls_model <- function(m, n, N){
  
  df <- data.frame(
    y = m,
    N = N,
    n = n
  )
  
  df_ols <- subset(df, y>0)
  df_ols$x1 <- log(df_ols$N)
  df_ols$x2 <- log(df_ols$n/df_ols$N)
  ols_fit <- lm(log(y) ~ x1 + x2 - 1, data = df_ols)
  estim_ols <- list(alpha = coef(ols_fit)[1], beta = coef(ols_fit)[2])
  
  df <- subset(df, y>0)
  nls_fit <- nls(y ~ N^alpha * (n/N)^beta, data = df, start = estim_ols)
  estim_nls <- coef(nls_fit)
  alpha_est <- estim_nls['alpha']
  beta_est <- estim_nls['beta']
  xi_est <- sum(N^alpha_est)      # target parameter - M estimator
  
  return(
    list(
      coefficients = c(alpha_est, beta_est, xi_est),
      nls_fit = nls_fit,
      ols_fit = ols_fit)
  )

}
