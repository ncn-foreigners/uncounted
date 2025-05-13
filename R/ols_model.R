ols_model <- function(m, n, N){
  
  df <- data.frame( 
    y = m,
    x1 = log(N),
    x2 = log(n/N)
  )
  
  df_ols <- subset(df, y>0)                          # exclude zeros if they occur
  ols_fit <- lm(log(y) ~ x1 + x2 - 1, data = df_ols)  
  estim <- coef(ols_fit)
  alpha_est <- estim[1]
  beta_est <- estim[2]
  xi_est <- sum(N^alpha_est)      # target parameter - M estimator
  
  return(
         list(coefficients = c(alpha_est, beta_est,xi_est), 
              conf_int_xi = c(sum(N^confint(ols_fit)[1,1]), sum(N^confint(ols_fit)[1,2])))
          # confidence intervals like Zhang 2008
         )
  
}
