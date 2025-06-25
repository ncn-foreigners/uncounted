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

  df_nls <- subset(df, y>0)
  nls_fit <- nls(y ~ N^alpha * (n/N)^beta, data = df_nls, start = estim_ols)
  estim_nls <- coef(nls_fit)
  alpha_est <- estim_nls['alpha']
  beta_est <- estim_nls['beta']
  xi_est <- sum(N^alpha_est)

  #confidence intervals for alpha
  confint_alpha <- confint(nls_fit)[1,]
  #confidence intervals for xi - M estimate
  confint_xi <- c(sum(N^confint_alpha[1]), sum(N^confint_alpha[2]))

  estimates <- c(xi = xi_est, alpha_est, beta_est)

  return(
    list(estimates = estimates,
         confint_alpha = setNames(confint_alpha, c("lower", "upper")),
         confint_xi = setNames(confint_xi, c("lower", "upper")),
         nls_fit = nls_fit
         )
    )

}
