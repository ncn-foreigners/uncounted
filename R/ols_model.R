ols_model <- function(m, n, N){

  df <- data.frame(
    y = m,
    x1 = log(N),
    x2 = log(n/N)
  )

  df_ols <- subset(df, y>0)       # exclude zeros if they occur
  ols_fit <- lm(log(y) ~ x1 + x2 - 1, data = df_ols)
  alpha_est <- unname(coef(ols_fit)[1])
  beta_est <- unname(coef(ols_fit)[2])
  xi_est <- sum(N^alpha_est)

  #confidence intervals for alpha
  confint_alpha <- confint(ols_fit)[1,]
  #confidence intervals for xi - M estimate
  confint_xi <- c(sum(N^confint_alpha[1]), sum(N^confint_alpha[2]))

  estimates <- c(xi = xi_est, alpha = alpha_est, beta = beta_est)

  return(
    list(estimates = estimates,
         confint_alpha = setNames(confint_alpha, c("lower", "upper")),
         confint_xi = setNames(confint_xi, c("lower", "upper")),
         ols_fit = ols_fit
         )
    )
}

