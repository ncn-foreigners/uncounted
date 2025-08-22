
#'@export
print.hidden <- function(results){

  cat('Number of observations:', length(results$m), '\n')
  #if (!is.null(results$uniq_county)){
  #  cat('Number of unique countries:', results$n_countries, '\n\n')   ### need to change the main function
  #}

  cat('Total observed m:', sum(results$m), '\n')
  cat('Total observed n:', sum(results$n), '\n')
  cat('Total observed N:', sum(results$N), '\n\n')

  cat('Target parameter estimate:', results$xi_est, '\n')
  cat('Target parameter confidence interval: \n')
  print(results$conf_int_xi)
  cat('\n')

  cat('Coefficients:\n')
  print(results$coefficients)

}

#'@export
summary.hidden <- function(results) {

  summary_list <- list(
    method = results$method,
    coefficients = results$coefficients,
    xi_est = results$xi_est,
    conf_int_xi = results$conf_int_xi,
    conf_int_alpha = results$conf_int_alpha,
    conf_int_beta = results$conf_int_beta,
    vcov_method = results$vcov_method,
    vcov = results$vcov,
    se_coef = results$se_coef,
    summary_stat = results$summary_stat,
    aic = results$aic,
    bic = results$bic,
    residuals = results$residuals,
    fitted = results$fitted,
    m = results$m,
    n = results$n,
    N = results$N
  )

  class(summary_list) <- 'summary.hidden'
  return(summary_list)
}


#'@export
print.summary.hidden <- function(results){

  cat('Estimation method:', results$method, '\n\n')

  cat('Number of observations:', length(results$m), '\n')
#  cat('Number of unique countries:', results$n_countries, '\n\n')

  cat('Total observed m:', sum(results$m), '\n')
  cat('Total observed n:', sum(results$n), '\n')
  cat('Total observed N:', sum(results$N), '\n\n')

  cat('Target parameter estimate:', results$xi_est, '\n')
  cat('Target parameter confidence interval:', '\n')
  print(results$conf_int_xi, row.names = FALSE)
  cat('\n')

  if (results$method == 'mle'){
    coef <- c(results$coefficients$alpha, results$coefficients$beta)
  } else {
    coef <- results$coefficients
  }

  coef_table <- data.frame(
    name = names(coef),
    Estimate = as.numeric(coef)
  )
  conf_all <- rbind(results$conf_int_alpha, results$conf_int_beta)
  coef_table <- merge(coef_table, conf_all, by = 'name')
  coef_table <- merge(coef_table, results$se_coef, by = 'name')

  coef_table$t.val <- coef_table$Estimate / coef_table$Std.error
  coef_table$p.val <- 2 * (1 - pnorm(abs(coef_table$t.val)))

  rownames(coef_table) <- coef_table$name
  coef_table$name <- NULL
  colnames(coef_table) <- c('Estimate', 'Lower', 'Upper', 'Std. Error', 't value', 'Pr(>|t|)')

  cat('Coefficients:\n')
  print(coef_table)
  cat('\n')

  # cat('Covariance matrix estimation method:', results$vcov_method,'\n')
  # cat('Covariance matrix:\n')
  # print(results$vcov)
  # cat('\n')

#  cat('Standard errors:\n')
#  print(results$se)
#  cat('\n')


  if(results$method == 'ols'){

    cat('Residual standard error:', round(results$summary_stat$resid_se, 3), 'on', results$summary_stat$df_resid, 'degrees od freedom \n')
    cat('Multiple R-squared: ', round(results$summary_stat$r_squared, 4),
        ',\t Adjusted R-squared:', round(results$summary_stat$adj_r_squared, 4), '\n')
    f_stat <- results$summary_stat$f_stat
    cat('F-statistic: ', round(f_stat[1], 2), 'on', f_stat[2], 'and', f_stat[3], 'DF, ',
        'p-value: ', pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE), '\n')

  } else if(results$method == 'nls'){

    cat('Residual standard error:', round(results$summary_stat$resid_se, 3), 'on', results$summary_stat$df, 'degrees of freedom\n')
    cat('Number of iterations to convergence:', results$summary_stat$iter, '\n')
    cat('Achieved convergence:', results$summary_stat$convergence, '\n')

  }  else if(results$method == 'glm - Poisson'){

    # null deviance
    cat('Null deviance:', round(results$summary_stat$null_deviance,4), 'on', results$summary_stat$df_null, 'degrees of freedom\n')
    # residual deviance
    cat('Residual deviance:', round(results$summary_stat$resid_deviance, 4), 'on', results$summary_stat$df_resid, 'degrees of freedom \n')
    # number of Fisher scoring iterations
    cat('Number of Fisher scoring iterations:', results$summary_stat$iter, '\n')

  } else if(results$method == 'mle'){

    if (results$summary_stat$convergence == 0) {
      cat('Convergence achieved.\n')
    } else {
      cat('Warning: convergence not achieved.\n')
    }
  }

  cat('AIC:', results$aic, '\t BIC:', results$bic, '\n\n')

}


