
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
    se = results$se,
    iter = results$iter,
    convergence = results$convergence,
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
    estimate = as.numeric(coef)
  )
  conf_all <- rbind(results$conf_int_alpha, results$conf_int_beta)
  coef_table <- merge(coef_table, conf_all, by = 'name')
  rownames(coef_table) <- coef_table$name
  coef_table$name <- NULL

  cat('Coefficients:\n')
  print(coef_table)
  cat('\n')

  cat('Covariance matrix estimation method:', results$vcov_method,'\n')
  cat('Covariance matrix:\n')
  print(results$vcov)
  cat('\n')

#  cat('Standard errors:\n')
#  print(results$se)
#  cat('\n')

  if (!is.null(results$iter)){
    cat('Number of iterations to convergence:', results$iter, '\n')
    cat('Achieved convergence:', results$convergence, '\n\n')
  }

  cat('AIC:', results$aic, '\t BIC:', results$bic, '\n\n')

}


