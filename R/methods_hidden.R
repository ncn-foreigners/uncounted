#' @title Print method for \code{hidden} object
#'
#' @description Print method for the object of class \code{hidden}.
#'
#' @param x an object of class \code{hidden} returned by the main estimation function \code{estimate_hidden_pop()}
#'
#' @examples
#' \dontrun{
#' # Load Polish irregular migration data
#' data(foreigners_pl)
#'
#' model_data <- subset(foreigners_pl, year == 2018 & half == 1)
#'
#' # Basic Zhang model estimation using GLM (recommended)
#' result_zhang <- estimate_hidden_pop(
#'   data = model_data,
#'   observed = ~ border,
#'   auxiliary = ~ police,
#'   reference_pop = ~ pesel,
#'   method = "mle",
#'   family = "poisson"
#' )
#'
#' print(result_zhang)
#'
#' }
#'
#' @export
print.hidden <- function(x){

  varname <- function(name) {
    if (inherits(name, 'call')) {
      return(all.vars(name)[1])
    } else if (is.character(name)) {
      return(name)
    }
  }
  m <- x$data[,varname(x$call$observed)]
  n <- x$data[,varname(x$call$auxiliary)]
  N <- x$data[,varname(x$call$reference_pop)]
  countries <- if(!is.null(x$call$countries)) x$data[,varname(x$call$countries)] else NULL

  cat('Number of observations:', length(m), '\n')
  if (!is.null(countries)){
   cat('Number of unique countries:', length(unique(countries)), '\n\n')
  }

  cat('Total observed m:', sum(m), '\n')
  cat('Total observed n:', sum(n), '\n')
  cat('Total observed N:', sum(N), '\n\n')

  cat('Target parameter estimate:', x$xi_est, '\n')
  if(!is.null(x$xi_bias_corr)){
    cat('Bias-corrected target parameter estimate: \n')
    print(x$xi_bias_corr)
    cat('\n')
  }
  cat('Target parameter standard error:', x$se_xi, '\n')
  cat('Target parameter confidence interval: \n')
  print(x$conf_int_xi, row.names = FALSE)
  cat('\n')

  cat('Coefficients:\n')
  print(x$coefficients)

  if (!is.null(x$by_covariates)){
    cat('\n')
    cat('Results by covariates: \n')
    print(x$by_covariates)
  }

}



#' @title Summary method for class 'hidden'
#'
#' @description The function constructs a summary list for objects of class \code{hidden}.
#'
#' @param object an object of class \code{hidden} returned by the main estimation function \code{estimate_hidden_pop()}
#'
#' @return An object of class \code{summary.hidden}, which is a list containing:
#' \itemize{
#' \item \code{method} - estimation method used
#' \item \code{coefficients} - estimation coefficients
#' \item \code{xi_est} - estimate of the target parameter
#' \item \code{conf_int_xi}, \code{conf_int_alpha}, \code{conf_int_beta} - confidence intervals
#' \item \code{vcov}, \code{vcov_method} - estimated variance-covariance matrix and method used
#' \item \code{se_coef} - standard errors for coefficients
#' \item \code{se_xi} - standard error for target parameter
#' \item \code{summary stat} - list of model-specific statistics
#' \item \code{aic}, \code{bic} - model selection criteria values
#' \item \code{residuals}, \code{fitted}, \code{m}, \code{n}, \code{N} - fitted values, residuals and observed data
#' }
#'
#' @details This function is typically used to prepare a structured summary that can then be printed using \code{print.summary.hidden()}.
#'Depending on the chosen estimation method in \code{estimate_hidden_pop()} function the following model statistics are calculated:
#'\itemize{
#' \item \strong{'ols'}: residual standard error, degrees of freedom, R-squared, adjusted R-squared, F-statistic and its p-value;
#' \item \strong{'nls'}: residual standard error, degrees of freedom, number of iterations, convergence status;
#' \item \strong{'glm'} - Poisson: null deviance, residual deviance, degrees of freedom, number of Fisher scoring iterations;
#' \item \strong{'mle'} - Zhang model: convergence status.
#' }
#' Other elements such as estimated coefficients, confidence intervals, standard errors, variance-covariance matrix, AIC, BIC, fitted values and residuals
#' are also included, regardless of the estimation method.
#'
#' @examples
#' \dontrun{
#' # Load Polish irregular migration data
#' data(foreigners_pl)
#'
#' model_data <- subset(foreigners_pl, year == 2018 & half == 1)
#'
#' # Basic Zhang model estimation using GLM (recommended)
#' result_zhang <- estimate_hidden_pop(
#'   data = model_data,
#'   observed = ~ border,
#'   auxiliary = ~ police,
#'   reference_pop = ~ pesel,
#'   method = "mle",
#'   family = "poisson"
#' )
#'
#' summary(result_zhang)
#'
#' }
#'
#' @export
summary.hidden <- function(object) {

  coef_table <- data.frame(name = names(object$coefficients),
                           Estimate = as.numeric(object$coefficients))

  ci_coef <- object$conf_int_coef
  ci_coef$name <- rownames(ci_coef)
  se_coef <- object$se_coef
  se_coef$name <- rownames(se_coef)

  coef_table <- merge(coef_table, ci_coef, by = 'name')
  coef_table <- merge(coef_table, se_coef, by = 'name')

  coef_table$Estimate <- round(coef_table$Estimate, 5)
  coef_table$Lower <- round(coef_table$Lower, 5)
  coef_table$Upper <- round(coef_table$Upper, 5)
  coef_table$Std.error <- round(coef_table$Std.error, 5)
  coef_table$t.val <- round(coef_table$Estimate / coef_table$Std.error,3)
  coef_table$p.val <- format.pval(2 * (1 - pnorm(abs(coef_table$t.val))),  digits = 4, eps = 2e-16)

  rownames(coef_table) <- coef_table$name
  coef_table$name <- NULL
  colnames(coef_table) <- c('Estimate', 'Lower', 'Upper', 'Std. Error', 't value', 'Pr(>|t|)')
  coefficients <- coef_table

  coef_alpha <- coefficients[startsWith(rownames(coefficients),'alpha'),]

  coef_beta <- coefficients[startsWith(rownames(coefficients),'beta'),]

  coef_phi <- if ('phi'== tail(rownames(coefficients), 1)) coefficients[rownames(coefficients)=='phi',] else NULL

  summary_list <- list(
    coef_alpha = coef_alpha,
    coef_beta = coef_beta,
    coef_phi = coef_phi,
    xi_est = object$xi_est,
    xi_est_bc = object$xi_est_bc,
    conf_int_xi = object$conf_int_xi,
    vcov = object$vcov,
    se_xi = object$se_xi,
    summary_stat = object$summary_stat,
    residuals = object$residuals,
    fitted = object$fitted,
    by_nationality = object$by_nationality,
    by_covariates = object$by_covariates,
    call = object$call,
    data = object$data
  )

  class(summary_list) <- 'summary.hidden'
  return(summary_list)
}




#' @title Print method for summary.hidden objects
#'
#' @description Print method for the \code{summary.hidden} object.
#'
#' @param x An object of class \code{summary.hidden} returned by \code{summary.hidden()}.
#'
#' @details
#' The function displays a formatted summary of the estimation results.
#' The output content depends on the estimation method used
#' in the \code{estimate_hidden_pop()} function and subsequently on the structure created by \code{summary.hidden}.
#'
#' @examples
#' \dontrun{
#' # Load Polish irregular migration data
#' data(foreigners_pl)
#'
#' model_data <- subset(foreigners_pl, year == 2018 & half == 1)
#'
#' # Basic Zhang model estimation using GLM (recommended)
#' result_zhang <- estimate_hidden_pop(
#'   data = model_data,
#'   observed = ~ border,
#'   auxiliary = ~ police,
#'   reference_pop = ~ pesel,
#'   method = "mle",
#'   family = "poisson"
#' )
#'
#' s <- summary(result_zhang)
#' print(s)
#'
#' }
#'
#' @export
print.summary.hidden <- function(x){

  varname <- function(name) {
    if (inherits(name, 'call')) {
      return(all.vars(name)[1])
    } else if (is.character(name)) {
      return(name)
    }
  }
  m <- x$data[,varname(x$call$observed)]
  n <- x$data[,varname(x$call$auxiliary)]
  N <- x$data[,varname(x$call$reference_pop)]
  countries <- if(!is.null(x$call$countries)) x$data[,varname(x$call$countries)] else NULL

  method <- if(is.null(x$call$method)) 'ols' else x$call$method
  family <- if(is.null(x$call$family)) 'poisson' else x$call$family
  est_method <- if(method == 'mle') paste0(method, ' - ', family) else method

  cat('Estimation method:', est_method, '\n\n')

  cat('Number of observations:', length(m), '\n')
  if (!is.null(countries)){
    cat('Number of unique countries:', length(unique(countries)), '\n\n')
  }

  cat('Total observed m:', sum(m), '\n')
  cat('Total observed n:', sum(n), '\n')
  cat('Total observed N:', sum(N), '\n\n')

  cat('Target parameter estimate:', x$xi_est, '\n')
  if(!is.null(x$xi_est_bc)){
    cat('Bias-corrected target parameter estimate:', x$xi_est_bc, '\n')
  }
  cat('Target parameter standard error:', x$se_xi, '\n')
  cat('Target parameter confidence interval: \n')
  print(x$conf_int_xi, row.names = FALSE)
  cat('\n')

  if (!is.null(x$by_covariates)){
    cat('Results by covariates: \n')
    print(x$by_covariates)
    cat('\n')
  }

  cat('Coefficients (population / alpha):\n')
  print(x$coef_alpha)
  cat('\n')

  cat('Coefficients (detection / beta):\n')
  print(x$coef_beta)
  cat('\n')

  if (!is.null(x$coef_phi)){
    cat('Coefficients (heterogeneity / phi):\n')
    print(x$coef_phi)
    cat('\n')
  }

  if(method == 'ols'){

    cat('Residual standard error:', round(x$summary_stat$resid_se, 3), 'on', x$summary_stat$df_resid, 'degrees od freedom \n')
    cat('Multiple R-squared: ', round(x$summary_stat$r_squared, 4),
        ',\t Adjusted R-squared:', round(x$summary_stat$adj_r_squared, 4), '\n')
    f_stat <- x$summary_stat$f_stat
    cat('F-statistic: ', round(f_stat[1], 2), 'on', f_stat[2], 'and', f_stat[3], 'DF, ',
        'p-value: ', pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE), '\n')

  } else if(method == 'nls'){

    cat('Residual standard error:', round(x$summary_stat$resid_se, 3), 'on', x$summary_stat$df, 'degrees of freedom\n')
    cat('Number of iterations to convergence:', x$summary_stat$iter, '\n')
    cat('Achieved convergence:', x$summary_stat$convergence, '\n')

  }  else if(method == 'glm'){

    # null deviance
    cat('Null deviance:', round(x$summary_stat$null_deviance,4), 'on', x$summary_stat$df_null, 'degrees of freedom\n')
    # residual deviance
    cat('Residual deviance:', round(x$summary_stat$resid_deviance, 4), 'on', x$summary_stat$df_resid, 'degrees of freedom \n')
    # number of Fisher scoring iterations
    cat('Number of Fisher scoring iterations:', x$summary_stat$iter, '\n')

  } else if(method == 'mle'){

    if (x$summary_stat$convergence == 0) {
      cat('Convergence achieved.\n')
    } else {
      cat('Warning: convergence not achieved.\n')
    }
  }

  cat('AIC:', x$summary_stat$aic, '\t BIC:', x$summary_stat$bic, '\n\n')

}



#' @export
extract <- function(object, what) {
  UseMethod("extract")
}

#' @export
extract.hidden <- function(object, what = c('all', 'se', 'confint', 'vcov')) {

  what_selected <- match.arg(what)

  ext <- switch(what_selected,
                'all' = data.frame(value = c(object$xi_est, object$coefficients),
                                   Std.error = c(object$se_xi, object$se_coef$Std.error),
                                   Lower = c(object$conf_int_xi$Lower, object$conf_int_coef$Lower),
                                   Upper = c(object$conf_int_xi$Upper, object$conf_int_coef$Upper),
                                   row.names = c('xi', rownames(object$se_coef))),
                'se' = data.frame(Std.error = c(object$se_xi, object$se_coef$Std.error),
                                  row.names = c('xi', rownames(object$se_coef))),
                'confint' = data.frame(Lower = c(object$conf_int_xi$Lower, object$conf_int_coef$Lower),
                                       Upper = c(object$conf_int_xi$Upper, object$conf_int_coef$Upper),
                                       row.names = c('xi', rownames(object$se_coef))),
                'vcov' = object$vcov)

  return(ext)
}


# coef(summary)
#' @export
coef.summary.hidden <- function(x){
  method <- if(is.null(x$call$method)) 'ols' else x$call$method
  coef_table <- if (method == 'mle') rbind(x$coef_alpha, x$coef_beta, x$coef_phi) else rbind(x$coef_alpha, x$coef_beta)
  return(coef_table)
}


#' @export
vcov.hidden <- function(x){
  return(x$vcov)
}


# AIC
#' @export
AIC.hidden <- function(x){
  return(x$summary_stat$aic)
}


# BIC
#' @export
BIC.hidden <- function(x){
  return(x$summary_stat$bic)
}

