
#' @export
dfbeta <- function(object) {
  UseMethod('dfbeta')
}

#' @export
dfbeta.hidden <- function(object){

  C <- nrow(object$data)
  coeff_full <- object$coefficients
  p <- length(coeff_full)

  refit <- function(x, remove) {
    args <- x$call

    new_data <- x$data[-remove, , drop=FALSE]
    args$data <- new_data
    args$bias_corr <- FALSE

    eval(args)
  }

  res <- matrix(NA, nrow = C, ncol = p)
  colnames(res) <- names(coeff_full)

  for (k in 1:C) {
    suppressMessages(
      fit_k <- refit(object, remove = k)
    )
    coeff_k <- fit_k$coefficients
    res[k, ] <- coeff_full - coeff_k
  }

  res
}


#' @export
dfpopsize <- function(object, dfbeta) {
  UseMethod('dfpopsize')
}

#' @export
dfpopsize.hidden <- function(object, dfbeta = NULL){

  if(is.null(dfbeta)){
    dfbeta <- dfbeta_test(object)
  }

  varname <- function(name) {
    if (inherits(name, 'call')) {
      return(all.vars(name)[1])
    } else if (is.character(name)) {
      return(name)
    }
  }

  N <- object$data[,varname(object$call$reference_pop)]
  C <- length(N)
  if(!is.null(object$call$cov_alpha)) X <- model.matrix(as.formula(object$call$cov_alpha), object$data)

  xi_full <- object$xi_est
  coef <- object$coefficients

  # we only select param values for alpha
  alpha_id <- grepl('^alpha', names(coef))
  alpha_names <- names(coef)[alpha_id]
  alpha_full <- coef[alpha_names]
  # differences from dfbeta function
  res_alpha <- dfbeta[, alpha_names, drop=FALSE]

  res <- rep(NA, C)

  for (k in 1:C) {
    # we calculate leave-one-out alpha value using dfbeta
    loo_alpha <- alpha_full - unname(res_alpha[k, ])
    # we calculate leave-one-out value for xi:
    # we remove kth row from dataset and use loo_alpha value
    alpha_k <- if(!is.null(object$call$cov_alpha)) as.vector(X[-k, , drop=FALSE] %*% loo_alpha) else loo_alpha
    xi_k <- sum(N[-k]^alpha_k)
    res[k] <- xi_full - xi_k
  }

  res
}
