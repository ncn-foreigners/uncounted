
estimate_hidden_pop <- function(
  data,
  observed,           # ~ m  formula for m (allow for one)
  auxiliary,          # ~ n  allow for one
  reference_pop,      # ~ N  reference population size
  cov_alpha = NULL,  # ~ x1 + x2  allow for covariates
  cov_beta  = NULL,  # ~ x1 + x2  allow for covariates
  method = 'ols',         # various method for fit / ols (linear regression), nls, glm
  vcov = 'hessian',     # or robust - (method to estimate standard errors)
  family = 'gaussian'     # poisson, nb etc
){

  if (!method %in% c("ols", "nls", "glm", "zhang")) {
    stop("Invalid estimation method. Choose 'ols', 'nls', 'glm'.")
  }

  if (!vcov %in% c("hessian", "robust")) {
    stop("Invalid vcov method. Choose 'hessian' or 'robust'.")
  }

  variable <- function(name) {

    if (inherits(name, 'formula')) {
      all.vars(name)[1]
    } else if (is.character(name)) {
      name
    } else {
      stop('Incorrect argument was given: observed, auxiliary, reference_pop must be a formula or a character string (variable name).')
    }
  }

  m_var <- variable(observed)
  n_var <- variable(auxiliary)
  N_var <- variable(reference_pop)

  vars_needed <- c(m_var, n_var, N_var)
  vars_missing <- vars_needed[!vars_needed %in% colnames(data)]

  if (length(vars_missing) > 0) {
    stop("The following variable(s) are missing in 'data': ", paste(vars_missing, collapse = ", "))
  }

  m <- data[[m_var]]
  n <- data[[n_var]]
  N <- data[[N_var]]


  if (is.list(m) || is.list(n) || is.list(N)) {
    stop(sprintf("Columns %s, %s, and %s cannot be lists.", m_var, n_var, N_var))
  }
  if (!is.numeric(m) || !is.numeric(n) || !is.numeric(N)) {
    stop(sprintf("Columns %s, %s, and %s must be numeric. Check input data types.", m_var, n_var, N_var))
  }

  # covariates
  X <- if (is.null(cov_alpha)==FALSE){model.matrix(cov_alpha, data)} else NULL
  Z <- if (is.null(cov_beta)==FALSE){model.matrix(cov_beta, data)} else NULL


  results <- switch(method, 'ols' = ols_model(m, n, N, vcov = vcov),
                      'nls' = nls_model(m, n, N, vcov = vcov),
                      'glm' = zhang_model_cov(m, n, N, X, Z, vcov = vcov),
                      'zhang' = zhang_model_cov(m, n, N, X, Z, vcov = vcov))


  if (method == 'ols'){
    if (is.null(X)==FALSE){
      message("Covariates won't be included in the 'ols' estimation method.")
    }
    if(results$estimates[2] < 0 || results$estimates[2] > 1 ){
      message('Estimated alpha parameter lies out of the expected interval (0,1).')
    }
    if(results$estimates[3] < 0 || results$estimates[3] > 1 ){
      message('Estimated beta parameter lies out of the expected interval (0,1).')
    }
  }

  if (method == 'nls'){
    if (is.null(X)==FALSE){
      message("Covariates won't be included in the 'nls' estimation method.")
    }
    if(results$estimates[2] < 0 || results$estimates[2] > 1 ){
      message('Estimated alpha parameter lies out of the expected interval (0,1).')
    }
    if(results$estimates[3] < 0 || results$estimates[3] > 1 ){
      message('Estimated beta parameter lies out of the expected interval (0,1).')
    }
  }

  if (method == 'glm' || method == 'zhang') {
    if (any(results$estimates$alpha < 0 | results$estimates$alpha > 1)) {
      message('The estimated parameter(s) alpha are outside the expected interval (0,1).')
      if (length(results$estimates$alpha) > 1) {
        if (sum(results$estimates$alpha) > 0 && sum(results$estimates$alpha) < 1) {
          message('The sum of the estimated parameters alpha lies within the expected interval (0,1).')
        } else {
          message('The sum of the estimated parameters alpha lies outside the expected interval (0,1).')
        }
      }
    }
    if (any(results$estimates$beta < 0 | results$estimates$beta > 1)) {
      message('The estimated parameter(s) beta are outside the expected interval (0,1).')
      if (length(results$estimates$beta) > 1) {
        if (sum(results$estimates$beta) > 0 && sum(results$estimates$beta) < 1) {
          message('The sum of the estimated parameters beta lies within the expected interval (0,1).')
        } else {
          message('The sum of the estimated parameters beta lies outside the expected interval (0,1).')
        }
      }
    }
  }

  return(results)

}

