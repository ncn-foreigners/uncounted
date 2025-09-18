#' @importFrom stats model.matrix
#' @importFrom stats lm
#' @importFrom stats nls
#' @importFrom stats glm
#' @importFrom stats poisson
#' @importFrom stats vcov
#' @importFrom stats coef
#' @importFrom stats optim
#'
#' @title Estimate hidden population size using non-linear count regression model
#'
#' @description This function estimates hidden population sizes using the Zhang (2008) methodology
#' for irregular migration populations. It implements a non-linear count regression
#' model specifically designed for estimating irregular residents based on administrative
#' data from border control, police records, and population registers
#'
#' @details The function implements the methodology developed by Zhang (2008) and adapted by
#' Beręsewicz & Pawlukiewicz (2020) for estimating irregular foreigner populations.
#' The core model assumes that observed apprehensions follow a Poisson distribution:
#'
#' \deqn{m_{it} \sim \text{Poisson}(\lambda_{it})}
#'
#' where \eqn{\lambda_{it} = \mu_{it} u_{it}} and:
#'
#' \deqn{\mu_{it} = N_{it}^{\alpha} \left(\frac{n_{it}}{N_{it}}\right)^{\beta}}
#'
#' The target parameter (theoretical size of irregular residents) is:
#' \deqn{\xi_t = \sum_{i=1}^{C} E(M_{it}|N_{it}) = \sum_{i=1}^{C} N_{it}^{\alpha}}
#'
#' where \eqn{\alpha} and \eqn{\beta} are parameters that should theoretically lie in \eqn{(0,1)}.
#'
#' @param data A `data.frame` containing the variables specified in the other parameters
#' @param observed A formula (e.g., `~ m_var`) or character string specifying the
#'   variable name for observed irregular residents apprehended by border guards.
#'   This represents \eqn{m_{it}} in the Zhang model -- the number of irregular foreigners
#'   detected through border control operations.
#' @param auxiliary A formula (e.g., `~ n_var`) or character string specifying the
#'   variable name for auxiliary information from police records. This represents
#'   \eqn{n_{it}} in the Zhang model -- the number of foreigners appearing in police
#'   administrative records (legally staying).
#' @param reference_pop A formula (e.g., `~ N_var`) or character string specifying
#'   the variable name for reference population size from population register.
#'   This represents \eqn{N_{it}} in the Zhang model - the number of foreigners registered
#'   in the central population register (PESEL in Poland).
#' @param cov_alpha A formula (e.g., `~ x1 + x2`) specifying covariates for the
#'   \eqn{\alpha} parameter. Allows \eqn{\alpha} to vary by groups (e.g., country, sex, time).
#'   Only used with `method = "mle"`. Default is `NULL`.
#' @param cov_beta A formula (e.g., `~ x1 + x2`) specifying covariates for the
#'   \eqn{\beta} parameter. Allows \eqn{\beta} to vary by groups (e.g., country, sex, time).
#'   Only used with `method = "mle"`. Default is `NULL`.
#' @param method Character string specifying the estimation method. Options are:
#'   \itemize{
#'     \item `"ols"` -- Ordinary least squares on linearized model
#'     \item `"nls"` -- Non-linear least squares estimation
#'     \item `"mle"` -- Generalized linear model with Poisson family (Zhang method)
#'   }
#'   Default is `"ols"`. The `"mle"` methods are recommended for
#'   actual estimation as they properly account for the Poisson distribution.
#' @param vcov Character string specifying the variance-covariance estimation method.
#'   Options are `"hessian"` (default) or `"robust"` for robust standard errors.
#' @param family Character string specifying the error distribution family for GLM
#'   methods. For the Zhang model, should be `"poisson"` (default: `"gaussian"`).
#'
#' @return A list containing estimation results with the following structure:
#'   \itemize{
#'     \item For `method = "ols"` or `method = "nls"`:
#'       \itemize{
#'         \item `estimates` --
#'         \item `vcov` - Variance-covariance matrix
#'         \item `se` - Standard errors
#'       }
#'     \item For `method = "glm"`:
#'       \itemize{
#'         \item `estimates` - List with `alpha` and `beta` parameter vectors
#'         \item `target_estimate` - Estimated total hidden population \eqn{\hat{\xi}}
#'         \item `vcov` - Variance-covariance matrix
#'         \item `se` - Standard errors
#'       }
#'   }
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

#' # Linearized model for assumption verification
#' result_linear <- estimate_hidden_pop(
#'   data = model_data,
#'   observed = ~ border,
#'   auxiliary = ~ police,
#'   reference_pop = ~ pesel,
#'   method = "ols"
#' )
#'
#' # Access results
#' print(result_zhang$target_estimate)  # Total estimated irregular population
#' print(result_zhang$estimates)        # Parameter estimates
#' }
#'
#' @references
#' Beręsewicz, M., & Pawlukiewicz, K. (2020). Estimation of the number of irregular
#' foreigners in Poland using non-linear count regression models. *arXiv preprint*
#' arXiv:2008.09407.
#'
#' Zhang, L.-C. (2008). Developing methods for determining the number of unauthorized
#' foreigners in Norway. Statistics Norway (SSB), Division for Statistical Methods
#' and Standards.
#'
#' Beręsewicz, M., Gudaszewski, G., and Szymkowiak, M. (2019). Estymacja liczby
#' cudzoziemców w Polsce z wykorzystaniem metody capture-recapture. *Wiadomości
#' Statystyczne. The Polish Statistician*, 64(10), 7-35.
#'
#' @seealso
#'
#' \code{\link{foreigners_pl}} for the example dataset.
#'
#' @export
estimate_hidden_pop <- function(
    data,
    observed,           # ~ m  formula for m (allow for one)
    auxiliary,          # ~ n  allow for one
    reference_pop,      # ~ N  reference population size
    cov_alpha = NULL,   # ~ x1 + x2  allow for covariates
    cov_beta  = NULL,   # ~ x1 + x2  allow for covariates
    method = 'ols',     # various method for fit / ols (linear regression), nls, glm (poisson regression), mle
    vcov = 'hessian',     # or robust - (method to estimate standard errors)
    family = 'gaussian',     # poisson, nb etc
    countries = NULL   # string column name with country of origin names, optional, but useful in plots
){

  if (!method %in% c("ols", "nls", 'glm', "mle")) {
    stop("Invalid estimation method. Choose 'ols', 'nls', 'glm', 'mle'.")
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
  countries <- if(!is.null(countries)) data[[countries]] else NULL

  # covariates
  X <- if (is.null(cov_alpha)==FALSE){model.matrix(cov_alpha, data)} else NULL
  Z <- if (is.null(cov_beta)==FALSE){model.matrix(cov_beta, data)} else NULL
  # columns to use in estimates by covariates
  cov_vars_alpha <-  if(!is.null(cov_alpha)) all.vars(cov_alpha) else character(0)
  cov_vars_beta <-  if(!is.null(cov_beta)) all.vars(cov_beta) else character(0)
  covariate_vars <- unique(c(cov_vars_alpha , cov_vars_beta))
  df_cov <- if(length(covariate_vars) > 0) df_cov <- data[, covariate_vars, drop = FALSE] else NULL



  if (any(is.na(m)) || any(is.na(n)) || any(is.na(N)) || any(is.na(X)) || any(is.na(Z))) {
    stop('Input data contains missing values (NA). Please handle missing data before estimation.')
  }

  if (any(N <= 0)) {
    stop('Variable ', N_var, ' contains zero or negative values, which are not allowed.')
  }

  if (any(n <= 0)) {
    stop('Variable ', n_var, ' contains zero or negative values, which are not allowed.')
  }

  if (any(N < n)) {
    stop('There are observations where ', N_var,  ' < ', n_var, ' , which is not allowed.')
  }


  if (is.list(m) || is.list(n) || is.list(N)) {
    stop(sprintf("Columns %s, %s, and %s cannot be lists.", m_var, n_var, N_var))
  }
  if (!is.numeric(m) || !is.numeric(n) || !is.numeric(N)) {
    stop(sprintf("Columns %s, %s, and %s must be numeric. Check input data types.", m_var, n_var, N_var))
  }


  results <- switch(method,
                    'ols' = ols_model(m, n, N, vcov = vcov, countries),
                    'nls' = nls_model(m, n, N, vcov = vcov, countries),
                    'glm' = glm_model(m, n, N, vcov = vcov, countries),
                    'mle' = zhang_model_cov(m, n, N, X, Z, vcov = vcov, countries, df_cov))


  if (method %in% c('ols', 'nls', 'glm')){
    if (is.null(X)==FALSE){
      message("Covariates won't be included in the ", method," estimation method.")
    }
    if(results$coefficients[1] < 0 || results$coefficients[1] > 1 ){
      message('Estimated alpha parameter lies out of the expected interval (0,1).')
    }
    if(results$coefficients[2] < 0 || results$coefficients[2] > 1 ){
      message('Estimated beta parameter lies out of the expected interval (0,1).')
    }
  }

  if (method == 'mle') {

    alpha_vals <- results$coefficients[startsWith(names(results$coefficients), 'alpha')]
    beta_vals <- results$coefficients[startsWith(names(results$coefficients), 'beta')]

    if (any(alpha_vals < 0 | alpha_vals > 1)) {
      message('The estimated parameter(s) alpha are outside the expected interval (0,1).')
      if (length(alpha_vals) > 1) {
        if (sum(alpha_vals) > 0 && sum(alpha_vals) < 1) {
          message('The sum of the estimated parameters alpha lies within the expected interval (0,1).')
        } else {
          message('The sum of the estimated parameters alpha lies outside the expected interval (0,1).')
        }
      }
    }
    if (any(beta_vals < 0 | beta_vals > 1)) {
      message('The estimated parameter(s) beta are outside the expected interval (0,1).')
      if (length(beta_vals) > 1) {
        if (sum(beta_vals) > 0 && sum(beta_vals) < 1) {
          message('The sum of the estimated parameters beta lies within the expected interval (0,1).')
        } else {
          message('The sum of the estimated parameters beta lies outside the expected interval (0,1).')
        }
      }
    }
  }


  class(results) <- 'hidden'

  return(results)

}
