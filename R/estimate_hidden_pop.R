#' Estimate hidden population size using non-linear count regression model
#'
#' This function estimates hidden population sizes using the Zhang (2008) methodology
#' for irregular migration populations. It implements a non-linear count regression
#' model specifically designed for estimating irregular residents based on administrative
#' data from border control, police records, and population registers.
#'
#' @description
#' The function implements the methodology developed by Zhang (2008) and adapted by
#' Beręsewicz & Pawlukiewicz (2020) for estimating irregular foreigner populations.
#' The core model assumes that observed apprehensions follow a Poisson distribution:
#'
#' $$m_{it} \sim \text{Poisson}(\lambda_{it})$$
#'
#' where $\lambda_{it} = \mu_{it} u_{it}$ and:
#'
#' $$\mu_{it} = N_{it}^{\alpha} \left(\frac{n_{it}}{N_{it}}\right)^{\beta}$$
#'
#' The target parameter (theoretical size of irregular residents) is:
#' $$\xi_t = \sum_{i=1}^{C} E(M_{it}|N_{it}) = \sum_{i=1}^{C} N_{it}^{\alpha}$$
#'
#' where $\alpha$ and $\beta$ are parameters that should theoretically lie in $(0,1)$.
#'
#' @param data A `data.frame` containing the variables specified in the other parameters
#' @param observed A formula (e.g., `~ m_var`) or character string specifying the
#'   variable name for observed irregular residents apprehended by border guards.
#'   This represents $m_{it}$ in the Zhang model -- the number of irregular foreigners
#'   detected through border control operations.
#' @param auxiliary A formula (e.g., `~ n_var`) or character string specifying the
#'   variable name for auxiliary information from police records. This represents
#'   $n_{it}$ in the Zhang model -- the number of foreigners appearing in police
#'   administrative records (legally staying).
#' @param reference_pop A formula (e.g., `~ N_var`) or character string specifying
#'   the variable name for reference population size from population register.
#'   This represents $N_{it}$ in the Zhang model - the number of foreigners registered
#'   in the central population register (PESEL in Poland).
#' @param cov_alpha A formula (e.g., `~ x1 + x2`) specifying covariates for the
#'   $\alpha$ parameter. Allows $\alpha$ to vary by groups (e.g., country, sex, time).
#'   Only used with `method = "mle"`. Default is `NULL`.
#' @param cov_beta A formula (e.g., `~ x1 + x2`) specifying covariates for the
#'   $\beta$ parameter. Allows $\beta$ to vary by groups (e.g., country, sex, time).
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
#'         \item `estimates` - Vector: $[\hat{\xi}, \hat{\alpha}, \hat{\beta}]$ where $\hat{\xi}$ is estimated hidden population
#'         \item `vcov` - Variance-covariance matrix
#'         \item `se` - Standard errors
#'       }
#'     \item For `method = "glm"`:
#'       \itemize{
#'         \item `estimates` - List with `alpha` and `beta` parameter vectors
#'         \item `target_estimate` - Estimated total hidden population $\hat{\xi}$
#'         \item `vcov` - Variance-covariance matrix
#'         \item `se` - Standard errors
#'       }
#'   }
#'
#' @details
#' **The Zhang Methodology**: This approach is specifically designed for irregular
#' migration estimation using administrative data. Unlike traditional capture-recapture
#' methods, it uses a functional relationship between observed apprehensions,
#' auxiliary police data, and reference population.
#'
#' **Model Assumptions**:
#' 1. $M_t$ = size of unauthorized resident population (random variable)
#' 2. $N_t$ = size of known reference population (fixed, known)
#' 3. Target parameter: $\xi_t = E(M_t|N_t)$ (theoretical size)
#' 4. Detection probability structure: $\omega_i = E(p_{it}|n_{it}, N_{it}) = (n_{it}/N_{it})^{\beta}$
#'
#' **Data Requirements**: The model requires that $m_{tij} > 0$, $n_{tij} > 0$ and
#' $n_{tij}/N_{tij} < 1$ for all observations included in estimation.
#'
#' **Model Verification**: The linearized model for assumption checking:
#' $$\log\left(\frac{m_i}{N_i}\right) = (\alpha - 1)\log N_i + \beta \log\left(\frac{n_i}{N_i}\right) + \epsilon_i$$
#'
#' Expected relationships: negative coefficient for $\log N_i$, positive for $\log(n_i/N_i)$.
#'
#' **Covariates Extension**: When covariates are included:
#' $$\alpha_i = \exp(X_i^T \gamma_\alpha), \quad \beta_i = \exp(Z_i^T \gamma_\beta)$$
#'
#' **Parameter Interpretation**:
#' - $\alpha$: Population scaling parameter (should be in $(0,1)$)
#' - $\beta$: Detection probability parameter (should be in $(0,1)$)
#' - Values outside $(0,1)$ may indicate model misspecification
#'
#' @section Data Sources (Polish Context):
#' The method was developed using three administrative data sources:
#' \itemize{
#'   \item **Border Guard data**: Apprehensions of irregular foreigners ($m$)
#'   \item **Police data**: Foreign nationals in police records ($n$)
#'   \item **PESEL register**: Central population register ($N$)
#' }
#'
#' @section Warnings:
#' The function provides warnings when:
#' \itemize{
#'   \item Covariates are specified for OLS or NLS methods (covariates are ignored)
#'   \item Estimated $\alpha$ or $\beta$ parameters lie outside $(0,1)$
#'   \item For multiple parameters, when individual or summed parameters are outside $(0,1)$
#' }
#'
#' @section Mathematical Foundation:
#' The methodology is based on the relationship:
#' $$E(M_{it}|N_{it}) = N_{it}^{\alpha}$$
#'
#' and detection probability:
#' $$E(p_{it}|n_{it}, N_{it}) = \left(\frac{n_{it}}{N_{it}}\right)^{\beta}$$
#'
#' The total irregular population estimate is:
#' $$\hat{\xi} = \sum_{i=1}^{C} N_i^{\hat{\alpha}}$$
#'
#' @importFrom stats model.matrix
#' @importFrom stats lm nls glm poisson
#' @importFrom stats vcov coef optim
#'
#' @examples
#' \dontrun{
#' # Load Polish irregular migration data
#' data(foreigners_pl)
#'
#' model_data <- subset(foreigners_pl, year == 2018 & half == 1)
#' model_data <- aggregate(cbind(border, police, pesel) ~ country, model_data, sum)
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
#' \code{\link{ols_model}}, \code{\link{nls_model}}, \code{\link{zhang_model_cov}}
#' for the underlying estimation functions.
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
    method = 'ols',     # various method for fit / ols (linear regression), nls, mle
    vcov = 'hessian',     # or robust - (method to estimate standard errors)
    family = 'gaussian'     # poisson, nb etc
){

  if (!method %in% c("ols", "nls", "mle")) {
    stop("Invalid estimation method. Choose 'ols', 'nls', 'mle'.")
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


  results <- switch(method,
                    'ols' = ols_model(m, n, N, vcov = vcov),
                    'nls' = nls_model(m, n, N, vcov = vcov),
                    'mle' = zhang_model_cov(m, n, N, X, Z, vcov = vcov))


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

  if (method == 'mle') {
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
