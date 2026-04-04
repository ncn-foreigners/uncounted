# ---- Prediction method ----

#' Predict from an uncounted model
#'
#' Returns predicted values from a fitted \code{uncounted} model, optionally
#' for new data. This method enables compatibility with the
#' \code{marginaleffects} package.
#'
#' @param object An \code{uncounted} object.
#' @param newdata Optional data frame for predictions. If \code{NULL},
#'   returns fitted values from the original data.
#' @param type Character: \code{"response"} (default) returns predicted counts
#'   \eqn{\hat{\mu}_i = \exp(\mathbf{z}_i'\hat{\boldsymbol{\theta}})};
#'   \code{"link"} returns the linear predictor
#'   \eqn{\mathbf{z}_i'\hat{\boldsymbol{\theta}}}.
#' @param ... Additional arguments (ignored).
#' @return A numeric vector of predicted values.
#'
#' @examples
#' data(irregular_migration)
#' d <- irregular_migration[irregular_migration$year == "2019", ]
#' fit <- estimate_hidden_pop(d, ~ m, ~ n, ~ N, method = "poisson",
#'                            gamma = 0.005)
#' # Fitted values
#' head(predict(fit))
#'
#' # Predictions for new data
#' head(predict(fit, newdata = d[1:5, ]))
#'
#' # Linear predictor (log scale)
#' head(predict(fit, type = "link"))
#'
#' @export
predict.uncounted <- function(object, newdata = NULL,
                               type = c("response", "link"), ...) {
  type <- match.arg(type)

  if (is.null(newdata)) {
    if (type == "response") return(object$fitted.values)
    return(object$log_mu)
  }

  # Reconstruct design matrix for new data
  call_args <- object$call_args
  observed_var <- all.vars(call_args$observed)
  auxiliary_var <- all.vars(call_args$auxiliary)
  refpop_var <- all.vars(call_args$reference_pop)

  N_new <- newdata[[refpop_var]]
  n_new <- newdata[[auxiliary_var]]
  ratio_new <- n_new / N_new

  gamma_val <- object$gamma
  if (!is.null(gamma_val)) {
    log_rate_new <- log(gamma_val + ratio_new)
  } else {
    log_rate_new <- log(ratio_new)
  }
  log_N_new <- log(N_new)

  # Build covariate design matrices
  cov_alpha_formula <- object$call$cov_alpha
  cov_beta_formula <- object$call$cov_beta

  if (!is.null(cov_alpha_formula)) {
    cov_alpha_formula <- eval(cov_alpha_formula)
    X_alpha_new <- .build_model_matrix(cov_alpha_formula, newdata)
  } else {
    X_alpha_new <- matrix(1, nrow = nrow(newdata), ncol = 1)
    colnames(X_alpha_new) <- "alpha"
  }

  if (!is.null(cov_beta_formula)) {
    cov_beta_formula <- eval(cov_beta_formula)
    X_beta_new <- .build_model_matrix(cov_beta_formula, newdata)
  } else {
    X_beta_new <- matrix(1, nrow = nrow(newdata), ncol = 1)
    colnames(X_beta_new) <- "beta"
  }

  # Build full design matrix
  Z_new <- cbind(X_alpha_new * log_N_new, X_beta_new * log_rate_new)

  # Handle constrained models
  alpha_coefs <- object$alpha_coefs
  beta_coefs <- object$beta_coefs
  is_constr <- isTRUE(object$constrained)

  eta_alpha <- as.numeric(X_alpha_new %*% alpha_coefs)
  eta_beta <- as.numeric(X_beta_new %*% beta_coefs)

  if (is_constr) {
    alpha_vals <- .inv_logit(eta_alpha)
    beta_vals <- exp(eta_beta)
  } else {
    alpha_vals <- eta_alpha
    beta_vals <- eta_beta
  }

  log_mu <- alpha_vals * log_N_new + beta_vals * log_rate_new
  if (type == "link") return(as.numeric(log_mu))
  as.numeric(exp(log_mu))
}
