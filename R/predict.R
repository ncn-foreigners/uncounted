# ---- Prediction method ----

# Build design matrix for prediction, aligned with training columns.
# Uses the training data's model.frame to extract terms/xlev for
# consistent factor encoding (same reference level, same dummies).
#' @noRd
.predict_design <- function(cov_formula, newdata, X_train, prefix,
                            train_data) {
  train_cols <- colnames(X_train)
  n <- nrow(newdata)

  if (is.null(cov_formula)) {
    X <- matrix(1, nrow = n, ncol = 1)
    colnames(X) <- train_cols
    return(X)
  }

  # Strip prefix from training column names to get raw model.matrix names
  raw_train_cols <- sub(paste0("^", prefix, ":"), "", train_cols)

  # Build terms + xlev from training data for consistent encoding
  mf_train <- model.frame(cov_formula, data = train_data)
  tt <- terms(mf_train)
  xlev <- .getXlevels(tt, mf_train)

  X <- tryCatch(
    model.matrix(tt, data = newdata, xlev = xlev),
    error = function(e) NULL
  )
  if (is.null(X)) {
    stop("Cannot build design matrix from newdata. Ensure all covariate ",
         "levels present in newdata match the training data.", call. = FALSE)
  }

  # Align: match raw model.matrix columns to raw training columns
  raw_new_cols <- colnames(X)
  result <- matrix(0, nrow = n, ncol = length(train_cols))
  colnames(result) <- train_cols
  for (j in seq_along(train_cols)) {
    idx <- match(raw_train_cols[j], raw_new_cols)
    if (!is.na(idx)) {
      result[, j] <- X[, idx]
    }
  }
  result
}


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
#' head(predict(fit))
#' head(predict(fit, newdata = d[1:5, ]))
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

  # Extract variable names from the original call
  call_args <- object$call_args
  refpop_var <- all.vars(call_args$reference_pop)
  auxiliary_var <- all.vars(call_args$auxiliary)

  N_new <- newdata[[refpop_var]]
  n_new <- newdata[[auxiliary_var]]
  ratio_new <- n_new / N_new

  # Compute gamma for new data
  cov_gamma_f <- if (!is.null(object$call$cov_gamma)) eval(object$call$cov_gamma) else NULL
  if (!is.null(cov_gamma_f) && !is.null(object$gamma_coefs) &&
      isTRUE(object$has_cov_gamma)) {
    X_gamma_new <- .predict_design(cov_gamma_f, newdata, object$X_gamma, "gamma", object$data)
    gamma_vals <- exp(as.numeric(X_gamma_new %*% object$gamma_coefs))
  } else {
    gamma_val <- object$gamma
    gamma_vals <- if (!is.null(gamma_val)) rep(gamma_val, length(ratio_new)) else NULL
  }
  rate_new <- .rate_from_gamma(ratio_new, gamma_vals)
  log_N_new <- log(N_new)

  # Build covariate design matrices aligned with training encoding
  cov_alpha_f <- if (!is.null(object$call$cov_alpha)) eval(object$call$cov_alpha) else NULL
  cov_beta_f <- if (!is.null(object$call$cov_beta)) eval(object$call$cov_beta) else NULL

  X_alpha_new <- .predict_design(cov_alpha_f, newdata, object$X_alpha, "alpha", object$data)
  X_beta_new <- .predict_design(cov_beta_f, newdata, object$X_beta, "beta", object$data)

  # Compute linear predictors
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

  log_mu <- .compute_log_mu(alpha_vals, log_N_new, beta_vals, rate_new,
                            link_rho = object$link_rho)
  if (type == "link") return(as.numeric(log_mu))
  as.numeric(exp(log_mu))
}
