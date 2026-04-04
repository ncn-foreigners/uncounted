#' Check if the column space of A is a subspace of B
#' @noRd
.is_col_subspace <- function(A, B) {
  if (ncol(A) > ncol(B)) return(FALSE)
  # Project each column of A onto column space of B
  Q <- qr(B)
  proj <- qr.fitted(Q, A)
  max(abs(A - proj)) < 1e-8
}

#' Compare multiple uncounted models
#'
#' Produces a side-by-side comparison table for two or more fitted
#' \code{uncounted} models. The table includes the following metrics for each
#' model:
#'
#' \describe{
#'   \item{\code{Method}}{Estimation method (OLS, NLS, POISSON, NB).}
#'   \item{\code{Constrained}}{Whether the model was fitted with parameter
#'     constraints.}
#'   \item{\code{n_par}}{Number of estimated parameters.}
#'   \item{\code{logLik}}{Log-likelihood (true for count models, Gaussian
#'     pseudo-log-likelihood for OLS/NLS).}
#'   \item{\code{AIC, BIC}}{Akaike and Bayesian information criteria.}
#'   \item{\code{Deviance}}{Model deviance (see
#'     \code{\link{deviance.uncounted}}).}
#'   \item{\code{Pearson_X2}}{Pearson chi-squared statistic
#'     \eqn{\sum(m_i - \hat\mu_i)^2 / V(\hat\mu_i)}.}
#'   \item{\code{RMSE}}{Root mean squared error of response residuals.}
#'   \item{\code{R2_cor}}{Squared correlation between observed and fitted
#'     values: \eqn{\mathrm{cor}(m, \hat\mu)^2}.}
#'   \item{\code{R2_D}}{Explained deviance pseudo \eqn{R^2}:
#'     \eqn{1 - D(\text{model}) / D(\text{null})}, where the null model is the
#'     same specification without covariates (single \eqn{\alpha}, single
#'     \eqn{\beta}). Only for Poisson/NB. See Cameron and Trivedi (2013).}
#'   \item{\code{R2_CW}}{Cameron--Windmeijer (1996) pseudo \eqn{R^2}:
#'     \eqn{1 - \sum (m_i - \hat\mu_i)^2 / V(\hat\mu_i) \big/
#'     \sum (m_i - \bar m)^2 / V(\bar m)}.
#'     Uses the model-implied variance function \eqn{V(\mu)}.
#'     Only for Poisson/NB.}
#' }
#'
#' A warning is issued when comparing OLS/NLS pseudo-log-likelihoods with true
#' count-model likelihoods, since AIC/BIC values are not comparable across
#' these families.
#'
#' @param ... Two or more objects of class \code{"uncounted"}. Named arguments
#'   are used as model labels.
#' @param sort_by Column to sort by: \code{"AIC"} (default), \code{"BIC"}, or
#'   \code{"loglik"}.
#' @return An object of class \code{"uncounted_comparison"} containing a
#'   \code{data.frame} in \code{$table} and the model objects in \code{$models}.
#'
#' @examples
#' set.seed(42)
#' df <- data.frame(
#'   N = rep(1000, 30),
#'   n = rpois(30, lambda = 50)
#' )
#' df$m <- rpois(30, lambda = df$N^0.5 * (df$n / df$N)^0.8)
#'
#' fit_po <- estimate_hidden_pop(data = df, observed = ~m, auxiliary = ~n,
#'                               reference_pop = ~N, method = "poisson")
#' fit_nb <- estimate_hidden_pop(data = df, observed = ~m, auxiliary = ~n,
#'                               reference_pop = ~N, method = "nb")
#'
#' comp <- compare_models(Poisson = fit_po, NB = fit_nb, sort_by = "AIC")
#' comp
#'
#' @export
compare_models <- function(..., sort_by = c("AIC", "BIC", "loglik")) {
  sort_by <- match.arg(sort_by)
  models <- list(...)

  if (length(models) < 2) stop("Need at least 2 models to compare")

  # Check for mixed likelihood types
  methods <- sapply(models, function(x) x$method)
  has_ols <- any(methods %in% c("ols", "nls"))
  has_ml <- any(methods %in% c("poisson", "nb"))
  if (has_ols && has_ml) {
    warning("Comparing OLS/NLS pseudo-loglik with true likelihood is not meaningful for AIC/BIC.")
  }

  # Build table
  tab <- do.call(rbind, lapply(seq_along(models), function(i) {
    x <- models[[i]]
    ll <- logLik(x)
    n_par <- attr(ll, "df")
    n <- x$n_obs
    dev <- tryCatch(deviance(x), error = function(e) NA_real_)

    # Pearson chi-squared
    v <- .variance_function(x)
    pearson <- sum((x$m - x$fitted.values)^2 / v)

    rmse <- sqrt(mean(x$residuals^2))

    model_name <- if (!is.null(names(models)[i]) && nchar(names(models)[i]) > 0) {
      names(models)[i]
    } else {
      paste0("Model ", i)
    }

    # Pseudo R^2 measures
    r2_cor <- cor(x$m, x$fitted.values)^2
    r2_dev <- .pseudo_r2_deviance(x)
    r2_cw <- .pseudo_r2_cameron_windmeijer(x)

    data.frame(
      Model = model_name,
      Method = toupper(x$method),
      Constrained = isTRUE(x$constrained),
      n_par = n_par,
      logLik = as.numeric(ll),
      AIC = -2 * as.numeric(ll) + 2 * n_par,
      BIC = -2 * as.numeric(ll) + log(n) * n_par,
      Deviance = dev,
      Pearson_X2 = pearson,
      RMSE = rmse,
      R2_cor = r2_cor,
      R2_D = r2_dev,
      R2_CW = r2_cw,
      stringsAsFactors = FALSE
    )
  }))

  # Sort
  ord <- switch(sort_by,
    "AIC" = order(tab$AIC),
    "BIC" = order(tab$BIC),
    "loglik" = order(tab$logLik, decreasing = TRUE)
  )
  tab <- tab[ord, , drop = FALSE]
  rownames(tab) <- NULL

  structure(list(table = tab, models = models), class = "uncounted_comparison")
}

#' @export
print.uncounted_comparison <- function(x, ...) {
  cat("Model comparison\n")
  cat(strrep("-", 60), "\n")
  tab <- x$table
  tab$logLik <- round(tab$logLik, 2)
  tab$AIC <- round(tab$AIC, 2)
  tab$BIC <- round(tab$BIC, 2)
  tab$Deviance <- round(tab$Deviance, 2)
  tab$Pearson_X2 <- round(tab$Pearson_X2, 2)
  tab$RMSE <- round(tab$RMSE, 2)
  tab$R2_cor <- round(tab$R2_cor, 4)
  tab$R2_D <- round(tab$R2_D, 4)
  tab$R2_CW <- round(tab$R2_CW, 4)
  print(tab, row.names = FALSE)
  invisible(x)
}

#' Likelihood ratio test for nested uncounted models
#'
#' Performs a likelihood ratio test between two nested \code{uncounted} models
#' fitted by maximum likelihood (Poisson or Negative Binomial). The test
#' statistic is
#'
#' \deqn{LR = 2\bigl(\ell_2 - \ell_1\bigr)}
#'
#' where \eqn{\ell_1} and \eqn{\ell_2} are the maximised log-likelihoods of
#' the simpler and more complex models, respectively. Under \eqn{H_0} the
#' statistic follows a \eqn{\chi^2} distribution with degrees of freedom equal
#' to the difference in the number of estimated parameters.
#'
#' \strong{Boundary correction for Poisson vs NB.}
#' When comparing a Poisson model (H0) against a Negative Binomial model (H1),
#' the dispersion parameter \eqn{\theta} lies on the boundary of its parameter
#' space under \eqn{H_0} (\eqn{\theta \to \infty}). The standard
#' \eqn{\chi^2} approximation is invalid in this setting. Following Self &
#' Liang (1987), the p-value is corrected as
#'
#' \deqn{p = 0.5\,\Pr(\chi^2_1 > LR)}
#'
#' which corresponds to a 50:50 mixture of a point mass at zero and a
#' \eqn{\chi^2_1}.
#'
#' This function is intended for comparing count models only (Poisson and NB).
#' OLS/NLS models are not supported because they use a pseudo-log-likelihood.
#' The models must be fitted to the same data (identical number of observations).
#'
#' @param object1 An \code{"uncounted"} model (Poisson or NB).
#' @param object2 An \code{"uncounted"} model (Poisson or NB). The function
#'   automatically identifies the simpler model regardless of argument order.
#' @return An object of class \code{"uncounted_lrtest"} with components
#'   \code{statistic}, \code{df}, \code{p_value}, \code{boundary}, model
#'   descriptions, and log-likelihoods.
#'
#' @references
#' Self, S. G. and Liang, K.-Y. (1987). Asymptotic properties of maximum
#' likelihood estimators and likelihood ratio tests under nonstandard
#' conditions. \emph{Journal of the American Statistical Association},
#' 82(398), 605--610.
#'
#' @examples
#' set.seed(42)
#' df <- data.frame(
#'   N = rep(1000, 30),
#'   n = rpois(30, lambda = 50)
#' )
#' df$m <- rpois(30, lambda = df$N^0.5 * (df$n / df$N)^0.8)
#'
#' fit_po <- estimate_hidden_pop(data = df, observed = ~m, auxiliary = ~n,
#'                               reference_pop = ~N, method = "poisson")
#' fit_nb <- estimate_hidden_pop(data = df, observed = ~m, auxiliary = ~n,
#'                               reference_pop = ~N, method = "nb")
#'
#' # Test for overdispersion (Poisson vs NB)
#' lrtest(fit_po, fit_nb)
#'
#' @export
lrtest <- function(object1, object2) {
  if (!inherits(object1, "uncounted") || !inherits(object2, "uncounted"))
    stop("Both objects must be of class 'uncounted'")

  if (is.null(object1$loglik) || is.null(object2$loglik))
    stop("Both models must have a log-likelihood (not OLS/NLS)")

  if (object1$n_obs != object2$n_obs)
    stop("Models must be fitted to the same data (different n_obs)")

  # Ensure object2 is the more complex model (more parameters)
  p1 <- .count_params(object1)
  p2 <- .count_params(object2)
  if (p1 > p2) {
    tmp <- object1; object1 <- object2; object2 <- tmp
  }

  ll1 <- object1$loglik
  ll2 <- object2$loglik
  lr_stat <- max(0, 2 * (ll2 - ll1))  # clamp to 0 if simpler model fits better
  df <- .count_params(object2) - .count_params(object1)

  if (df == 0) {
    stop("Models have the same number of parameters -- not nested.")
  }

  # Boundary correction for Poisson vs NB (Self & Liang, 1987)
  is_boundary <- (object1$method == "poisson" && object2$method == "nb") ||
                 (object1$method == "nb" && object2$method == "poisson")

  # Warn for different methods (only Poisson vs NB is a valid cross-method LR)
  if (object1$method != object2$method && !is_boundary) {
    warning("Models use different methods (", object1$method, " vs ",
            object2$method, "). LR test may not be valid.", call. = FALSE)
  }

  # Warn for different constrained settings
  if (!identical(isTRUE(object1$constrained), isTRUE(object2$constrained))) {
    warning("Models have different constraint settings. ",
            "LR test may not be valid.", call. = FALSE)
  }

  # Check nesting via column-space inclusion
  nested_alpha <- .is_col_subspace(object1$X_alpha, object2$X_alpha)
  nested_beta  <- .is_col_subspace(object1$X_beta, object2$X_beta)
  gamma_ok <- identical(object1$gamma, object2$gamma) ||
              (is.numeric(object1$gamma) && isTRUE(object2$gamma_estimated)) ||
              (isTRUE(object1$gamma_estimated) && isTRUE(object2$gamma_estimated))
  if (!(nested_alpha && nested_beta && gamma_ok)) {
    warning("Models may not be nested: covariate or gamma specifications differ. ",
            "LR test is only valid for nested models.", call. = FALSE)
  }

  if (is_boundary && df == 1) {
    p_value <- 0.5 * pchisq(lr_stat, df = 1, lower.tail = FALSE)
  } else {
    p_value <- pchisq(lr_stat, df = df, lower.tail = FALSE)
  }

  result <- list(
    statistic = lr_stat,
    df = df,
    p_value = p_value,
    boundary = is_boundary,
    model1 = paste(toupper(object1$method),
                   if(object1$constrained) "(constrained)" else ""),
    model2 = paste(toupper(object2$method),
                   if(object2$constrained) "(constrained)" else ""),
    ll1 = ll1, ll2 = ll2
  )
  class(result) <- "uncounted_lrtest"
  result
}

#' @export
print.uncounted_lrtest <- function(x, ...) {
  cat("Likelihood ratio test\n")
  cat(strrep("-", 40), "\n")
  cat("Model 1:", x$model1, " (logLik =", round(x$ll1, 2), ")\n")
  cat("Model 2:", x$model2, " (logLik =", round(x$ll2, 2), ")\n")
  cat("LR statistic:", round(x$statistic, 4), "on", x$df, "df\n")
  if (x$boundary) {
    cat("(Boundary-corrected: 0.5 * P(chi2 > LR), Self & Liang 1987)\n")
  }
  cat("p-value:", format.pval(x$p_value, digits = 4), "\n")
  invisible(x)
}


# ---- Pseudo R^2 helpers (internal) ----

#' Explained deviance R^2
#'
#' R^2_D = 1 - D(model) / D(null), where null = same model without covariates
#' (single alpha, single beta). Only for Poisson/NB models.
#' @noRd
.pseudo_r2_deviance <- function(object) {
  if (!object$method %in% c("poisson", "nb")) return(NA_real_)

  # Fit null model (no covariates) on the same data
  null_fit <- tryCatch({
    estimate_hidden_pop(
      data = object$data,
      observed = object$call_args$observed,
      auxiliary = object$call_args$auxiliary,
      reference_pop = object$call_args$reference_pop,
      method = object$method,
      gamma = if (object$gamma_estimated) "estimate" else object$gamma,
      constrained = isTRUE(object$constrained)
    )
  }, error = function(e) NULL)

  if (is.null(null_fit)) return(NA_real_)

  dev_model <- tryCatch(deviance(object), error = function(e) NA_real_)
  dev_null <- tryCatch(deviance(null_fit), error = function(e) NA_real_)

  if (is.na(dev_model) || is.na(dev_null) || dev_null == 0) return(NA_real_)
  1 - dev_model / dev_null
}

#' Cameron-Windmeijer (1996) pseudo R^2 for count data
#'
#' R^2_CW = 1 - sum((m - mu)^2 / V(mu)) / sum((m - mbar)^2 / V(mbar))
#' Uses the model-implied variance function V(mu).
#' @references Cameron, A.C. and Windmeijer, F.A.G. (1996). R-squared measures
#'   for count data regression models with applications to health-care
#'   utilization. \emph{Journal of Business & Economic Statistics}, 14(2),
#'   209--220.
#' @noRd
.pseudo_r2_cameron_windmeijer <- function(object) {
  if (!object$method %in% c("poisson", "nb")) return(NA_real_)

  m <- object$m
  mu <- object$fitted.values
  mbar <- mean(m)

  v_model <- .variance_function(object)
  # Variance at mean: use same V() but evaluated at mbar
  if (object$method == "poisson") {
    v_null <- rep(mbar, length(m))
  } else {
    theta <- object$theta
    v_null <- rep(mbar + mbar^2 / theta, length(m))
  }

  num <- sum((m - mu)^2 / v_model)
  den <- sum((m - mbar)^2 / v_null)

  if (den == 0) return(NA_real_)
  1 - num / den
}
