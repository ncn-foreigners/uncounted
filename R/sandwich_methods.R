# ---- Sandwich package S3 methods for uncounted objects ----
#
# These methods enable compatibility with the sandwich package,
# allowing vcovHC(), vcovCL(), and sandwich() to work on uncounted objects.
#
# Methods already defined elsewhere:
#   nobs.uncounted        -> diagnostics.R
#   residuals.uncounted   -> diagnostics.R (includes type = "working")
#   vcov.uncounted        -> estimate_hidden_pop.R

#' @importFrom sandwich bread estfun
#' @importFrom stats nobs residuals hatvalues model.matrix weights dfbeta
#'   cor deviance dnbinom dpois lm lm.fit lm.wfit logLik median optim
#'   optimize pchisq pnbinom pnorm ppois printCoefmat pt qnorm qqline
#'   qqnorm quantile sd AIC BIC vcov terms model.frame .getXlevels qt
#'   weighted.mean
#' @importFrom grDevices adjustcolor devAskNewPage
#' @importFrom graphics abline arrows barplot legend lines mtext par
#'   plot.new rect text axis points segments
#' @importFrom utils head
NULL

#' Extract model matrix
#'
#' Returns the Jacobian/design matrix used in sandwich computation.
#' When gamma was estimated, the returned matrix includes the gamma column.
#'
#' @param object An uncounted object
#' @param ... Ignored
#' @export
model.matrix.uncounted <- function(object, ...) {
  object$model_matrix_full
}

#' Hat values (leverage)
#'
#' Computes diagonal of the hat matrix
#' \eqn{H = W^{1/2} Z (Z'WZ)^{-1} Z' W^{1/2}}
#' where Z is the model matrix and W is the diagonal weight matrix.
#'
#' @param model An uncounted object
#' @param ... Ignored
#' @export
hatvalues.uncounted <- function(model, ...) {
  .hat_values_wls(model$model_matrix_full, model$bread_weights)
}

#' Bread matrix for sandwich estimator
#'
#' Returns the bread component B such that the model-based variance is B/n.
#' Specifically: \eqn{B = n (Z' W Z)^{-1}} where Z is the model matrix and
#' \eqn{W = \mathrm{diag}(\mathrm{bread\_weights})}.
#'
#' @param x An uncounted object
#' @param ... Ignored
#'
#' @note For NB models, these methods operate on the mean-model parameters
#'   (alpha, beta, and optionally gamma) only, excluding theta. The stored
#'   \code{vcov()} on NB objects uses a dedicated theta-aware path instead.
#'   Calling \code{sandwich::vcovHC()} directly on an NB object gives
#'   theta-conditional standard errors, which differ from \code{vcov()}.
#'
#' @export
bread.uncounted <- function(x, ...) {
  Z <- x$model_matrix_full
  w <- x$bread_weights
  n <- x$n_obs
  wZ <- Z * sqrt(w)
  B <- n * .solve_safe(crossprod(wZ))
  rownames(B) <- colnames(B) <- colnames(Z)
  B
}

#' Empirical estimating functions (score contributions)
#'
#' Returns an n x p matrix of per-observation score contributions
#' for use with the sandwich package. Row i is the gradient of the
#' (quasi-)log-likelihood contribution of observation i with respect
#' to the parameter vector.
#'
#' @param x An uncounted object
#' @param ... Ignored
#' @export
estfun.uncounted <- function(x, ...) {
  Z <- x$model_matrix_full
  ef <- Z * x$score_residuals
  colnames(ef) <- colnames(Z)
  ef
}


#' Theta-Aware Full Variance-Covariance for NB Models
#'
#' Returns the full sandwich variance-covariance matrix for NB models,
#' including the theta (dispersion) parameter. For non-NB models, returns
#' \code{vcov(object)} unchanged.
#'
#' This is the explicit interface for the theta-aware covariance that
#' \code{vcov()} uses internally for NB fits. Use this when you need the
#' full matrix including theta, or when you want to be explicit about
#' which covariance you are requesting.
#'
#' @param object A fitted \code{uncounted} object.
#' @param vcov_type HC type for the sandwich (\code{"HC0"} or \code{"HC1"}).
#'   Default \code{"HC1"}. HC2+ are not available for the theta-aware path.
#' @param cluster Optional cluster vector for cluster-robust variance.
#' @return A square variance-covariance matrix. For NB models, dimensions
#'   are \code{p + 1} (or \code{p + 2} if gamma is estimated), where the
#'   extra row/column corresponds to \code{log(theta)}. For non-NB models,
#'   returns \code{vcov(object)}.
#'
#' @examples
#' data(irregular_migration)
#' d <- irregular_migration[irregular_migration$year == "2019", ]
#' fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "nb", gamma = 0.005)
#' vcov_nb(fit)  # includes theta row/column
#' vcov(fit)     # alpha/beta only (same as vcov_nb submatrix)
#'
#' @export
vcov_nb <- function(object, vcov_type = "HC1", cluster = NULL) {
  if (!inherits(object, "uncounted"))
    stop("'object' must be of class 'uncounted'")
  if (is.null(object$hessian_nll) || is.null(object$score_full)) {
    return(vcov(object))
  }
  .compute_nb_sandwich(
    hessian_nll = object$hessian_nll,
    score_full  = object$score_full,
    n_obs       = object$n_obs,
    vcov_type   = vcov_type,
    cluster     = cluster
  )
}
