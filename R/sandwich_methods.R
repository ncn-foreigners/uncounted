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
#'   qqnorm quantile sd
#' @importFrom grDevices adjustcolor devAskNewPage
#' @importFrom graphics abline arrows barplot legend lines mtext par
#'   plot.new rect text
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
