#' Log-Likelihood for uncounted models
#'
#' Extracts the log-likelihood from a fitted \code{uncounted} model.
#'
#' For Poisson and Negative Binomial models, the true log-likelihood is returned
#' as computed during maximum-likelihood estimation. For OLS and NLS models,
#' a Gaussian pseudo-log-likelihood is returned, computed on the log scale as
#'
#' \deqn{\ell = -\frac{n}{2}\bigl(\log(2\pi\hat\sigma^2) + 1\bigr)}
#'
#' where \eqn{\hat\sigma^2 = n^{-1}\sum(\log m_i - \log\hat\mu_i)^2}. This
#' pseudo-log-likelihood enables AIC/BIC comparison among OLS/NLS variants but
#' is not comparable to the true likelihood of count models.
#'
#' @param object An object of class \code{"uncounted"}.
#' @param ... Additional arguments (ignored).
#' @return An object of class \code{"logLik"} with attributes \code{df}
#'   (number of estimated parameters) and \code{nobs} (number of observations).
#' @references
#' McCullagh, P. and Nelder, J. A. (1989). \emph{Generalized Linear Models},
#' 2nd ed. Chapman & Hall.
#'
#' @examples
#' # Simulate synthetic data
#' set.seed(123)
#' df <- data.frame(
#'   N = rep(1000, 30),
#'   n = rpois(30, lambda = 50)
#' )
#' df$m <- rpois(30, lambda = df$N^0.5 * (df$n / df$N)^0.8)
#'
#' fit <- estimate_hidden_pop(data = df, observed = ~m, auxiliary = ~n,
#'                            reference_pop = ~N, method = "poisson")
#' logLik(fit)
#' AIC(fit)
#' BIC(fit)
#'
#' @export
logLik.uncounted <- function(object, ...) {
  n_par <- .count_params(object)

  if (object$method %in% c("poisson", "nb")) {
    ll <- object$loglik
  } else {
    # Gaussian pseudo-loglik on log scale (match transformation used in fitting)
    n <- object$n_obs
    log_m <- .log_m_transform(object$m)
    resid_log <- log_m - object$log_mu
    sigma2 <- sum(resid_log^2) / n
    ll <- -n/2 * (log(2 * pi * sigma2) + 1)
  }

  structure(ll, df = n_par, nobs = object$n_obs, class = "logLik")
}

#' Number of observations for uncounted models
#'
#' Returns the number of observations used to fit the model.
#'
#' @param object An object of class \code{"uncounted"}.
#' @param ... Additional arguments (ignored).
#' @return An integer scalar.
#' @export
nobs.uncounted <- function(object, ...) {
  object$n_obs
}

#' Deviance for uncounted models
#'
#' Computes the deviance (twice the difference between the saturated and fitted
#' log-likelihoods) for a fitted \code{uncounted} model.
#'
#' The formula depends on the estimation method:
#'
#' \strong{Poisson:}
#' \deqn{D = 2\sum_{i}\bigl[m_i\log(m_i/\hat\mu_i) - (m_i - \hat\mu_i)\bigr]}
#' with the convention \eqn{0 \log 0 = 0}.
#'
#' \strong{Negative Binomial:}
#' \deqn{D = 2\sum_{i}\bigl[m_i\log(m_i/\hat\mu_i) -
#'   (m_i+\theta)\log\bigl((m_i+\theta)/(\hat\mu_i+\theta)\bigr)\bigr]}
#'
#' \strong{OLS / NLS:}
#' \deqn{D = \sum_{i}(\log m_i - \log\hat\mu_i)^2}
#' i.e., the residual sum of squares on the log scale.
#'
#' @param object An object of class \code{"uncounted"}.
#' @param ... Additional arguments (ignored).
#' @return A numeric scalar giving the deviance.
#' @references
#' McCullagh, P. and Nelder, J. A. (1989). \emph{Generalized Linear Models},
#' 2nd ed. Chapman & Hall.
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(
#'   N = rep(1000, 30),
#'   n = rpois(30, lambda = 50)
#' )
#' df$m <- rpois(30, lambda = df$N^0.5 * (df$n / df$N)^0.8)
#'
#' fit_po <- estimate_hidden_pop(data = df, observed = ~m, auxiliary = ~n,
#'                               reference_pop = ~N, method = "poisson")
#' deviance(fit_po)
#'
#' @export
deviance.uncounted <- function(object, ...) {
  m <- object$m
  mu <- object$fitted.values

  switch(object$method,
    "poisson" = {
      d <- ifelse(m == 0, 2 * mu, 2 * (m * log(m / mu) - (m - mu)))
      sum(d)
    },
    "nb" = {
      theta <- object$theta
      d <- ifelse(m == 0,
                  2 * theta * log((theta + mu) / theta),
                  2 * (m * log(m / mu) - (m + theta) * log((m + theta) / (mu + theta))))
      sum(d)
    },
    {
      # OLS/NLS: RSS on log scale (match transformation used in fitting)
      log_m <- .log_m_transform(m)
      resid_log <- log_m - object$log_mu
      sum(resid_log^2)
    }
  )
}

#' Residuals for uncounted models
#'
#' Extracts residuals of the specified type from a fitted \code{uncounted} model.
#'
#' @param object An object of class \code{"uncounted"}.
#' @param type Type of residual: \code{"response"} (default), \code{"pearson"},
#'   \code{"deviance"}, or \code{"anscombe"}.
#' @param ... Additional arguments (ignored).
#' @return A numeric vector of residuals with the same length as the number of
#'   observations.
#'
#' @details
#' Four residual types are supported:
#'
#' \describe{
#'   \item{\strong{Response residuals}}{
#'     \deqn{r_i = m_i - \hat\mu_i}
#'   }
#'   \item{\strong{Pearson residuals}}{
#'     \deqn{r_i = \frac{m_i - \hat\mu_i}{\sqrt{V(\hat\mu_i)}}}
#'     where the variance function \eqn{V(\mu)} depends on the model:
#'     \itemize{
#'       \item Poisson: \eqn{V(\mu) = \mu}
#'       \item Negative Binomial: \eqn{V(\mu) = \mu + \mu^2/\theta}
#'       \item OLS/NLS: \eqn{V(\mu) = \mu^2 \hat\sigma^2} (delta-method
#'         approximation on the log scale)
#'     }
#'   }
#'   \item{\strong{Deviance residuals}}{
#'     Signed square root of the individual deviance contributions:
#'     \deqn{r_i = \mathrm{sign}(m_i - \hat\mu_i)\,\sqrt{d_i}}
#'     where \eqn{d_i} is the \eqn{i}-th term in the deviance sum (see
#'     \code{\link{deviance.uncounted}}). For OLS/NLS models, standardised
#'     log-scale residuals are returned instead.
#'   }
#'   \item{\strong{Anscombe residuals}}{
#'     Variance-stabilising residuals.
#'
#'     For Poisson (McCullagh & Nelder, 1989):
#'     \deqn{r_i = \frac{3\bigl(m_i^{2/3} - \hat\mu_i^{2/3}\bigr)}
#'       {2\,\hat\mu_i^{1/6}}}
#'
#'     For Negative Binomial (Hilbe, 2011):
#'     \deqn{r_i = \frac{3\bigl(m_i^{2/3} - \hat\mu_i^{2/3}\bigr)}
#'       {2\,\hat\mu_i^{1/6}\,(1 + \hat\mu_i/\theta)^{1/6}}}
#'
#'     For OLS/NLS: standardised log-scale residuals.
#'   }
#' }
#'
#' @references
#' McCullagh, P. and Nelder, J. A. (1989). \emph{Generalized Linear Models},
#' 2nd ed. Chapman & Hall.
#'
#' Hilbe, J. M. (2011). \emph{Negative Binomial Regression}, 2nd ed. Cambridge
#' University Press.
#'
#' Pierce, D. A. and Schafer, D. W. (1986). Residuals in generalized linear
#' models. \emph{Journal of the American Statistical Association}, 81(396),
#' 977--986.
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(
#'   N = rep(1000, 30),
#'   n = rpois(30, lambda = 50)
#' )
#' df$m <- rpois(30, lambda = df$N^0.5 * (df$n / df$N)^0.8)
#'
#' fit <- estimate_hidden_pop(data = df, observed = ~m, auxiliary = ~n,
#'                            reference_pop = ~N, method = "poisson")
#'
#' # Response residuals (default)
#' head(residuals(fit))
#'
#' # Pearson residuals
#' head(residuals(fit, type = "pearson"))
#'
#' # Deviance and Anscombe residuals
#' head(residuals(fit, type = "deviance"))
#' head(residuals(fit, type = "anscombe"))
#'
#' @export
residuals.uncounted <- function(object, type = c("response", "pearson", "deviance", "anscombe", "working"), ...) {
  type <- match.arg(type)
  m <- object$m
  mu <- object$fitted.values

  switch(type,
    "response" = object$residuals,
    "working" = object$score_residuals,
    "pearson" = {
      v <- .variance_function(object)
      (m - mu) / sqrt(v)
    },
    "deviance" = {
      switch(object$method,
        "poisson" = {
          d <- ifelse(m == 0, 2 * mu, 2 * (m * log(m / mu) - (m - mu)))
          sign(m - mu) * sqrt(pmax(d, 0))
        },
        "nb" = {
          theta <- object$theta
          d <- ifelse(m == 0,
                      2 * theta * log((theta + mu) / theta),
                      2 * (m * log(m / mu) - (m + theta) * log((m + theta) / (mu + theta))))
          sign(m - mu) * sqrt(pmax(d, 0))
        },
        {
          # OLS/NLS: standardized log residuals
          log_m <- .log_m_transform(m)
          resid_log <- log_m - object$log_mu
          resid_log / sd(resid_log)
        }
      )
    },
    "anscombe" = {
      switch(object$method,
        "poisson" = {
          # McCullagh & Nelder (1989), Zhang (2008)
          3 * (m^(2/3) - mu^(2/3)) / (2 * mu^(1/6))
        },
        "nb" = {
          # Hilbe (2011) approximation
          theta <- object$theta
          3 * (m^(2/3) - mu^(2/3)) / (2 * mu^(1/6) * (1 + mu/theta)^(1/6))
        },
        {
          # OLS/NLS: standardized log residuals
          log_m <- .log_m_transform(m)
          resid_log <- log_m - object$log_mu
          resid_log / sd(resid_log)
        }
      )
    }
  )
}

#' Variance function for GLM-type models
#' @noRd
.variance_function <- function(object) {
  mu <- object$fitted.values
  switch(object$method,
    "poisson" = mu,
    "nb" = mu + mu^2 / object$theta,
    {
      # OLS/NLS on log scale: V(m) ~ mu^2 * sigma2 via delta method
      n <- object$n_obs
      p <- length(object$coefficients) + object$gamma_estimated
      log_m <- .log_m_transform(object$m)
      resid_log <- log_m - object$log_mu
      sigma2 <- sum(resid_log^2) / (n - p)
      mu^2 * sigma2
    }
  )
}

#' Log-transform m consistently: log(m) if all positive, log(m+1) if zeros exist
#' @noRd
.log_m_transform <- function(m) {
  if (any(m == 0)) log(m + 1) else log(m)
}

#' Count total estimated parameters
#' @noRd
.count_params <- function(object) {
  p <- object$p_alpha + object$p_beta
  if (object$gamma_estimated) p <- p + 1
  if (!is.null(object$theta)) p <- p + 1  # NB theta
  if (object$method %in% c("ols", "nls")) p <- p + 1  # sigma2
  p
}

# ---- dfbeta / dfpopsize ----

#' Leave-One-Out Influence on Coefficients
#'
#' Computes the change in estimated coefficients when each observation (or
#' country) is dropped. This is a convenience wrapper around \code{\link{loo}}.
#'
#' @param model A fitted \code{uncounted} object.
#' @param by Character: \code{"obs"} (default) drops one observation at a time,
#'   \code{"country"} drops all observations for one country at a time
#'   (requires \code{countries} in the original fit).
#' @param ... Additional arguments passed to \code{\link{loo}}.
#' @return A matrix with one row per dropped unit and one column per coefficient.
#'   Each entry is the change in the coefficient estimate relative to the full
#'   model (\eqn{\hat\beta_{(-i)} - \hat\beta}).
#'
#' @seealso \code{\link{loo}}, \code{\link{dfpopsize}}
#'
#' @examples
#' data(irregular_migration)
#' d <- irregular_migration[irregular_migration$year == 2019 & irregular_migration$sex == "m", ]
#' fit <- estimate_hidden_pop(d, ~ m, ~ n, ~ N, method = "poisson",
#'                            gamma = 0.001, countries = ~ country_code)
#' db <- dfbeta(fit)
#' head(db)
#'
#' @export
dfbeta.uncounted <- function(model, by = c("obs", "country"), ...) {
  by <- match.arg(by)
  loo_res <- loo(model, by = by, ...)
  loo_res$dfbeta
}

#' Leave-One-Out Influence on Population Size
#'
#' Computes the change in the total estimated population size
#' \eqn{\hat\xi = \sum_i N_i^{\hat\alpha_i}} when each observation (or country)
#' is dropped. This is a convenience wrapper around \code{\link{loo}}.
#'
#' @param model A fitted \code{uncounted} object.
#' @param by Character: \code{"obs"} (default) or \code{"country"}.
#' @param ... Additional arguments passed to \code{\link{loo}}.
#' @return A named numeric vector. Each element is the change in \eqn{\hat\xi}
#'   when that unit is dropped (\eqn{\hat\xi_{(-i)} - \hat\xi}).
#'
#' @seealso \code{\link{loo}}, \code{\link{dfbeta.uncounted}}
#'
#' @examples
#' data(irregular_migration)
#' d <- irregular_migration[irregular_migration$year == 2019 & irregular_migration$sex == "m", ]
#' fit <- estimate_hidden_pop(d, ~ m, ~ n, ~ N, method = "poisson",
#'                            gamma = 0.001, countries = ~ country_code)
#' dp <- dfpopsize(fit)
#' head(dp)
#'
#' @export
dfpopsize <- function(model, ...) {
  UseMethod("dfpopsize")
}

#' @rdname dfpopsize
#' @export
dfpopsize.uncounted <- function(model, by = c("obs", "country"), ...) {
  by <- match.arg(by)
  loo_res <- loo(model, by = by, ...)
  dxi <- loo_res$dxi
  names(dxi) <- loo_res$dropped
  dxi
}
