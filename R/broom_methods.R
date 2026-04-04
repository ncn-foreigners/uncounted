# ---- broom/modelsummary integration ----
#
# Provides tidy() and glance() methods for uncounted objects,
# enabling compatibility with modelsummary, parameters, and performance.

#' @importFrom generics tidy glance
NULL

#' Tidy an uncounted model
#'
#' Extracts a tidy coefficient table from a fitted \code{uncounted} object.
#' Compatible with \code{modelsummary::modelsummary()}.
#'
#' @param x An \code{uncounted} object.
#' @param conf.int Logical; include confidence intervals? Default FALSE.
#' @param conf.level Confidence level for intervals (default 0.95).
#' @param ... Additional arguments (ignored).
#' @return A data frame with columns: \code{term}, \code{estimate},
#'   \code{std.error}, \code{statistic}, \code{p.value}, and optionally
#'   \code{conf.low}, \code{conf.high}.
#'
#' @examples
#' data(irregular_migration)
#' d <- irregular_migration[irregular_migration$year == "2019", ]
#' fit <- estimate_hidden_pop(d, ~ m, ~ n, ~ N, method = "poisson",
#'                            gamma = 0.005)
#' tidy(fit)
#' tidy(fit, conf.int = TRUE)
#'
#' @export
tidy.uncounted <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  coefs <- x$coefficients
  se <- sqrt(diag(x$vcov))
  stat <- coefs / se
  p <- if (x$method %in% c("ols", "nls")) {
    2 * pt(abs(stat), df = x$df.residual, lower.tail = FALSE)
  } else {
    2 * pnorm(abs(stat), lower.tail = FALSE)
  }
  result <- data.frame(
    term = names(coefs),
    estimate = as.numeric(coefs),
    std.error = as.numeric(se),
    statistic = as.numeric(stat),
    p.value = as.numeric(p),
    stringsAsFactors = FALSE
  )
  if (conf.int) {
    crit <- if (x$method %in% c("ols", "nls")) {
      qt((1 + conf.level) / 2, df = x$df.residual)
    } else {
      qnorm((1 + conf.level) / 2)
    }
    result$conf.low <- as.numeric(coefs - crit * se)
    result$conf.high <- as.numeric(coefs + crit * se)
  }
  result
}


#' Glance at an uncounted model
#'
#' Returns a one-row data frame summarizing model fit, compatible with
#' \code{modelsummary::modelsummary()}.
#'
#' @param x An \code{uncounted} object.
#' @param ... Additional arguments (ignored).
#' @return A one-row data frame with columns: \code{method}, \code{nobs},
#'   \code{logLik}, \code{AIC}, \code{BIC}, \code{deviance},
#'   \code{df.residual}, \code{gamma}, \code{theta}.
#'
#' @examples
#' data(irregular_migration)
#' d <- irregular_migration[irregular_migration$year == "2019", ]
#' fit <- estimate_hidden_pop(d, ~ m, ~ n, ~ N, method = "poisson",
#'                            gamma = 0.005)
#' glance(fit)
#'
#' @export
glance.uncounted <- function(x, ...) {
  ll <- logLik(x)
  data.frame(
    method = toupper(x$method),
    nobs = x$n_obs,
    logLik = as.numeric(ll),
    AIC = AIC(x),
    BIC = BIC(x),
    deviance = deviance(x),
    df.residual = x$df.residual,
    gamma = if (!is.null(x$gamma)) x$gamma else NA_real_,
    theta = if (!is.null(x$theta)) x$theta else NA_real_,
    stringsAsFactors = FALSE
  )
}
