# ---- broom/modelsummary integration ----
#
# Provides tidy() and glance() methods for uncounted objects,
# enabling compatibility with modelsummary, parameters, and performance.

#' @importFrom generics tidy glance
#' @export
generics::tidy
#' @export
generics::glance

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

  # Add gamma row(s) if estimated
  # When cov_gamma: gamma coefs are already in x$coefficients
  # and appear in the main table. Only add separate gamma row for scalar case.
  if (isTRUE(x$gamma_estimated) && !is.null(x$gamma) &&
      !isTRUE(x$has_cov_gamma)) {
    # vcov_full["gamma","gamma"] is the variance of log(gamma) (optimizer scale).
    # Delta method: SE(gamma) = gamma * SE(log_gamma)
    gamma_se_log <- NA_real_
    if (!is.null(x$vcov_full)) {
      gamma_idx <- which(colnames(x$vcov_full) == "gamma")
      if (length(gamma_idx) == 1) {
        gamma_se_log <- sqrt(max(0, x$vcov_full[gamma_idx, gamma_idx]))
      }
    }
    gamma_se <- if (is.finite(gamma_se_log)) x$gamma * gamma_se_log else NA_real_
    # No natural null hypothesis for gamma (positive nuisance parameter),
    # so statistic and p-value are not reported.
    gamma_stat <- NA_real_
    gamma_p <- NA_real_
    gamma_row <- data.frame(
      term = "gamma",
      estimate = x$gamma,
      std.error = gamma_se,
      statistic = gamma_stat,
      p.value = gamma_p,
      stringsAsFactors = FALSE
    )
    if (conf.int) {
      # CI via exp of log-scale CI: gamma * exp(+/- crit * SE_log)
      gamma_row$conf.low <- if (is.finite(gamma_se_log)) {
        x$gamma * exp(-crit * gamma_se_log)
      } else NA_real_
      gamma_row$conf.high <- if (is.finite(gamma_se_log)) {
        x$gamma * exp(crit * gamma_se_log)
      } else NA_real_
    }
    result <- rbind(result, gamma_row)
  }

  # Add theta row for NB
  if (!is.null(x$theta)) {
    theta_row <- data.frame(
      term = "theta",
      estimate = x$theta,
      std.error = if (!is.null(x$theta_se)) x$theta_se else NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      stringsAsFactors = FALSE
    )
    if (conf.int) {
      theta_row$conf.low <- NA_real_
      theta_row$conf.high <- NA_real_
    }
    result <- rbind(result, theta_row)
  }

  rownames(result) <- NULL
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
    stringsAsFactors = FALSE
  )
}
