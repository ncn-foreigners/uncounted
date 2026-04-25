#' Bayesian Hidden-Population Models via brms
#'
#' Fit Bayesian Poisson or negative-binomial hidden-population models using
#' \pkg{brms}. The model mirrors the count-model part of
#' \code{\link{estimate_hidden_pop}}:
#' \deqn{\log(\mu_i) = \alpha_i \log(N_i) + \log\rho_i,}
#' where \eqn{\rho_i = h(\beta_i \log(\gamma_i + n_i / N_i))}. Posterior
#' draws of \eqn{\alpha_i} are propagated to population-size draws
#' \eqn{\xi = \sum_i N_i^{\alpha_i}} by \code{\link{popsize}}.
#'
#' @param data Data frame.
#' @param observed One-sided formula naming the observed count.
#' @param auxiliary One-sided formula naming the auxiliary count.
#' @param reference_pop One-sided formula naming the reference population.
#' @param method Count likelihood. Currently \code{"poisson"} and \code{"nb"}
#'   are supported.
#' @param cov_alpha Optional one-sided formula for alpha covariates.
#' @param cov_beta Optional one-sided formula for beta covariates.
#' @param gamma Numeric fixed gamma, \code{"estimate"} for one bounded scalar
#'   gamma, or \code{NULL} to use \eqn{n_i / N_i} directly.
#' @param cov_gamma Reserved for a future release. It must be \code{NULL}.
#' @param gamma_bounds Length-two numeric bounds used when
#'   \code{gamma = "estimate"}.
#' @param link_rho Detection link: \code{"power"}, \code{"cloglog"},
#'   \code{"logit"}, or \code{"probit"}.
#' @param constrained Logical. If \code{TRUE}, alpha is mapped through an
#'   inverse-logit transform and beta through an exponential transform.
#' @param weights Optional numeric observation weights.
#' @param countries Optional one-sided formula for a country/group label stored
#'   on the fitted object.
#' @param seed,chains,iter,warmup,cores,control Passed to \code{brms::brm()}.
#' @param prior Optional \pkg{brms} prior specification. If \code{NULL}, weakly
#'   informative defaults from \code{default_priors_uncounted_bayes()} are used.
#' @param backend Stan backend passed to \code{brms::brm()}. Defaults to
#'   \code{"rstan"}; \code{"cmdstanr"} is available when installed.
#' @param ... Additional arguments passed to \code{brms::brm()}.
#'
#' @return An object of class \code{"uncounted_bayes"} containing the
#'   \code{brmsfit} and the design metadata needed by \code{\link{popsize}}.
#'
#' @examples
#' \dontrun{
#' data(irregular_migration)
#' d <- irregular_migration
#' d$ukr <- as.integer(d$country_code == "UKR")
#'
#' fit_bayes <- estimate_hidden_pop_bayes(
#'   data = d, observed = ~m, auxiliary = ~n, reference_pop = ~N,
#'   method = "poisson",
#'   cov_alpha = ~ year * ukr + sex,
#'   cov_beta = ~ year,
#'   gamma = 0.005,
#'   chains = 2, iter = 1000, seed = 1
#' )
#'
#' ps <- popsize(fit_bayes, by = ~ year + country_code, total = TRUE)
#' hypotheses_popsize(
#'   ps,
#'   "xi[year == 2024 & country_code == 'UKR'] < 15000"
#' )
#' }
#'
#' @export
estimate_hidden_pop_bayes <- function(data,
                                      observed,
                                      auxiliary,
                                      reference_pop,
                                      method = c("poisson", "nb"),
                                      cov_alpha = NULL,
                                      cov_beta = NULL,
                                      gamma = "estimate",
                                      cov_gamma = NULL,
                                      gamma_bounds = c(1e-10, 0.5),
                                      link_rho = c("power", "cloglog", "logit", "probit"),
                                      constrained = FALSE,
                                      weights = NULL,
                                      countries = NULL,
                                      seed = NULL,
                                      chains = 4,
                                      iter = 2000,
                                      warmup = floor(iter / 2),
                                      cores = getOption("mc.cores", 1L),
                                      control = NULL,
                                      prior = NULL,
                                      backend = c("rstan", "cmdstanr"),
                                      ...) {
  .require_namespace("brms")

  method <- match.arg(method)
  backend <- match.arg(backend)
  link_rho <- .normalize_link_rho(link_rho)

  if (!is.null(cov_gamma)) {
    stop("'cov_gamma' is not currently supported by estimate_hidden_pop_bayes().",
         call. = FALSE)
  }
  if (identical(backend, "cmdstanr")) {
    .require_namespace("cmdstanr")
  } else {
    .require_namespace("rstan")
  }

  call <- match.call()
  call$data <- data
  call$observed <- observed
  call$auxiliary <- auxiliary
  call$reference_pop <- reference_pop
  call$method <- method
  call$gamma <- gamma
  call$gamma_bounds <- gamma_bounds
  call$link_rho <- link_rho
  call$constrained <- constrained
  call$backend <- backend
  if (!is.null(cov_alpha)) call$cov_alpha <- cov_alpha
  if (!is.null(cov_beta)) call$cov_beta <- cov_beta
  if (!is.null(weights)) call$weights <- weights
  if (!is.null(countries)) call$countries <- countries

  spec <- .bayes_make_brms_spec(
    data = data,
    observed = observed,
    auxiliary = auxiliary,
    reference_pop = reference_pop,
    method = method,
    cov_alpha = cov_alpha,
    cov_beta = cov_beta,
    gamma = gamma,
    gamma_bounds = gamma_bounds,
    link_rho = link_rho,
    constrained = constrained,
    weights = weights,
    prior = prior
  )

  fit <- brms::brm(
    formula = spec$formula,
    data = spec$data,
    family = spec$family,
    prior = spec$prior,
    stanvars = spec$stanvars,
    backend = backend,
    seed = seed,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    control = control,
    ...
  )

  out <- list(
    fit = fit,
    engine = "brms",
    backend = backend,
    method = method,
    estimator = "bayes",
    family = spec$family,
    formula = spec$formula,
    prior = spec$prior,
    stanvars = spec$stanvars,
    gamma = spec$gamma,
    gamma_estimated = spec$gamma_estimated,
    gamma_bounds = gamma_bounds,
    constrained = constrained,
    link_rho = link_rho,
    call = call,
    data = data,
    m = spec$m,
    N = spec$N,
    n_aux = spec$n_aux,
    ratio = spec$ratio,
    X_alpha = spec$X_alpha,
    X_beta = spec$X_beta,
    p_alpha = ncol(spec$X_alpha),
    p_beta = ncol(spec$X_beta),
    cov_alpha_vars = spec$cov_alpha_vars,
    countries_var = if (!is.null(countries)) eval(countries[[2]], data) else NULL,
    obs_weights = spec$obs_weights,
    n_obs = length(spec$m),
    alpha_draws = NULL
  )
  class(out) <- "uncounted_bayes"
  out
}

#' Default Priors for Bayesian Hidden-Population Models
#'
#' @inheritParams estimate_hidden_pop_bayes
#'
#' @return A \pkg{brms} prior specification.
#'
#' @export
default_priors_uncounted_bayes <- function(method = c("poisson", "nb"),
                                           constrained = FALSE,
                                           gamma = "estimate",
                                           gamma_bounds = c(1e-10, 0.5)) {
  .require_namespace("brms")
  method <- match.arg(method)

  priors <- if (isTRUE(constrained)) {
    c(
      brms::set_prior("normal(0, 1.5)", nlpar = "alpha"),
      brms::set_prior("normal(0, 1)", nlpar = "beta")
    )
  } else {
    c(
      brms::set_prior("normal(0.7, 0.5)", nlpar = "alpha"),
      brms::set_prior("normal(1, 0.5)", nlpar = "beta")
    )
  }

  if (identical(gamma, "estimate")) {
    .validate_gamma_bounds(gamma_bounds)
    lb <- gamma_bounds[1]
    ub <- gamma_bounds[2]
    center <- stats::qlogis((mean(gamma_bounds) - lb) / (ub - lb))
    priors <- c(
      priors,
      brms::set_prior(sprintf("normal(%.8f, 2)", center), nlpar = "gamma")
    )
  }

  if (identical(method, "nb")) {
    priors <- c(priors, brms::set_prior("exponential(1)", class = "shape"))
  }

  priors
}

#' Extract the brmsfit from a Bayesian uncounted Model
#'
#' @param object An object.
#' @param ... Additional arguments, currently ignored.
#'
#' @return The underlying \code{brmsfit}.
#'
#' @export
as_brmsfit <- function(object, ...) {
  UseMethod("as_brmsfit")
}

#' @export
as_brmsfit.uncounted_bayes <- function(object, ...) {
  object$fit
}

#' Extract Posterior Population-Size Draws
#'
#' @param object A Bayesian fitted model from
#'   \code{\link{estimate_hidden_pop_bayes}} or a Bayesian population-size
#'   object from \code{\link{popsize}}.
#' @param by Optional grouping formula when \code{object} is an
#'   \code{"uncounted_bayes"} fit.
#' @param include_total Logical; if \code{TRUE}, append a \code{"Total"} draw
#'   column.
#' @param format Output format: \code{"matrix"} or \code{"data.frame"}.
#' @param ... Additional arguments passed to \code{\link{popsize}} when
#'   \code{object} is an \code{"uncounted_bayes"} fit.
#'
#' @return A draw-by-group matrix or wide data frame.
#'
#' @export
posterior_popsize_draws <- function(object, ...) {
  UseMethod("posterior_popsize_draws")
}

#' @export
posterior_popsize_draws.uncounted_bayes <- function(object, by = NULL,
                                                    include_total = FALSE,
                                                    format = c("matrix", "data.frame"),
                                                    ...) {
  ps <- popsize(object, by = by, total = include_total, ...)
  posterior_popsize_draws(ps, include_total = include_total, format = format)
}

#' @export
posterior_popsize_draws.uncounted_popsize_bayes <- function(object,
                                                            include_total = FALSE,
                                                            format = c("matrix", "data.frame"),
                                                            ...) {
  format <- match.arg(format)
  draws <- attr(object, "draws")
  if (is.null(draws)) {
    stop("This object does not contain posterior population-size draws.",
         call. = FALSE)
  }
  draws <- as.matrix(draws)

  if (include_total) {
    total_draws <- attr(object, "total_draws")
    if (is.null(total_draws)) {
      total_draws <- rowSums(draws)
    }
    draws <- cbind(draws, Total = as.numeric(total_draws))
  }

  if (identical(format, "matrix")) {
    return(draws)
  }

  out <- as.data.frame(draws, check.names = FALSE)
  out <- cbind(.draw = seq_len(nrow(out)), out)
  rownames(out) <- NULL
  out
}

#' @export
popsize.uncounted_bayes <- function(object, by = NULL, level = 0.95,
                                    estimate = c("median", "mean"),
                                    total = FALSE, draw_ids = NULL, ...) {
  estimate <- match.arg(estimate)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number between 0 and 1.", call. = FALSE)
  }

  alpha_eta <- .bayes_alpha_eta_draws(object, draw_ids = draw_ids)
  if (ncol(alpha_eta) != length(object$N)) {
    stop("Posterior alpha draws do not match the fitted data.", call. = FALSE)
  }
  alpha <- if (isTRUE(object$constrained)) stats::plogis(alpha_eta) else alpha_eta

  xi_obs <- exp(sweep(alpha, 2L, log(object$N), `*`))
  group_info <- .popsize_group_info(object, by = by)
  group_idx <- group_info$index

  draws <- matrix(NA_real_, nrow = nrow(xi_obs), ncol = length(group_idx),
                  dimnames = list(NULL, names(group_idx)))
  observed <- numeric(length(group_idx))
  for (i in seq_along(group_idx)) {
    idx <- group_idx[[i]]
    draws[, i] <- rowSums(xi_obs[, idx, drop = FALSE])
    observed[i] <- sum(object$m[idx])
  }

  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  med <- apply(draws, 2L, stats::median, na.rm = TRUE)
  avg <- colMeans(draws, na.rm = TRUE)
  sds <- apply(draws, 2L, stats::sd, na.rm = TRUE)
  ci <- t(apply(draws, 2L, stats::quantile, probs = probs, na.rm = TRUE,
                names = FALSE))
  point <- if (identical(estimate, "median")) med else avg

  ps <- data.frame(
    group = names(group_idx),
    observed = observed,
    estimate = as.numeric(point),
    median = as.numeric(med),
    mean = as.numeric(avg),
    sd = as.numeric(sds),
    lower = as.numeric(ci[, 1]),
    upper = as.numeric(ci[, 2]),
    stringsAsFactors = FALSE
  )
  ps$share_pct <- ps$estimate / sum(ps$estimate) * 100

  if (isTRUE(total)) {
    total_draws <- rowSums(draws)
    total_ci <- stats::quantile(total_draws, probs = probs, na.rm = TRUE,
                                names = FALSE)
    total_med <- stats::median(total_draws, na.rm = TRUE)
    total_mean <- mean(total_draws, na.rm = TRUE)
    attr(ps, "total") <- list(
      estimate = if (identical(estimate, "median")) total_med else total_mean,
      median = total_med,
      mean = total_mean,
      sd = stats::sd(total_draws, na.rm = TRUE),
      lower = total_ci[1],
      upper = total_ci[2],
      n_draws = sum(is.finite(total_draws))
    )
    attr(ps, "total_draws") <- total_draws
  }

  attr(ps, "groups") <- group_info$groups
  attr(ps, "draws") <- draws
  attr(ps, "level") <- level
  attr(ps, "estimate_type") <- estimate
  class(ps) <- c("uncounted_popsize_bayes", "data.frame")
  ps
}

#' @export
print.uncounted_bayes <- function(x, ...) {
  cat("Bayesian hidden-population model\n")
  cat("  engine:", x$engine, "| backend:", x$backend,
      "| method:", x$method, "| link_rho:", x$link_rho, "\n")
  cat("  observations:", x$n_obs,
      "| gamma:", if (isTRUE(x$gamma_estimated)) "estimated" else format(x$gamma),
      "| constrained:", isTRUE(x$constrained), "\n")
  invisible(x)
}

#' @export
print.uncounted_popsize_bayes <- function(x, ...) {
  cat("Bayesian population-size posterior summary\n")
  print.data.frame(x, ...)
  total <- attr(x, "total")
  if (!is.null(total)) {
    cat("\nTotal: estimate =", format(round(total$estimate), big.mark = ","),
        " median =", format(round(total$median), big.mark = ","),
        " mean =", format(round(total$mean), big.mark = ","),
        " CrI [", format(round(total$lower), big.mark = ","),
        ", ", format(round(total$upper), big.mark = ","), "]\n", sep = "")
  }
  invisible(x)
}

#' Summarize a Bayesian Hidden-Population Model
#'
#' @param object An \code{"uncounted_bayes"} object.
#' @param total Logical; include a total population-size summary?
#' @param level Credible interval level.
#' @param diagnostics Logical; include HMC and posterior diagnostic summaries?
#' @param ... Additional arguments, currently ignored.
#'
#' @return Invisibly, a list with coefficient, population-size, and diagnostic
#'   summaries.
#'
#' @export
summary.uncounted_bayes <- function(object, total = FALSE, level = 0.95,
                                    diagnostics = TRUE, ...) {
  coef_tab <- .bayes_coef_table(object, level = level)
  ps <- popsize(object, total = total, level = level)
  diag <- if (isTRUE(diagnostics)) .bayes_diagnostics(object) else NULL

  cat("Bayesian unauthorized population estimation\n")
  cat("Method:", toupper(object$method),
      "| engine:", object$engine,
      "| backend:", object$backend,
      "| link_rho:", object$link_rho, "\n")
  cat("N obs:", object$n_obs, "\n")
  if (isTRUE(object$gamma_estimated)) {
    cat("Gamma: scalar estimated within bounds [",
        format(object$gamma_bounds[1]), ", ",
        format(object$gamma_bounds[2]), "]\n", sep = "")
  } else if (!is.null(object$gamma)) {
    cat("Gamma:", format(object$gamma), "(fixed)\n")
  } else {
    cat("Gamma: none; using auxiliary/reference_pop directly\n")
  }
  if (isTRUE(object$constrained)) {
    cat("Constraints: alpha = inv_logit(eta_alpha), beta = exp(eta_beta)\n")
  }

  if (!is.null(object$fit)) {
    fit_args <- .bayes_fit_args(object)
    if (length(fit_args) > 0L) {
      cat("Chains:", fit_args$chains,
          "| Iter:", fit_args$iter,
          "| Warmup:", fit_args$warmup, "\n")
    }
  }

  cat("\nPosterior coefficients:\n")
  print.data.frame(coef_tab, row.names = FALSE, digits = 4)

  if (!is.null(diag)) {
    cat("\nMCMC diagnostics:\n")
    if (!is.na(diag$divergences)) {
      cat("Divergences:", diag$divergences, "\n")
    }
    if (!is.na(diag$max_treedepth_hits)) {
      cat("Max treedepth hits:", diag$max_treedepth_hits, "\n")
    }
    cat("Max Rhat:", format(diag$max_rhat, digits = 4),
        "| Min bulk ESS:", format(diag$min_bulk_ess, digits = 4),
        "| Min tail ESS:", format(diag$min_tail_ess, digits = 4), "\n")
  }

  cat("\n-----------------------\n")
  cat("Population size posterior summary:\n")
  print.uncounted_popsize_bayes(ps)

  out <- list(
    coefficients = coef_tab,
    popsize = ps,
    diagnostics = diag
  )
  class(out) <- "summary_uncounted_bayes"
  invisible(out)
}

#' @export
coef.uncounted_bayes <- function(object,
                                 summary = c("median", "mean", "draws"),
                                 ...) {
  summary <- match.arg(summary)
  draws <- .bayes_coef_draws(object)
  if (identical(summary, "draws")) {
    return(draws$matrix)
  }
  tab <- .bayes_coef_table(object)
  vals <- tab[[summary]]
  stats::setNames(vals, tab$term)
}

#' @export
fitted.uncounted_bayes <- function(object,
                                   summary = c("mean", "median", "draws"),
                                   ...) {
  summary <- match.arg(summary)
  draws <- .bayes_epred_draws(object, newdata = NULL, ...)
  if (identical(summary, "draws")) {
    return(draws)
  }
  .bayes_summarize_prediction_draws(draws, summary)
}

#' Predict from a Bayesian uncounted Model
#'
#' @param object An \code{"uncounted_bayes"} object.
#' @param newdata Optional data frame for prediction.
#' @param type Prediction scale. \code{"response"} returns posterior expected
#'   observed counts, \code{"link"} returns the log mean, \code{"counts"}
#'   returns posterior predictive replicated counts, and \code{"draws"}
#'   returns response-scale posterior expected-count draws.
#' @param summary Logical or character. Use \code{FALSE} to return draws;
#'   \code{TRUE} or \code{"mean"} to return posterior means; \code{"median"}
#'   to return posterior medians.
#' @param ... Additional arguments passed to the underlying \pkg{brms}
#'   posterior prediction function.
#'
#' @return A numeric vector when summarized, otherwise a draw-by-observation
#'   matrix.
#'
#' @export
predict.uncounted_bayes <- function(object, newdata = NULL,
                                    type = c("response", "link", "counts", "draws"),
                                    summary = TRUE, ...) {
  type <- match.arg(type)
  pred_data <- .bayes_prediction_data(object, newdata)

  draws <- switch(type,
    response = .bayes_epred_draws(object, newdata = pred_data, ...),
    draws = .bayes_epred_draws(object, newdata = pred_data, ...),
    link = .bayes_linpred_draws(object, newdata = pred_data, ...),
    counts = .bayes_predictive_draws(object, newdata = pred_data, ...)
  )

  if (identical(type, "draws") || identical(summary, FALSE)) {
    return(draws)
  }
  summary <- .bayes_match_prediction_summary(summary)
  .bayes_summarize_prediction_draws(draws, summary)
}

#' Bayesian Residuals
#'
#' @param object An \code{"uncounted_bayes"} object.
#' @param type Residual type: response or Pearson.
#' @param summary Posterior summary used for fitted values.
#' @param ... Additional arguments passed to \code{\link{fitted.uncounted_bayes}}.
#'
#' @return A numeric vector of residuals.
#'
#' @export
residuals.uncounted_bayes <- function(object,
                                      type = c("response", "pearson"),
                                      summary = c("mean", "median"),
                                      ...) {
  type <- match.arg(type)
  summary <- match.arg(summary)
  mu <- stats::fitted(object, summary = summary, ...)
  res <- object$m - mu
  if (identical(type, "response")) {
    return(as.numeric(res))
  }

  var_mu <- if (identical(object$method, "nb")) {
    theta <- .bayes_theta_estimate(object, summary = summary)
    if (is.finite(theta) && theta > 0) mu + mu^2 / theta else mu
  } else {
    mu
  }
  as.numeric(res / sqrt(pmax(var_mu, .Machine$double.eps)))
}

#' @export
loo.uncounted_bayes <- function(object, ...) {
  .require_namespace("brms")
  brms::loo(object$fit, ...)
}

#' @export
tidy.uncounted_bayes <- function(x, conf.int = FALSE, conf.level = 0.95,
                                 ...) {
  tab <- .bayes_coef_table(x, level = conf.level)
  out <- data.frame(
    term = tab$term,
    estimate = tab$median,
    std.error = tab$sd,
    statistic = NA_real_,
    p.value = NA_real_,
    stringsAsFactors = FALSE
  )
  if (conf.int) {
    out$conf.low <- tab$conf.low
    out$conf.high <- tab$conf.high
  }
  out
}

#' @export
glance.uncounted_bayes <- function(x, loo = FALSE, ...) {
  diag <- .bayes_diagnostics(x)
  out <- data.frame(
    method = toupper(x$method),
    estimator = "BAYES",
    engine = x$engine,
    backend = x$backend,
    link_rho = x$link_rho,
    nobs = x$n_obs,
    gamma = if (is.null(x$gamma)) NA_real_ else x$gamma,
    gamma_estimated = isTRUE(x$gamma_estimated),
    divergences = diag$divergences,
    max_treedepth_hits = diag$max_treedepth_hits,
    max_rhat = diag$max_rhat,
    min_bulk_ess = diag$min_bulk_ess,
    min_tail_ess = diag$min_tail_ess,
    stringsAsFactors = FALSE
  )
  if (isTRUE(loo)) {
    loo_obj <- tryCatch(loo.uncounted_bayes(x, ...), error = function(e) NULL)
    out$elpd_loo <- if (!is.null(loo_obj) && "estimates" %in% names(loo_obj)) {
      loo_obj$estimates["elpd_loo", "Estimate"]
    } else NA_real_
    out$p_loo <- if (!is.null(loo_obj) && "estimates" %in% names(loo_obj)) {
      loo_obj$estimates["p_loo", "Estimate"]
    } else NA_real_
    out$looic <- if (!is.null(loo_obj) && "estimates" %in% names(loo_obj)) {
      loo_obj$estimates["looic", "Estimate"]
    } else NA_real_
  }
  out
}

.bayes_make_brms_spec <- function(data,
                                  observed,
                                  auxiliary,
                                  reference_pop,
                                  method = c("poisson", "nb"),
                                  cov_alpha = NULL,
                                  cov_beta = NULL,
                                  gamma = "estimate",
                                  gamma_bounds = c(1e-10, 0.5),
                                  link_rho = c("power", "cloglog", "logit", "probit"),
                                  constrained = FALSE,
                                  weights = NULL,
                                  prior = NULL) {
  .require_namespace("brms")
  method <- match.arg(method)
  link_rho <- .normalize_link_rho(link_rho)

  m <- eval(observed[[2]], data)
  n_aux <- eval(auxiliary[[2]], data)
  N <- eval(reference_pop[[2]], data)

  if (length(m) != length(N) || length(m) != length(n_aux)) {
    stop("'observed', 'auxiliary', and 'reference_pop' must have equal lengths.",
         call. = FALSE)
  }
  if (any(!is.finite(m)) || any(m < 0) ||
      any(abs(m - round(m)) > sqrt(.Machine$double.eps))) {
    stop("'observed' must evaluate to non-negative whole counts.",
         call. = FALSE)
  }
  if (any(!is.finite(n_aux)) || any(n_aux < 0)) {
    stop("'auxiliary' must evaluate to non-negative finite values.",
         call. = FALSE)
  }
  if (any(!is.finite(N)) || any(N <= 0)) {
    stop("'reference_pop' must evaluate to positive finite values.",
         call. = FALSE)
  }

  ratio <- n_aux / N
  spec_data <- as.data.frame(data)
  spec_data$.uncounted_observed <- as.integer(round(m))
  spec_data$.uncounted_log_N <- log(N)

  estimate_gamma <- identical(gamma, "estimate")
  gamma_fixed <- is.numeric(gamma) && length(gamma) == 1L && is.finite(gamma)
  no_gamma <- is.null(gamma)
  gamma_value <- NULL

  if (estimate_gamma) {
    .validate_gamma_bounds(gamma_bounds)
    spec_data$.uncounted_ratio <- ratio
    spec_data$.uncounted_gamma_lb <- gamma_bounds[1]
    spec_data$.uncounted_gamma_ub <- gamma_bounds[2]
  } else if (gamma_fixed) {
    if (gamma < 0) {
      stop("'gamma' must be non-negative.", call. = FALSE)
    }
    rate <- gamma + ratio
    if (any(rate <= 0)) {
      stop("'gamma + auxiliary/reference_pop' must be strictly positive.",
           call. = FALSE)
    }
    spec_data$.uncounted_log_rate <- log(rate)
    gamma_value <- gamma
  } else if (no_gamma) {
    if (any(ratio <= 0)) {
      stop("When 'gamma = NULL', all auxiliary/reference_pop ratios must be ",
           "strictly positive.", call. = FALSE)
    }
    spec_data$.uncounted_log_rate <- log(ratio)
  } else {
    stop("'gamma' must be numeric, 'estimate', or NULL.", call. = FALSE)
  }

  obs_weights <- .bayes_eval_weights(weights, data)
  lhs <- ".uncounted_observed"
  if (!is.null(obs_weights)) {
    spec_data$.uncounted_weights <- obs_weights
    lhs <- ".uncounted_observed | weights(.uncounted_weights)"
  }

  alpha_expr <- if (isTRUE(constrained)) "inv_logit(alpha)" else "alpha"
  beta_expr <- if (isTRUE(constrained)) "exp(beta)" else "beta"
  eta_expr <- if (estimate_gamma) {
    paste0(beta_expr,
           " * log(.uncounted_gamma_lb + ",
           "(.uncounted_gamma_ub - .uncounted_gamma_lb) * ",
           "inv_logit(gamma) + .uncounted_ratio)")
  } else {
    paste0(beta_expr, " * .uncounted_log_rate")
  }
  log_rho_expr <- switch(link_rho,
    power = eta_expr,
    logit = paste0("uncounted_log_logit(", eta_expr, ")"),
    probit = paste0("uncounted_log_probit(", eta_expr, ")"),
    cloglog = paste0("uncounted_log_cloglog(", eta_expr, ")")
  )

  main_formula <- stats::as.formula(
    paste(lhs, "~", alpha_expr, "* .uncounted_log_N +", log_rho_expr)
  )
  bf_args <- list(
    formula = main_formula,
    alpha = if (is.null(cov_alpha)) ~ 1 else cov_alpha,
    beta = if (is.null(cov_beta)) ~ 1 else cov_beta,
    nl = TRUE
  )
  if (estimate_gamma) {
    bf_args$gamma <- ~ 1
  }
  brms_formula <- do.call(brms::bf, bf_args)

  brms_family <- switch(method,
    poisson = brms::brmsfamily("poisson", link = "log"),
    nb = brms::brmsfamily("negbinomial", link = "log")
  )
  brms_prior <- if (is.null(prior)) {
    default_priors_uncounted_bayes(method = method, constrained = constrained,
                                   gamma = gamma, gamma_bounds = gamma_bounds)
  } else {
    prior
  }

  X_alpha <- .build_model_matrix(cov_alpha, data)
  X_beta <- .build_model_matrix(cov_beta, data)
  if (!is.null(cov_alpha)) {
    colnames(X_alpha) <- paste0("alpha:", colnames(X_alpha))
  } else {
    colnames(X_alpha) <- "alpha"
  }
  if (!is.null(cov_beta)) {
    colnames(X_beta) <- paste0("beta:", colnames(X_beta))
  } else {
    colnames(X_beta) <- "beta"
  }

  list(
    formula = brms_formula,
    family = brms_family,
    data = spec_data,
    prior = brms_prior,
    stanvars = .bayes_stanvars(link_rho),
    m = as.integer(round(m)),
    n_aux = n_aux,
    N = N,
    ratio = ratio,
    gamma = gamma_value,
    gamma_estimated = estimate_gamma,
    obs_weights = obs_weights,
    X_alpha = X_alpha,
    X_beta = X_beta,
    cov_alpha_vars = if (!is.null(cov_alpha)) {
      data[, all.vars(cov_alpha), drop = FALSE]
    } else {
      NULL
    }
  )
}

.bayes_stanvars <- function(link_rho) {
  if (identical(link_rho, "power")) {
    return(NULL)
  }
  .require_namespace("brms")
  code <- "
real uncounted_log_logit(real eta) {
  return log_inv_logit(eta);
}
real uncounted_log_probit(real eta) {
  return std_normal_lcdf(eta);
}
real uncounted_log_cloglog(real eta) {
  if (eta > 20) return 0;
  return log1m_exp(-exp(eta));
}
"
  brms::stanvar(scode = code, block = "functions")
}

.bayes_eval_weights <- function(weights, data) {
  if (is.null(weights)) {
    return(NULL)
  }
  w <- if (inherits(weights, "formula")) {
    eval(weights[[2]], data)
  } else {
    weights
  }
  if (!is.numeric(w) || length(w) != nrow(data) || any(!is.finite(w)) ||
      any(w < 0)) {
    stop("'weights' must be a non-negative finite numeric vector with one ",
         "entry per row of 'data'.", call. = FALSE)
  }
  as.numeric(w)
}

.bayes_alpha_eta_draws <- function(object, draw_ids = NULL) {
  if (!is.null(object$alpha_draws)) {
    out <- as.matrix(object$alpha_draws)
    if (!is.null(draw_ids)) {
      out <- out[draw_ids, , drop = FALSE]
    }
    return(out)
  }
  .require_namespace("brms")
  out <- brms::posterior_linpred(
    object$fit,
    nlpar = "alpha",
    transform = FALSE,
    draw_ids = draw_ids
  )
  if (length(dim(out)) == 3L) {
    out <- out[, , 1L, drop = TRUE]
  }
  as.matrix(out)
}

.bayes_coef_draws <- function(object) {
  .require_namespace("posterior")
  draws <- if (!is.null(object$coef_draws)) {
    posterior::as_draws_df(as.matrix(object$coef_draws))
  } else {
    .require_namespace("brms")
    posterior::as_draws_df(object$fit)
  }
  vars <- posterior::variables(draws)
  vars <- grep("^(b_alpha_|b_beta_|b_gamma_|shape$)", vars, value = TRUE)
  if (length(vars) == 0L) {
    stop("No alpha, beta, gamma, or negative-binomial shape posterior draws ",
         "were found.", call. = FALSE)
  }
  mat <- as.matrix(as.data.frame(draws)[vars])
  colnames(mat) <- vars
  list(draws = draws, vars = vars, matrix = mat)
}

.bayes_coef_table <- function(object, level = 0.95) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number between 0 and 1.", call. = FALSE)
  }
  cd <- .bayes_coef_draws(object)
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  mat <- cd$matrix

  tab <- data.frame(
    variable = colnames(mat),
    term = vapply(colnames(mat), .bayes_clean_term, character(1)),
    mean = colMeans(mat, na.rm = TRUE),
    median = apply(mat, 2L, stats::median, na.rm = TRUE),
    sd = apply(mat, 2L, stats::sd, na.rm = TRUE),
    conf.low = apply(mat, 2L, stats::quantile, probs = probs[1],
                     na.rm = TRUE, names = FALSE),
    conf.high = apply(mat, 2L, stats::quantile, probs = probs[2],
                      na.rm = TRUE, names = FALSE),
    stringsAsFactors = FALSE
  )

  diag <- tryCatch({
    s <- posterior::summarise_draws(
      posterior::subset_draws(cd$draws, variable = cd$vars),
      "rhat", "ess_bulk", "ess_tail"
    )
    as.data.frame(s)
  }, error = function(e) NULL)

  if (!is.null(diag)) {
    tab$rhat <- diag$rhat[match(tab$variable, diag$variable)]
    tab$ess_bulk <- diag$ess_bulk[match(tab$variable, diag$variable)]
    tab$ess_tail <- diag$ess_tail[match(tab$variable, diag$variable)]
  } else {
    tab$rhat <- NA_real_
    tab$ess_bulk <- NA_real_
    tab$ess_tail <- NA_real_
  }

  tab <- tab[, c("term", "mean", "median", "sd", "conf.low", "conf.high",
                 "rhat", "ess_bulk", "ess_tail"), drop = FALSE]
  rownames(tab) <- NULL
  tab
}

.bayes_clean_term <- function(variable) {
  if (identical(variable, "shape")) {
    return("theta")
  }
  term <- sub("^b_", "", variable)
  term <- sub("^alpha_", "alpha:", term)
  term <- sub("^beta_", "beta:", term)
  term <- sub("^gamma_Intercept$", "gamma", term)
  term <- sub("^gamma_", "gamma:", term)
  term <- sub(":Intercept$", ":(Intercept)", term)
  term
}

.bayes_diagnostics <- function(object) {
  coef_tab <- tryCatch(.bayes_coef_table(object), error = function(e) NULL)
  max_rhat <- if (!is.null(coef_tab)) {
    suppressWarnings(max(coef_tab$rhat, na.rm = TRUE))
  } else NA_real_
  min_bulk <- if (!is.null(coef_tab)) {
    suppressWarnings(min(coef_tab$ess_bulk, na.rm = TRUE))
  } else NA_real_
  min_tail <- if (!is.null(coef_tab)) {
    suppressWarnings(min(coef_tab$ess_tail, na.rm = TRUE))
  } else NA_real_
  if (!is.finite(max_rhat)) max_rhat <- NA_real_
  if (!is.finite(min_bulk)) min_bulk <- NA_real_
  if (!is.finite(min_tail)) min_tail <- NA_real_

  divergences <- max_treedepth_hits <- NA_integer_
  if (!is.null(object$fit) && requireNamespace("brms", quietly = TRUE)) {
    np <- tryCatch(brms::nuts_params(object$fit), error = function(e) NULL)
    if (!is.null(np)) {
      divergences <- sum(np$Parameter == "divergent__" & np$Value == 1,
                         na.rm = TRUE)
      max_treedepth_hits <- sum(np$Parameter == "treedepth__" &
                                  np$Value >= 10, na.rm = TRUE)
    }
  }

  list(
    divergences = divergences,
    max_treedepth_hits = max_treedepth_hits,
    max_rhat = max_rhat,
    min_bulk_ess = min_bulk,
    min_tail_ess = min_tail
  )
}

.bayes_fit_args <- function(object) {
  out <- tryCatch({
    sim <- object$fit$fit@sim
    list(
      chains = sim$chains,
      iter = sim$iter,
      warmup = if (!is.null(sim$warmup2)) sim$warmup2[1] else NA_integer_
    )
  }, error = function(e) list())
  out
}

.bayes_epred_draws <- function(object, newdata = NULL, ...) {
  if (is.null(newdata) && !is.null(object$fitted_draws)) {
    return(as.matrix(object$fitted_draws))
  }
  .require_namespace("brms")
  as.matrix(brms::posterior_epred(object$fit, newdata = newdata, ...))
}

.bayes_linpred_draws <- function(object, newdata = NULL, ...) {
  .require_namespace("brms")
  as.matrix(brms::posterior_linpred(object$fit, newdata = newdata,
                                    transform = FALSE, ...))
}

.bayes_predictive_draws <- function(object, newdata = NULL, ...) {
  .require_namespace("brms")
  as.matrix(brms::posterior_predict(object$fit, newdata = newdata, ...))
}

.bayes_prediction_data <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    return(NULL)
  }
  if (is.null(object$call)) {
    stop("This Bayesian fit does not contain the original call, so newdata ",
         "predictions cannot be prepared.", call. = FALSE)
  }

  out <- as.data.frame(newdata)
  refpop_var <- all.vars(object$call$reference_pop)
  auxiliary_var <- all.vars(object$call$auxiliary)
  observed_var <- all.vars(object$call$observed)
  if (!all(c(refpop_var, auxiliary_var) %in% names(out))) {
    stop("'newdata' must contain the auxiliary and reference-population ",
         "variables used to fit the model.", call. = FALSE)
  }
  N <- out[[refpop_var]]
  n_aux <- out[[auxiliary_var]]
  ratio <- n_aux / N
  out$.uncounted_observed <- if (observed_var %in% names(out)) {
    out[[observed_var]]
  } else {
    0L
  }
  out$.uncounted_log_N <- log(N)

  if (isTRUE(object$gamma_estimated)) {
    out$.uncounted_ratio <- ratio
    out$.uncounted_gamma_lb <- object$gamma_bounds[1]
    out$.uncounted_gamma_ub <- object$gamma_bounds[2]
  } else {
    rate <- .rate_from_gamma(ratio, if (!is.null(object$gamma)) object$gamma else NULL)
    if (any(rate <= 0)) {
      stop("Prediction data imply non-positive detection-rate terms.",
           call. = FALSE)
    }
    out$.uncounted_log_rate <- log(rate)
  }
  if (!is.null(object$obs_weights)) {
    out$.uncounted_weights <- 1
  }
  out
}

.bayes_match_prediction_summary <- function(summary) {
  if (identical(summary, TRUE)) {
    return("mean")
  }
  if (!is.character(summary) || length(summary) != 1L) {
    stop("'summary' must be TRUE, FALSE, 'mean', or 'median'.", call. = FALSE)
  }
  match.arg(summary, c("mean", "median"))
}

.bayes_summarize_prediction_draws <- function(draws, summary) {
  if (identical(summary, "mean")) {
    return(as.numeric(colMeans(draws, na.rm = TRUE)))
  }
  as.numeric(apply(draws, 2L, stats::median, na.rm = TRUE))
}

.bayes_theta_estimate <- function(object, summary = c("mean", "median")) {
  summary <- match.arg(summary)
  tab <- tryCatch(.bayes_coef_table(object), error = function(e) NULL)
  if (is.null(tab) || !("theta" %in% tab$term)) {
    return(NA_real_)
  }
  tab[[summary]][match("theta", tab$term)]
}

.validate_gamma_bounds <- function(gamma_bounds) {
  if (!is.numeric(gamma_bounds) || length(gamma_bounds) != 2L ||
      any(!is.finite(gamma_bounds)) || gamma_bounds[1] < 0 ||
      gamma_bounds[1] >= gamma_bounds[2]) {
    stop("'gamma_bounds' must be two finite numbers with 0 <= lower < upper.",
         call. = FALSE)
  }
  invisible(gamma_bounds)
}

.require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for this operation. ",
         "Install it and try again.", call. = FALSE)
  }
  invisible(TRUE)
}
