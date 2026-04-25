#' Hypothesis Tests for Population Size Estimates
#'
#' Tests scalar hypotheses and contrasts involving population size estimates
#' returned by \code{\link{popsize}} or \code{\link{bootstrap_popsize}}.
#'
#' @param object An \code{"uncounted_popsize"} object from \code{\link{popsize}}
#'   or an \code{"uncounted_boot"} object from \code{\link{bootstrap_popsize}}.
#' @param hypothesis Numeric null value(s), a character hypothesis, or a
#'   function. Character hypotheses can refer to population size estimates with
#'   \code{xi[...]} filters, such as \code{"xi[year == 2024] > 200000"} or
#'   \code{"xi[year == 2024] - xi[year == 2019] > 0"}. Positional aliases
#'   \code{b1}, \code{b2}, ... are also available.
#' @param estimate Required character string naming the estimate to test.
#'   For \code{popsize()} objects, use \code{"estimate"} or
#'   \code{"estimate_bc"}. For \code{bootstrap_popsize()} objects, use
#'   \code{"plugin"}, \code{"boot_median"}, or \code{"boot_mean"}.
#' @param level Confidence level for Wald intervals on the tested contrast.
#' @param df Degrees of freedom for p values and confidence intervals. The
#'   default \code{Inf} uses the standard normal distribution.
#' @param hypothesis_side Character string. \code{"alternative"} (default)
#'   treats directional character hypotheses such as \code{"xi[...] < 15000"}
#'   as the alternative hypothesis to support. \code{"null"} treats the
#'   directional character hypothesis as the null hypothesis to test against.
#'   Equality hypotheses and numeric/function hypotheses are always interpreted
#'   as null-equality tests.
#' @param include_total Logical; if \code{TRUE}, append a \code{"Total"} xi
#'   target when multiple groups are available.
#' @param ... Additional arguments, currently ignored.
#'
#' @details
#' For analytical \code{popsize()} results, tests use Wald inference with
#' delta-method standard errors for nonlinear functions of the fitted alpha
#' coefficients. For bootstrap results, the same contrast is evaluated in each
#' fractional weighted bootstrap draw, and the standard error is the bootstrap
#' standard deviation of those contrasts. The p value is a Wald p value, not an
#' empirical exceedance probability; use \code{\link{exceedance_popsize}} for
#' bootstrap tail areas.
#'
#' @return A data frame with one row per hypothesis and columns
#'   \code{hypothesis}, \code{null_hypothesis},
#'   \code{alternative_hypothesis}, \code{estimate}, \code{null.value},
#'   \code{contrast}, \code{std.error}, \code{statistic}, \code{p.value},
#'   \code{s.value}, \code{conf.low}, \code{conf.high}, \code{alternative},
#'   \code{method}, \code{estimate_type}, and \code{n_draws}.
#'
#' @examples
#' data(irregular_migration)
#' d <- irregular_migration[irregular_migration$year %in% c("2019", "2024"), ]
#' fit <- estimate_hidden_pop(
#'   data = d, observed = ~m, auxiliary = ~n, reference_pop = ~N,
#'   method = "poisson", gamma = 0.005
#' )
#' ps_year <- popsize(fit, by = ~year)
#' hypotheses_popsize(ps_year, "xi[year == 2024] > 200000",
#'                    estimate = "estimate_bc")
#' hypotheses_popsize(ps_year, "xi[year == 2024] - xi[year == 2019] > 0",
#'                    estimate = "estimate_bc")
#'
#' @export
hypotheses_popsize <- function(object, hypothesis, estimate, level = 0.95,
                               df = Inf,
                               hypothesis_side = c("alternative", "null"),
                               include_total = FALSE, ...) {
  if (missing(estimate)) {
    stop("'estimate' must be supplied explicitly.", call. = FALSE)
  }
  hypothesis_side <- match.arg(hypothesis_side)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(df) || length(df) != 1L || is.na(df) || df <= 0) {
    stop("'df' must be a single positive number or Inf.", call. = FALSE)
  }

  source <- .pop_hyp_source(object, estimate = estimate,
                            include_total = include_total)

  out <- if (is.numeric(hypothesis)) {
    .pop_hyp_numeric(source, hypothesis, level = level, df = df)
  } else if (is.character(hypothesis)) {
    .pop_hyp_character(source, hypothesis, level = level, df = df,
                       hypothesis_side = hypothesis_side)
  } else if (is.function(hypothesis)) {
    .pop_hyp_function(source, hypothesis, level = level, df = df)
  } else {
    stop("'hypothesis' must be numeric, character, or a function.",
         call. = FALSE)
  }

  rownames(out) <- NULL
  class(out) <- c("uncounted_popsize_hypotheses", "data.frame")
  out
}

#' @export
print.uncounted_popsize_hypotheses <- function(x, ...) {
  cat("Population-size hypothesis tests\n")
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}

.pop_hyp_source <- function(object, estimate, include_total = FALSE) {
  if (inherits(object, "uncounted_popsize")) {
    .pop_hyp_source_popsize(object, estimate, include_total)
  } else if (inherits(object, "uncounted_boot")) {
    .pop_hyp_source_boot(object, estimate, include_total)
  } else {
    stop("'object' must be an 'uncounted_popsize' or 'uncounted_boot' object.",
         call. = FALSE)
  }
}

.pop_hyp_source_popsize <- function(object, estimate, include_total) {
  valid <- c("estimate", "estimate_bc")
  if (!(estimate %in% valid)) {
    stop("'estimate' must be one of: ", paste(valid, collapse = ", "),
         call. = FALSE)
  }
  x <- object[[estimate]]
  if (is.null(x) || any(!is.finite(x))) {
    stop("Selected estimate contains missing or non-finite values.", call. = FALSE)
  }
  groups <- .pop_hyp_groups_from_popsize(object)
  names(x) <- groups$.group

  V <- attr(object, paste0("vcov_", estimate))
  if (is.null(V)) {
    stop("This popsize object does not contain covariance metadata. ",
         "Recompute it with popsize() before calling hypotheses_popsize().",
         call. = FALSE)
  }
  V <- as.matrix(V)
  dimnames(V) <- list(groups$.group, groups$.group)

  if (include_total && length(x) > 1L) {
    total <- .pop_hyp_append_total(x, V, groups)
    x <- total$x
    V <- total$V
    groups <- total$groups
  }

  list(type = "popsize", x = x, V = V, draws = NULL, groups = groups,
       method = "Delta method Wald test", estimate_type = estimate)
}

.pop_hyp_source_boot <- function(object, estimate, include_total) {
  valid <- c("plugin", "boot_median", "boot_mean")
  if (!(estimate %in% valid)) {
    stop("'estimate' must be one of: ", paste(valid, collapse = ", "),
         call. = FALSE)
  }
  if (is.null(object$popsize_full) || !(estimate %in% names(object$popsize_full))) {
    stop("Bootstrap object does not contain the selected estimate.",
         call. = FALSE)
  }

  non_total <- object$popsize_full$group != "Total"
  ps_full <- object$popsize_full[non_total, , drop = FALSE]
  x <- ps_full[[estimate]]
  groups <- .pop_hyp_groups_from_boot(object, ps_full$group)
  names(x) <- groups$.group

  draws <- as.matrix(object$t)
  if (ncol(draws) != length(x)) {
    stop("Bootstrap draw matrix does not match the population-size table.",
         call. = FALSE)
  }
  colnames(draws) <- groups$.group

  if (include_total && length(x) > 1L) {
    total_draws <- rowSums(draws)
    total_value <- .pop_hyp_boot_total_value(object, estimate, x, total_draws)
    x <- c(x, Total = total_value)
    draws <- cbind(draws, Total = total_draws)
    groups <- .pop_hyp_append_total_group(groups)
  }

  list(type = "boot", x = x, V = NULL, draws = draws, groups = groups,
       method = "FWB Wald test", estimate_type = estimate)
}

.pop_hyp_numeric <- function(source, null, level, df) {
  x <- source$x
  null <- .pop_hyp_expand_null(null, x)
  labels <- paste0("xi[", names(x), "] = ", format(null, trim = TRUE))

  if (identical(source$type, "popsize")) {
    se <- sqrt(pmax(diag(source$V), 0))
    contrast <- x - null
    .pop_hyp_result(labels, estimate = x, null = null, contrast = contrast,
                    se = se, alternative = rep("two.sided", length(x)),
                    null_hypothesis = paste0("xi[", names(x), "] = ",
                                             format(null, trim = TRUE)),
                    alternative_hypothesis = paste0("xi[", names(x), "] != ",
                                                    format(null, trim = TRUE)),
                    method = source$method,
                    estimate_type = source$estimate_type,
                    n_draws = rep(NA_integer_, length(x)),
                    level = level, df = df)
  } else {
    rows <- lapply(seq_along(x), function(i) {
      contrast_draws <- source$draws[, i] - null[i]
      finite <- is.finite(contrast_draws)
      se <- stats::sd(contrast_draws[finite])
      .pop_hyp_result(labels[i], estimate = x[i], null = null[i],
                      contrast = x[i] - null[i], se = se,
                      alternative = "two.sided",
                      null_hypothesis = paste0("xi[", names(x)[i], "] = ",
                                               format(null[i], trim = TRUE)),
                      alternative_hypothesis = paste0("xi[", names(x)[i],
                                                      "] != ",
                                                      format(null[i], trim = TRUE)),
                      method = source$method,
                      estimate_type = source$estimate_type,
                      n_draws = sum(finite), level = level, df = df)
    })
    do.call(rbind, rows)
  }
}

.pop_hyp_character <- function(source, hypothesis, level, df,
                               hypothesis_side) {
  rows <- lapply(hypothesis, function(h) {
    parsed <- .pop_hyp_parse_comparison(h, hypothesis_side = hypothesis_side)

    if (identical(source$type, "popsize")) {
      lhs <- .pop_hyp_eval_value(parsed$lhs, source$x, source$groups)
      rhs <- .pop_hyp_eval_value(parsed$rhs, source$x, source$groups)
      contrast <- lhs - rhs
      se <- sqrt(max(0, as.numeric(t(contrast$weights) %*%
                                     source$V %*% contrast$weights)))
      .pop_hyp_result(h, estimate = lhs$value, null = rhs$value,
                      contrast = contrast$value, se = se,
                      alternative = parsed$alternative,
                      null_hypothesis = parsed$null_hypothesis,
                      alternative_hypothesis = parsed$alternative_hypothesis,
                      method = source$method,
                      estimate_type = source$estimate_type,
                      n_draws = NA_integer_, level = level, df = df)
    } else {
      lhs0 <- .pop_hyp_eval_numeric(parsed$lhs, source$x, source$groups)
      rhs0 <- .pop_hyp_eval_numeric(parsed$rhs, source$x, source$groups)
      contrast0 <- lhs0 - rhs0
      contrast_draws <- vapply(seq_len(nrow(source$draws)), function(i) {
        row <- source$draws[i, ]
        .pop_hyp_eval_numeric(parsed$lhs, row, source$groups) -
          .pop_hyp_eval_numeric(parsed$rhs, row, source$groups)
      }, numeric(1))
      finite <- is.finite(contrast_draws)
      se <- stats::sd(contrast_draws[finite])
      .pop_hyp_result(h, estimate = lhs0, null = rhs0, contrast = contrast0,
                      se = se, alternative = parsed$alternative,
                      null_hypothesis = parsed$null_hypothesis,
                      alternative_hypothesis = parsed$alternative_hypothesis,
                      method = source$method,
                      estimate_type = source$estimate_type,
                      n_draws = sum(finite), level = level, df = df)
    }
  })
  do.call(rbind, rows)
}

.pop_hyp_function <- function(source, hypothesis, level, df) {
  f0 <- .pop_hyp_call_fun(hypothesis, source$x, source$groups)
  labels <- names(f0)
  if (is.null(labels) || any(labels == "")) {
    labels <- paste0("function_", seq_along(f0))
  }

  if (identical(source$type, "popsize")) {
    J <- .pop_hyp_jacobian(hypothesis, source$x, source$groups, f0)
    V <- J %*% source$V %*% t(J)
    se <- sqrt(pmax(diag(V), 0))
    .pop_hyp_result(labels, estimate = f0, null = rep(0, length(f0)),
                    contrast = f0, se = se,
                    alternative = rep("two.sided", length(f0)),
                    null_hypothesis = paste0(labels, " = 0"),
                    alternative_hypothesis = paste0(labels, " != 0"),
                    method = source$method,
                    estimate_type = source$estimate_type,
                    n_draws = rep(NA_integer_, length(f0)),
                    level = level, df = df)
  } else {
    vals <- lapply(seq_len(nrow(source$draws)), function(i) {
      .pop_hyp_call_fun(hypothesis, source$draws[i, ], source$groups)
    })
    mat <- do.call(rbind, vals)
    if (ncol(mat) != length(f0)) {
      stop("Function hypothesis must return the same length for all draws.",
           call. = FALSE)
    }
    se <- apply(mat, 2, function(col) stats::sd(col[is.finite(col)]))
    n_draws <- apply(mat, 2, function(col) sum(is.finite(col)))
    .pop_hyp_result(labels, estimate = f0, null = rep(0, length(f0)),
                    contrast = f0, se = se,
                    alternative = rep("two.sided", length(f0)),
                    null_hypothesis = paste0(labels, " = 0"),
                    alternative_hypothesis = paste0(labels, " != 0"),
                    method = source$method,
                    estimate_type = source$estimate_type,
                    n_draws = n_draws, level = level, df = df)
  }
}

.pop_hyp_result <- function(hypothesis, estimate, null, contrast, se,
                            alternative, null_hypothesis,
                            alternative_hypothesis, method, estimate_type,
                            n_draws,
                            level, df) {
  statistic <- contrast / se
  statistic[!is.finite(statistic) & contrast == 0] <- NA_real_

  p <- .pop_hyp_p_value(statistic, alternative, df)
  s <- -log2(p)
  crit <- if (is.finite(df)) {
    stats::qt((1 + level) / 2, df = df)
  } else {
    stats::qnorm((1 + level) / 2)
  }

  data.frame(
    hypothesis = hypothesis,
    null_hypothesis = null_hypothesis,
    alternative_hypothesis = alternative_hypothesis,
    estimate = as.numeric(estimate),
    null.value = as.numeric(null),
    contrast = as.numeric(contrast),
    std.error = as.numeric(se),
    statistic = as.numeric(statistic),
    p.value = as.numeric(p),
    s.value = as.numeric(s),
    conf.low = as.numeric(contrast - crit * se),
    conf.high = as.numeric(contrast + crit * se),
    alternative = alternative,
    method = method,
    estimate_type = estimate_type,
    n_draws = as.integer(n_draws),
    stringsAsFactors = FALSE
  )
}

.pop_hyp_p_value <- function(statistic, alternative, df) {
  p <- numeric(length(statistic))
  for (i in seq_along(statistic)) {
    z <- statistic[i]
    if (!is.finite(z)) {
      p[i] <- if (is.infinite(z)) {
        if (alternative[i] == "greater") {
          if (z > 0) 0 else 1
        } else if (alternative[i] == "less") {
          if (z < 0) 0 else 1
        } else {
          0
        }
      } else {
        NA_real_
      }
      next
    }
    if (is.finite(df)) {
      p[i] <- switch(alternative[i],
        greater = stats::pt(z, df = df, lower.tail = FALSE),
        less = stats::pt(z, df = df),
        two.sided = 2 * stats::pt(abs(z), df = df, lower.tail = FALSE)
      )
    } else {
      p[i] <- switch(alternative[i],
        greater = stats::pnorm(z, lower.tail = FALSE),
        less = stats::pnorm(z),
        two.sided = 2 * stats::pnorm(abs(z), lower.tail = FALSE)
      )
    }
  }
  p
}

.pop_hyp_parse_comparison <- function(hypothesis,
                                      hypothesis_side = c("alternative", "null")) {
  hypothesis_side <- match.arg(hypothesis_side)
  op <- .pop_hyp_find_operator(hypothesis)
  if (is.null(op)) {
    stop("Character hypotheses must include one of =, ==, >, >=, <, or <=.",
         call. = FALSE)
  }
  lhs <- trimws(substr(hypothesis, 1, op$start - 1))
  rhs <- trimws(substr(hypothesis, op$end + 1, nchar(hypothesis)))
  if (!nzchar(lhs) || !nzchar(rhs)) {
    stop("Both sides of a character hypothesis must be non-empty.",
         call. = FALSE)
  }
  hypotheses <- .pop_hyp_h0_h1(lhs, rhs, op$operator, hypothesis_side)
  alternative <- hypotheses$alternative
  list(lhs = lhs, rhs = rhs, operator = op$operator,
       alternative = alternative,
       null_hypothesis = hypotheses$null_hypothesis,
       alternative_hypothesis = hypotheses$alternative_hypothesis)
}

.pop_hyp_h0_h1 <- function(lhs, rhs, operator, hypothesis_side) {
  if (operator %in% c("=", "==")) {
    return(list(
      alternative = "two.sided",
      null_hypothesis = paste(lhs, "=", rhs),
      alternative_hypothesis = paste(lhs, "!=", rhs)
    ))
  }

  lhs_greater <- operator %in% c(">", ">=")
  if (identical(hypothesis_side, "alternative")) {
    alternative <- if (lhs_greater) "greater" else "less"
    null_op <- if (lhs_greater) "<=" else ">="
    alt_op <- if (lhs_greater) ">" else "<"
  } else {
    alternative <- if (lhs_greater) "less" else "greater"
    null_op <- if (lhs_greater) ">=" else "<="
    alt_op <- if (lhs_greater) "<" else ">"
  }

  list(
    alternative = alternative,
    null_hypothesis = paste(lhs, null_op, rhs),
    alternative_hypothesis = paste(lhs, alt_op, rhs)
  )
}

.pop_hyp_find_operator <- function(x) {
  chars <- strsplit(x, "", fixed = TRUE)[[1]]
  depth <- 0L
  quote <- NULL
  escape <- FALSE
  n <- length(chars)

  for (i in seq_len(n)) {
    ch <- chars[i]
    if (!is.null(quote)) {
      if (escape) {
        escape <- FALSE
      } else if (ch == "\\") {
        escape <- TRUE
      } else if (ch == quote) {
        quote <- NULL
      }
      next
    }
    if (ch %in% c("'", "\"")) {
      quote <- ch
      next
    }
    if (ch %in% c("(", "[", "{")) {
      depth <- depth + 1L
      next
    }
    if (ch %in% c(")", "]", "}")) {
      depth <- max(0L, depth - 1L)
      next
    }
    if (depth == 0L) {
      two <- if (i < n) paste0(ch, chars[i + 1L]) else ch
      if (two %in% c(">=", "<=", "==")) {
        return(list(operator = two, start = i, end = i + 1L))
      }
      if (ch %in% c(">", "<", "=")) {
        return(list(operator = ch, start = i, end = i))
      }
    }
  }
  NULL
}

.pop_hyp_eval_value <- function(expr, x, groups) {
  n <- length(x)
  env <- new.env(parent = baseenv())
  env$.xi <- function(filter) {
    idx <- .pop_hyp_select_groups(filter, groups)
    w <- numeric(n)
    w[idx] <- 1
    .pop_hyp_value(sum(x[idx]), w)
  }
  for (i in seq_len(n)) {
    w <- numeric(n)
    w[i] <- 1
    assign(paste0("b", i), .pop_hyp_value(x[i], w), envir = env)
  }
  out <- eval(parse(text = .pop_hyp_rewrite_xi(expr)), envir = env)
  .pop_hyp_as_value(out, n)
}

.pop_hyp_eval_numeric <- function(expr, x, groups) {
  env <- new.env(parent = baseenv())
  env$.xi <- function(filter) {
    idx <- .pop_hyp_select_groups(filter, groups)
    sum(x[idx])
  }
  for (i in seq_along(x)) {
    assign(paste0("b", i), x[i], envir = env)
  }
  as.numeric(eval(parse(text = .pop_hyp_rewrite_xi(expr)), envir = env))
}

.pop_hyp_rewrite_xi <- function(expr) {
  expr <- gsub("\\bxi\\s+\\[", "xi[", expr)
  chars <- strsplit(expr, "", fixed = TRUE)[[1]]
  out <- character()
  i <- 1L
  n <- length(chars)

  while (i <= n) {
    starts_xi <- i + 2L <= n && paste0(chars[i], chars[i + 1L], chars[i + 2L]) == "xi["
    boundary <- i == 1L || !grepl("[[:alnum:]_.]", chars[i - 1L])
    if (starts_xi && boundary) {
      end <- .pop_hyp_find_matching_bracket(chars, i + 2L)
      content <- paste(chars[(i + 3L):(end - 1L)], collapse = "")
      out <- c(out, paste0(".xi(", encodeString(content, quote = "\""), ")"))
      i <- end + 1L
    } else {
      out <- c(out, chars[i])
      i <- i + 1L
    }
  }
  paste(out, collapse = "")
}

.pop_hyp_find_matching_bracket <- function(chars, start) {
  depth <- 0L
  quote <- NULL
  escape <- FALSE
  for (i in seq.int(start, length(chars))) {
    ch <- chars[i]
    if (!is.null(quote)) {
      if (escape) {
        escape <- FALSE
      } else if (ch == "\\") {
        escape <- TRUE
      } else if (ch == quote) {
        quote <- NULL
      }
      next
    }
    if (ch %in% c("'", "\"")) {
      quote <- ch
      next
    }
    if (ch == "[") {
      depth <- depth + 1L
    } else if (ch == "]") {
      depth <- depth - 1L
      if (depth == 0L) return(i)
    }
  }
  stop("Unmatched '[' in xi[...] expression.", call. = FALSE)
}

.pop_hyp_select_groups <- function(filter, groups) {
  filter <- trimws(filter)
  parsed <- try(parse(text = filter), silent = TRUE)
  if (inherits(parsed, "try-error")) {
    stop("Could not parse xi[...] filter: ", filter, call. = FALSE)
  }

  if (length(parsed) == 1L && .pop_hyp_is_constant(parsed[[1L]])) {
    key <- eval(parsed[[1L]], envir = baseenv())
    mask <- .pop_hyp_match_key(groups, key)
  } else {
    env <- list2env(as.list(groups), parent = baseenv())
    mask <- eval(parsed, envir = env)
  }

  if (!is.logical(mask)) {
    stop("xi[...] filter must evaluate to a logical vector.", call. = FALSE)
  }
  if (length(mask) == 1L) {
    mask <- rep(mask, nrow(groups))
  }
  if (length(mask) != nrow(groups)) {
    stop("xi[...] filter must have length 1 or the number of groups.",
         call. = FALSE)
  }
  mask[is.na(mask)] <- FALSE
  idx <- which(mask)
  if (length(idx) == 0L) {
    stop("xi[...] filter matched no groups. Available groups: ",
         paste(groups$.group, collapse = ", "), call. = FALSE)
  }
  idx
}

.pop_hyp_is_constant <- function(expr) {
  is.numeric(expr) || is.character(expr)
}

.pop_hyp_match_key <- function(groups, key) {
  vars <- .pop_hyp_group_vars(groups)
  if (length(vars) == 1L) {
    .pop_hyp_equal(groups[[vars]], key)
  } else {
    .pop_hyp_equal(groups$.group, key)
  }
}

.pop_hyp_group_vars <- function(groups) {
  vars <- setdiff(names(groups), ".group")
  if ("group" %in% vars &&
      identical(as.character(groups$group), as.character(groups$.group))) {
    vars <- setdiff(vars, "group")
  }
  vars
}

.pop_hyp_equal <- function(x, key) {
  out <- as.character(x) == as.character(key)
  nx <- suppressWarnings(as.numeric(as.character(x)))
  nk <- suppressWarnings(as.numeric(as.character(key)))
  if (is.finite(nk)) {
    out <- out | (!is.na(nx) & nx == nk)
  }
  out
}

.pop_hyp_value <- function(value, weights) {
  structure(list(value = as.numeric(value), weights = as.numeric(weights)),
            class = "uncounted_popsize_hyp_value")
}

.pop_hyp_as_value <- function(x, n) {
  if (inherits(x, "uncounted_popsize_hyp_value")) return(x)
  if (is.numeric(x) && length(x) == 1L) {
    return(.pop_hyp_value(x, numeric(n)))
  }
  stop("Hypothesis expression must evaluate to a single numeric value.",
       call. = FALSE)
}

#' @export
Ops.uncounted_popsize_hyp_value <- function(e1, e2) {
  unary <- nargs() == 1L
  if (unary) {
    x <- .pop_hyp_as_value(e1, length(e1$weights))
    if (.Generic == "+") return(x)
    if (.Generic == "-") return(.pop_hyp_value(-x$value, -x$weights))
    stop("Unsupported unary operator in hypothesis expression.", call. = FALSE)
  }

  n <- if (inherits(e1, "uncounted_popsize_hyp_value")) {
    length(e1$weights)
  } else {
    length(e2$weights)
  }
  a <- .pop_hyp_as_value(e1, n)
  b <- .pop_hyp_as_value(e2, n)

  switch(.Generic,
    "+" = .pop_hyp_value(a$value + b$value, a$weights + b$weights),
    "-" = .pop_hyp_value(a$value - b$value, a$weights - b$weights),
    "*" = .pop_hyp_value(a$value * b$value,
                         a$weights * b$value + b$weights * a$value),
    "/" = .pop_hyp_value(a$value / b$value,
                         (a$weights * b$value - b$weights * a$value) /
                           (b$value^2)),
    "^" = {
      if (any(b$weights != 0)) {
        stop("Powers with an estimated exponent are not supported.",
             call. = FALSE)
      }
      .pop_hyp_value(a$value^b$value,
                     b$value * a$value^(b$value - 1) * a$weights)
    },
    stop("Unsupported operator in hypothesis expression: ", .Generic,
         call. = FALSE)
  )
}

#' @export
Math.uncounted_popsize_hyp_value <- function(x, ...) {
  switch(.Generic,
    log = .pop_hyp_value(log(x$value), x$weights / x$value),
    exp = .pop_hyp_value(exp(x$value), exp(x$value) * x$weights),
    sqrt = .pop_hyp_value(sqrt(x$value), x$weights / (2 * sqrt(x$value))),
    stop("Unsupported function in hypothesis expression: ", .Generic,
         call. = FALSE)
  )
}

.pop_hyp_call_fun <- function(fun, x, groups) {
  args <- names(formals(fun))
  out <- if (!is.null(args) && "groups" %in% args) {
    fun(x, groups = groups)
  } else {
    fun(x)
  }
  if (!is.numeric(out)) {
    stop("Function hypothesis must return a numeric vector.", call. = FALSE)
  }
  if (is.null(names(out))) {
    names(out) <- paste0("function_", seq_along(out))
  }
  stats::setNames(as.numeric(out), names(out))
}

.pop_hyp_jacobian <- function(fun, x, groups, f0) {
  p <- length(x)
  q <- length(f0)
  J <- matrix(NA_real_, nrow = q, ncol = p,
              dimnames = list(names(f0), names(x)))
  step <- sqrt(.Machine$double.eps) * pmax(abs(x), 1)
  for (j in seq_len(p)) {
    xp <- xm <- x
    xp[j] <- xp[j] + step[j]
    xm[j] <- xm[j] - step[j]
    fp <- .pop_hyp_call_fun(fun, xp, groups)
    fm <- .pop_hyp_call_fun(fun, xm, groups)
    if (length(fp) != q || length(fm) != q) {
      stop("Function hypothesis must return a stable-length numeric vector.",
           call. = FALSE)
    }
    J[, j] <- (fp - fm) / (2 * step[j])
  }
  J
}

.pop_hyp_expand_null <- function(null, x) {
  if (length(null) == 1L && is.null(names(null))) {
    return(rep(null, length(x)))
  }
  if (!is.null(names(null))) {
    if (!all(names(x) %in% names(null))) {
      stop("Named numeric hypotheses must include all xi group names.",
           call. = FALSE)
    }
    return(as.numeric(null[names(x)]))
  }
  if (length(null) != length(x)) {
    stop("Numeric hypotheses must have length 1 or match the number of xi estimates.",
         call. = FALSE)
  }
  as.numeric(null)
}

.pop_hyp_groups_from_popsize <- function(object) {
  groups <- attr(object, "groups")
  if (is.null(groups)) {
    groups <- data.frame(.group = object$group, group = object$group,
                         stringsAsFactors = FALSE)
  }
  groups$.group <- as.character(groups$.group)
  groups
}

.pop_hyp_groups_from_boot <- function(object, labels) {
  groups <- object$groups
  if (is.null(groups)) {
    groups <- attr(object$popsize_full, "groups")
  }
  if (is.null(groups) || nrow(groups) != length(labels)) {
    groups <- data.frame(.group = labels, group = labels,
                         stringsAsFactors = FALSE)
  }
  groups$.group <- as.character(groups$.group)
  groups
}

.pop_hyp_append_total <- function(x, V, groups) {
  total_x <- sum(x)
  cov_total <- colSums(V)
  var_total <- sum(V)
  V2 <- rbind(cbind(V, Total = cov_total),
              Total = c(cov_total, Total = var_total))
  x2 <- c(x, Total = total_x)
  list(x = x2, V = V2, groups = .pop_hyp_append_total_group(groups))
}

.pop_hyp_append_total_group <- function(groups) {
  total_row <- groups[NA_integer_, , drop = FALSE][1, , drop = FALSE]
  total_row$.group <- "Total"
  if ("group" %in% names(total_row)) {
    if (is.factor(groups$group)) {
      groups$group <- factor(as.character(groups$group),
                             levels = unique(c(levels(groups$group), "Total")))
      total_row$group <- factor("Total", levels = levels(groups$group))
    } else {
      total_row$group <- "Total"
    }
  }
  rbind(groups, total_row)
}

.pop_hyp_boot_total_value <- function(object, estimate, x, total_draws) {
  total_row <- object$popsize_full[object$popsize_full$group == "Total", ,
                                   drop = FALSE]
  if (nrow(total_row) == 1L && estimate %in% names(total_row) &&
      is.finite(total_row[[estimate]])) {
    return(total_row[[estimate]])
  }
  switch(estimate,
    plugin = sum(x),
    boot_median = stats::median(total_draws[is.finite(total_draws)]),
    boot_mean = mean(total_draws[is.finite(total_draws)])
  )
}
