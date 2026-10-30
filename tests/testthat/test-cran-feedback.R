# ---- Detection-feedback sensitivity (CRAN) ----

feedback_fit <- function(...) {
  estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005, ...)
}

quiet_feedback <- function(...) suppressWarnings(feedback_sensitivity(...))

grid <- seq(-0.2, 0.2, by = 0.05)

test_that("anchor is required and arguments are validated", {
  fit <- feedback_fit()
  expect_error(feedback_sensitivity(fit), "'anchor' must be supplied")
  expect_error(feedback_sensitivity(fit, anchor = "other"))
  expect_error(feedback_sensitivity(fit, "geomean", delta_grid = c(-1, 0)), "greater than -1")
  expect_error(feedback_sensitivity(fit, "geomean", delta_grid = c(-1.5, 0)), "greater than -1")
  expect_error(feedback_sensitivity(fit, "geomean", delta_grid = c(0.1, 0.2)), "contain 0")
  expect_error(feedback_sensitivity(fit, "value"), "requires 'K'")
  expect_error(feedback_sensitivity(fit, "value", K = 0), "positive")
  expect_error(feedback_sensitivity(fit, "value", K = -1), "positive")
  expect_error(feedback_sensitivity(fit, "value", K = c(1, 2)), "scalar")
  expect_error(feedback_sensitivity(fit, "geomean", K = 10), "only used")
  expect_error(feedback_sensitivity(testdata, "geomean"), "uncounted")
})

test_that("delta = 0 reproduces the popsize() plug-in estimates", {
  specs <- list(
    list(),
    list(cov_alpha = ~sex),
    list(method = "nb", cov_alpha = ~year),
    list(cov_alpha = ~sex, constrained = TRUE),
    list(link_rho = "cloglog", cov_alpha = ~sex)
  )
  for (s in specs) {
    fit <- do.call(feedback_fit, s)
    for (a in c("geomean", "share")) {
      x <- quiet_feedback(fit, a, delta_grid = grid)
      p0 <- x$profile[x$profile$delta == 0, ]
      ps <- popsize(fit)
      expect_equal(p0$estimate, ps$estimate, tolerance = 1e-10)
      expect_identical(p0$group, as.character(ps$group))
    }
  }
})

test_that("a grid value within rounding of 0 counts as 0", {
  x <- quiet_feedback(feedback_fit(), "geomean")   # seq() leaves about 5e-17
  expect_true(0 %in% x$delta_grid)
  expect_equal(nrow(x$profile), length(x$delta_grid))
})

test_that("the transform preserves the fit", {
  fit <- feedback_fit(cov_alpha = ~sex)
  rho <- fit$rho_values
  for (a in c("geomean", "share")) {
    x <- quiet_feedback(fit, a, delta_grid = grid, keep_cells = TRUE)
    K <- x$anchor$K
    for (j in seq_along(grid)) {
      xi <- x$cells$xi[, j]
      implied <- rho * xi * (xi / K)^x$delta_grid[j]
      expect_equal(unname(implied / fitted(fit)), rep(1, length(xi)), tolerance = 1e-8)
    }
  }
})

test_that("a cell at the anchor does not move", {
  fit <- feedback_fit()
  xi0 <- fit$N^fit$alpha_values
  x <- quiet_feedback(fit, "value", K = xi0[7], delta_grid = grid, keep_cells = TRUE)
  expect_equal(unname(x$cells$xi[7, ]), rep(xi0[7], length(grid)), tolerance = 1e-12)
})

test_that("with a constant anchor the slope on log N is alpha / (1 + delta)", {
  fit <- feedback_fit()
  x <- quiet_feedback(fit, "geomean", delta_grid = grid, keep_cells = TRUE)
  for (j in seq_along(grid)) {
    slope <- unname(coef(lm(log(x$cells$xi[, j]) ~ log(fit$N)))[2])
    expect_equal(slope, unname(coef(fit)["alpha"]) / (1 + grid[j]), tolerance = 1e-8)
    expect_equal(unname(x$cells$alpha[1, j]), unname(coef(fit)["alpha"]) / (1 + grid[j]))
  }
  s <- quiet_feedback(fit, "share", delta_grid = grid, keep_cells = TRUE)
  expect_equal(unname(s$cells$alpha[1, 1]), 1 - (1 - unname(coef(fit)["alpha"])) / 0.8)
  v <- quiet_feedback(fit, "value", K = 10, delta_grid = grid, keep_cells = TRUE)
  expect_true(all(is.na(v$cells$alpha)))
})

test_that("the geomean total is non-increasing and the derivative is right", {
  fit <- feedback_fit(cov_alpha = ~sex)
  x <- quiet_feedback(fit, "geomean", by = ~sex, total = TRUE,
                      delta_grid = seq(-0.2, 0.3, by = 0.01))
  tot <- x$profile[x$profile$group == "Total", ]
  expect_true(all(diff(tot$estimate) <= 0))

  ## The derivative matches finite differences of the total
  h <- 1e-6
  g <- function(d) sum(exp(x$log_K + (x$log_xi0 - x$log_K) / (1 + d)))
  for (d in c(-0.1, 0, 0.1)) {
    j <- which.min(abs(x$diagnostics$delta - d))
    expect_equal(x$diagnostics$d_total[j],
                 (g(x$diagnostics$delta[j] + h) - g(x$diagnostics$delta[j] - h)) / (2 * h),
                 tolerance = 1e-5)
  }
})

test_that("the anchor can reverse the direction of the total", {
  fit <- estimate_hidden_pop(irregular_migration, ~m, ~n, ~N, gamma = 0.005,
                             countries = ~country_code)
  tot <- function(a) {
    p <- quiet_feedback(fit, a, delta_grid = c(0, 0.2))$profile
    p$estimate[p$delta == 0.2] - p$estimate[p$delta == 0]
  }
  expect_lt(tot("geomean"), 0)
  expect_gt(tot("share"), 0)
})

test_that("A1 binds at max(alpha) - 1 to grid precision", {
  ## Fine grid, and no other constraint binds first
  fit <- feedback_fit(cov_alpha = ~sex)
  d_grid <- seq(-0.5, 0.2, by = 0.001)
  x <- quiet_feedback(fit, "geomean", delta_grid = d_grid, p_max = Inf)
  limit <- max(fit$alpha_values) - 1
  expect_lt(abs(x$admissible_range[1] - limit), 0.001 + 1e-9)
  expect_gt(x$admissible_range[1], limit)
  low <- x$diagnostics[x$diagnostics$delta < limit, ]
  expect_true(all(grepl("A1", low$binding)))

  ## Not applied to the share anchor, to a varying K, or when switched off
  expect_true(all(is.na(quiet_feedback(fit, "share", delta_grid = grid)$diagnostics$A1)))
  expect_true(all(is.na(quiet_feedback(fit, "value", K = fit$N, delta_grid = grid)$diagnostics$A1)))
  expect_false(anyNA(quiet_feedback(fit, "value", K = 50, delta_grid = grid)$diagnostics$A1))
  expect_true(all(is.na(quiet_feedback(fit, "geomean", delta_grid = grid,
                                       sublinear = FALSE)$diagnostics$A1)))
})

test_that("A2 uses the implied detection and leaves out baseline violators", {
  fit <- feedback_fit()
  x0 <- fit$N^fit$alpha_values
  p0 <- fitted(fit) / x0

  ## A cap below the largest baseline detection excludes the cells above it
  cap <- sort(p0, decreasing = TRUE)[3]
  expect_warning(x <- feedback_sensitivity(fit, "geomean", delta_grid = grid, p_max = cap),
                 "left out of constraint A2")
  expect_equal(sort(x$a2_excluded$row), sort(which(p0 > cap)))
  expect_true(x$diagnostics$A2[x$diagnostics$delta == 0])
  expect_false(anyNA(x$admissible_range))

  ## The rest of the cells still bind
  x <- quiet_feedback(fit, "value", K = 1e6, delta_grid = c(-0.5, 0, 0.5),
                      p_max = max(p0), sublinear = FALSE)
  expect_false(x$diagnostics$A2[x$diagnostics$delta == -0.5])
  expect_true(x$diagnostics$A2[x$diagnostics$delta == 0.5])
  expect_equal(x$diagnostics$max_p[x$diagnostics$delta == 0.5],
               max(fitted(fit) / (1e6 * (x0 / 1e6)^(1 / 1.5))))
})

test_that("on the bundled data delta = 0 stays admissible under the power link", {
  ## PNG 2023 Male (N = 1, n = 1) has fitted detection above 1 at delta = 0
  d <- irregular_migration
  d$year <- factor(d$year)
  d$ukr <- as.integer(d$country_code == "UKR")
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, cov_alpha = ~year * ukr + sex, cov_beta = ~year,
                             countries = ~country_code)
  expect_warning(x <- feedback_sensitivity(fit, "geomean"), "PNG")
  expect_equal(nrow(x$a2_excluded), 1L)
  expect_equal(x$admissible_range[1], -0.09, tolerance = 1e-8)
  expect_equal(x$admissible_range[2], 0.30, tolerance = 1e-8)
  expect_equal(x$diagnostics$binding[which.min(abs(x$diagnostics$delta + 0.10))], "A2")
})

test_that("A3 compares group totals with the observed counts", {
  fit <- feedback_fit(cov_alpha = ~sex)
  x <- quiet_feedback(fit, "value", K = 1e-3, by = ~sex, delta_grid = c(0, 5))
  m_group <- rowsum(testdata$m, testdata$sex)
  p5 <- x$profile[x$profile$delta == 5, ]
  expect_false(all(p5$estimate >= m_group))
  expect_false(x$diagnostics$A3[x$diagnostics$delta == 5])
  expect_true(x$diagnostics$A3[x$diagnostics$delta == 0])
})

test_that("zero being inadmissible, and gaps, give warnings", {
  fit <- feedback_fit()
  ## Force alpha >= 1 so that A1 fails at 0
  fit$alpha_values[1] <- 1.2
  expect_warning(x <- feedback_sensitivity(fit, "geomean", delta_grid = grid, p_max = Inf),
                 "inadmissible")
  expect_true(all(is.na(x$admissible_range)))

  ## Contiguity is checked on the run around 0
  expect_warning(r <- .feedback_admissible_range(c(-2, -1, 0, 1, 2),
                                                 c(TRUE, FALSE, TRUE, TRUE, FALSE)),
                 "not contiguous")
  expect_equal(r, c(0, 1))
})

test_that("anchor_by computes the anchor within groups", {
  fit <- feedback_fit(cov_alpha = ~sex)
  x <- quiet_feedback(fit, "geomean", anchor_by = ~sex, delta_grid = grid)
  lx <- x$log_xi0
  for (s in unique(testdata$sex)) {
    i <- testdata$sex == s
    expect_equal(unique(round(x$log_K[i], 8)), round(mean(lx[i]), 8))
  }
  s <- quiet_feedback(fit, "share", anchor_by = ~sex, delta_grid = grid)
  i <- testdata$sex == "M"
  expect_equal(s$log_K[i] - log(testdata$N[i]),
               rep(mean(lx[i] - log(testdata$N[i])), sum(i)))
  expect_error(feedback_sensitivity(fit, "value", K = 5, anchor_by = ~sex), "anchor_by")
})

test_that("bootstrap replicates are transformed with the anchor held fixed", {
  skip_if_not_installed("fwb")
  fit <- feedback_fit(cov_alpha = ~sex)
  b <- bootstrap_popsize(fit, R = 20, seed = 1, verbose = FALSE, by = ~sex, total = TRUE)
  x <- quiet_feedback(fit, "geomean", by = ~sex, total = TRUE, boot = b, delta_grid = grid)

  ## At delta = 0 the replicates are the bootstrap's own totals
  p0 <- x$profile[x$profile$delta == 0, ]
  expect_equal(p0$boot_median[p0$group != "Total"], b$popsize_full$boot_median[b$popsize_full$group != "Total"])
  expect_equal(p0$boot_median[p0$group == "Total"], b$total$median)
  expect_equal(p0$lower[p0$group == "Total"], b$total$lower)
  expect_equal(p0$upper[p0$group == "Total"], b$total$upper)

  ## The same K in every replicate
  d <- 0.1
  a1 <- b$boot_params[1, grep("alpha", colnames(b$boot_params))]
  al <- as.numeric(fit$X_alpha %*% a1)
  by_hand <- sum(exp(x$log_K + (al * log(fit$N) - x$log_K) / 1.1))
  expect_equal(colSums(exp(x$log_K + (x$replicates - x$log_K) / 1.1))[1], by_hand)
  expect_true(all(is.finite(x$profile$lower)))
})

test_that("a bootstrap of another fit is refused", {
  skip_if_not_installed("fwb")
  fit <- feedback_fit()
  other <- feedback_fit(cov_alpha = ~sex)
  b <- bootstrap_popsize(other, R = 3, seed = 1, verbose = FALSE, ci_type = "cheap")
  expect_error(feedback_sensitivity(fit, "geomean", boot = b), "alpha coefficients")
  fit2 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.01)
  b2 <- bootstrap_popsize(fit2, R = 3, seed = 1, verbose = FALSE, ci_type = "cheap")
  expect_error(feedback_sensitivity(fit, "geomean", boot = b2), "same fit")
  expect_error(feedback_sensitivity(fit, "geomean", boot = list()), "bootstrap_popsize")
})

test_that("feedback works across methods, links and constraints with a bootstrap", {
  skip_if_not_installed("fwb")
  specs <- list(
    list(method = "nb"),
    list(constrained = TRUE, cov_alpha = ~sex),
    list(link_rho = "cloglog"),
    list(cov_alpha = ~sex * year)
  )
  for (s in specs) {
    fit <- do.call(estimate_hidden_pop,
                   c(list(testdata, ~m, ~n, ~N, gamma = 0.005, countries = ~country), s))
    b <- suppressWarnings(bootstrap_popsize(fit, R = 12, cluster = ~country, seed = 1,
                                            verbose = FALSE, by = ~year, total = TRUE))
    x <- quiet_feedback(fit, "geomean", by = ~year, total = TRUE, boot = b, delta_grid = grid)
    p0 <- x$profile[x$profile$delta == 0 & x$profile$group == "Total", ]
    expect_equal(p0$estimate, sum(popsize(fit)$estimate))
    expect_equal(p0$boot_median, b$total$median)
  }
})

test_that("the Monte Carlo version uses one delta per replicate", {
  skip_if_not_installed("fwb")
  fit <- feedback_fit()
  b <- bootstrap_popsize(fit, R = 20, seed = 1, verbose = FALSE)
  x <- quiet_feedback(fit, "geomean", boot = b, delta_grid = grid,
                      delta_prior = function(R) rep(0, R))
  ## A point mass at 0 reproduces the bootstrap
  expect_equal(x$monte_carlo$median, b$popsize_full$boot_median)
  expect_equal(x$monte_carlo$lower, b$popsize_full$lower)

  set.seed(1)
  u <- quiet_feedback(fit, "geomean", boot = b, delta_grid = grid, delta_prior = TRUE)
  d <- attr(u$monte_carlo, "delta")
  expect_length(d, 20L)
  expect_true(all(d >= u$admissible_range[1] & d <= u$admissible_range[2]))
  expect_error(quiet_feedback(fit, "geomean", delta_prior = TRUE), "requires 'boot'")
  expect_error(quiet_feedback(fit, "geomean", boot = b, delta_prior = function(R) rep(-2, R)),
               "above -1")
  expect_error(quiet_feedback(fit, "geomean", boot = b, delta_prior = "x"), "function")
})

test_that("print and plot run", {
  skip_if_not_installed("fwb")
  fit <- feedback_fit(cov_alpha = ~sex)
  x <- quiet_feedback(fit, "geomean", by = ~sex, total = TRUE)
  out <- capture.output(print(x))
  expect_true(any(grepl("Anchor: geomean", out)))
  expect_true(any(grepl("Admissible range", out)))
  expect_true(any(grepl("d=0.1", out)))

  b <- bootstrap_popsize(fit, R = 12, seed = 1, verbose = FALSE, by = ~sex, total = TRUE)
  xb <- quiet_feedback(fit, "share", by = ~sex, total = TRUE, boot = b)
  expect_true(any(grepl("\\[", capture.output(print(xb)))))

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(x))
  expect_no_error(plot(xb, group = c("F", "M"), plausible = c(-0.1, 0.1)))
  expect_error(plot(x, group = "nope"), "Unknown group")
})

# ---- feedback_identified_set() ----

test_that("the identified set finds an interior minimum", {
  ## Value anchor, one cell at 10K and sixty at K/10: over [-0.2, 0.2] the minimum is
  ## interior, at 10^(1 / (1 + d)) = sqrt(60), with total 2 sqrt(60) K
  K <- 3
  x <- quiet_feedback(feedback_fit(), "value", K = K, delta_grid = grid, p_max = Inf)
  xi0 <- c(10 * K, rep(K / 10, 60))
  x$log_xi0 <- log(xi0)
  x$log_K <- rep(log(K), 61)
  x$group_id <- rep(1L, 61)
  x$labels <- "(all)"
  set <- feedback_identified_set(x, delta_range = c(-0.2, 0.2))

  expect_equal(set$delta_at_min, 1 / log10(sqrt(60)) - 1, tolerance = 1e-5)
  expect_equal(set$L, 2 * sqrt(60) * K, tolerance = 1e-8)
  end_value <- K * (10^(1 / 1.2) + 60 * 10^(-1 / 1.2))
  expect_lt(set$L, end_value)
  expect_equal(end_value / K, 15.62, tolerance = 1e-3)
})

test_that("the bounds match a brute-force search", {
  fit <- feedback_fit(cov_alpha = ~sex)
  x <- quiet_feedback(fit, "value", K = 40, by = ~sex, total = TRUE, delta_grid = grid)
  D <- c(-0.15, 0.25)
  set <- feedback_identified_set(x, delta_range = D)
  dense <- seq(D[1], D[2], length.out = 4001)
  tot <- sapply(dense, function(d) sum(exp(x$log_K + (x$log_xi0 - x$log_K) / (1 + d))))
  row <- set[set$group == "Total", ]
  expect_equal(row$L, min(tot), tolerance = 1e-7)
  expect_equal(row$U, max(tot), tolerance = 1e-7)
  expect_equal(nrow(set), 3L)
  expect_true(all(is.na(set$ci_lower)))
  expect_output(print(set), "Identified set")
})

test_that("the default range is the admissible range", {
  x <- quiet_feedback(feedback_fit(), "geomean")
  set <- feedback_identified_set(x)
  expect_equal(attr(set, "delta_range"), x$admissible_range)
  expect_error(feedback_identified_set(x, delta_range = c(0.2, 0.1)), "delta range")
  expect_error(feedback_identified_set(list()), "feedback_sensitivity")
})

test_that("the Imbens-Manski constant stays between its limits", {
  expect_equal(.im_critical(0, 0.5, 0.95), qnorm(0.975), tolerance = 1e-6)
  expect_equal(.im_critical(100, 0.5, 0.95), qnorm(0.95), tolerance = 1e-6)
  for (r in c(0.01, 0.1, 0.5, 1, 2, 3)) {
    cc <- .im_critical(r, 1, 0.95)
    expect_gt(cc, qnorm(0.95))
    expect_lt(cc, qnorm(0.975))
    expect_equal(pnorm(cc + r) - pnorm(-cc), 0.95, tolerance = 1e-6)
  }
  cs <- sapply(c(0.1, 0.5, 1, 2), function(r) .im_critical(r, 1, 0.95))
  expect_true(all(diff(cs) < 0))
  expect_equal(.im_critical(1, 0, 0.95), qnorm(0.95))
  expect_equal(.im_critical(0, 0, 0.95), qnorm(0.975))
})

test_that("identified-set intervals follow their definitions", {
  skip_if_not_installed("fwb")
  fit <- feedback_fit()
  b <- bootstrap_popsize(fit, R = 30, seed = 1, verbose = FALSE)
  x <- quiet_feedback(fit, "geomean", boot = b)
  D <- x$admissible_range
  im <- feedback_identified_set(x)
  u <- feedback_identified_set(x, method = "union")

  bounds <- apply(x$replicates, 2, function(a) {
    tot <- sapply(seq(D[1], D[2], length.out = 41),
                  function(d) sum(exp(x$log_K + (a - x$log_K) / (1 + d))))
    c(min(tot), max(tot))
  })
  expect_equal(im$se_L, sd(log(bounds[1, ])), tolerance = 1e-4)
  expect_equal(im$se_U, sd(log(bounds[2, ])), tolerance = 1e-4)
  expect_equal(im$c, .im_critical(im$Delta, max(im$se_L, im$se_U), 0.95))
  expect_equal(im$ci_lower, exp(log(im$L) - im$c * im$se_L))
  expect_equal(im$ci_upper, exp(log(im$U) + im$c * im$se_U))
  expect_equal(u$ci_lower, unname(quantile(bounds[1, ], 0.025)), tolerance = 1e-4)
  expect_equal(u$ci_upper, unname(quantile(bounds[2, ], 0.975)), tolerance = 1e-4)

  ## A one-point range is point identification: the two-sided constant
  pt <- feedback_identified_set(x, delta_range = c(0, 0))
  expect_equal(pt$Delta, 0)
  expect_equal(pt$c, qnorm(0.975), tolerance = 1e-6)

  iq <- feedback_identified_set(x, se = "iqr")
  expect_equal(iq$se_L, IQR(log(bounds[1, ])) / 1.349, tolerance = 1e-4)
})

# ---- feedback_tipping() ----

tipping_fit <- function() {
  feedback_fit(cov_alpha = ~year, countries = ~country)
}

test_that("a threshold crossing is found exactly", {
  fit <- feedback_fit(cov_alpha = ~sex)
  x <- quiet_feedback(fit, "geomean", by = ~sex, total = TRUE)
  tot <- function(d) sum(exp(x$log_K + (x$log_xi0 - x$log_K) / (1 + d)))
  ## The total is decreasing, so a higher threshold is crossed below 0
  thr <- tot(-0.07)
  tp <- feedback_tipping(x, threshold = thr)
  expect_equal(unname(tp$point["negative"]), -0.07, tolerance = 1e-5)
  expect_true(is.na(tp$point["positive"]))
  expect_equal(tp$alternative, "less")
  expect_equal(tp$baseline, tot(0) - thr)
  expect_match(tp$statement, "Total")

  low <- feedback_tipping(x, threshold = tot(0.11))
  expect_equal(unname(low$point["positive"]), 0.11, tolerance = 1e-5)
  expect_equal(low$alternative, "greater")
  expect_true(is.na(low$supported_at_zero))

  none <- feedback_tipping(x, threshold = 1)
  expect_true(all(is.na(none$point)))

  one <- feedback_tipping(x, threshold = tot(0) / 2, group = "M")
  expect_match(one$statement, "\\[M\\]")
})

test_that("a contrast uses the hypotheses_popsize syntax", {
  fit <- tipping_fit()
  x <- quiet_feedback(fit, "geomean", by = ~year)
  ctr <- function(d) {
    z <- tapply(exp(x$log_K + (x$log_xi0 - x$log_K) / (1 + d)), testdata$year, sum)
    unname(z["2023"] - z["2019"])
  }
  h <- "xi[year == 2023] - xi[year == 2019] > 0"
  tp <- feedback_tipping(x, hypothesis = h)
  expect_equal(tp$baseline, ctr(0))
  expect_equal(tp$alternative, "greater")
  expect_equal(tp$curve$contrast[tp$curve$delta == 0], ctr(0))
  for (r in tp$point[is.finite(tp$point)]) expect_equal(ctr(r), 0, tolerance = 1e-3)
  expect_equal(feedback_tipping(x, hypothesis = "xi[year == 2023] < 100")$alternative, "less")
})

test_that("the test stops supporting a statement at the first failing grid value", {
  skip_if_not_installed("fwb")
  fit <- feedback_fit()
  b <- bootstrap_popsize(fit, R = 30, seed = 1, verbose = FALSE)
  x <- quiet_feedback(fit, "geomean", boot = b)
  low <- min(x$profile$lower[x$profile$delta == 0]) / 2
  tp <- feedback_tipping(x, threshold = low)
  expect_true(tp$supported_at_zero)
  expect_equal(tp$curve$lower[tp$curve$delta == 0],
               unname(quantile(colSums(exp(x$replicates)), 0.05)) - low)

  ## The first grid value above 0 where the 5% quantile is not above the threshold
  above <- tp$curve[tp$curve$delta > 0, ]
  first <- above$delta[which(!above$supported)[1]]
  expect_equal(unname(tp$test["positive"]), first)
  expect_output(print(tp), "stops supporting")

  ## Not supported at 0: nothing to overturn
  mid <- feedback_tipping(x, threshold = x$profile$estimate[x$profile$delta == 0] * 0.999)
  expect_false(mid$supported_at_zero)
  expect_true(all(is.na(mid$test)))
  expect_output(print(mid), "no significant result to overturn")
})

test_that("feedback_tipping() validates its input", {
  x <- quiet_feedback(feedback_fit(cov_alpha = ~sex), "geomean", by = ~sex)
  expect_error(feedback_tipping(x), "either")
  expect_error(feedback_tipping(x, threshold = 1, hypothesis = "xi[TRUE] > 0"), "either")
  expect_error(feedback_tipping(x, threshold = 1), "'group'")
  expect_error(feedback_tipping(x, threshold = 1, group = "nope"), "'group'")
  expect_error(feedback_tipping(x, threshold = c(1, 2), group = "F"), "single")
  expect_error(feedback_tipping(list(), threshold = 1), "feedback_sensitivity")
})

# ---- feedback_calibrate() ----

test_that("feedback_calibrate() returns the share and its mean", {
  fit <- feedback_fit(cov_alpha = ~sex)
  cal <- feedback_calibrate(fit)
  xi0 <- fit$N^fit$alpha_values
  s <- xi0 / (fit$N + xi0)
  expect_equal(cal$s, s)
  expect_equal(cal$mean_s, mean(s))
  expect_equal(cal$weighted_mean_s, sum(s * xi0) / sum(xi0))
  expect_equal(cal$congestion_elasticity, -cal$weighted_mean_s)
  expect_output(print(cal), "congestion elasticity")
  expect_error(feedback_calibrate(list()), "uncounted")
})

test_that("a subset of the identified set still prints", {
  x <- quiet_feedback(feedback_fit(cov_alpha = ~sex), "geomean", by = ~sex, total = TRUE)
  set <- feedback_identified_set(x)
  expect_output(print(set[3, c("group", "L")]), "Total")
})
