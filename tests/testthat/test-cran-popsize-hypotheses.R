# ---- Population-size hypothesis tests (CRAN) ----

.hyp_year_popsize <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      countries <- unique(testdata$country)[1:8]
      d <- testdata[testdata$country %in% countries &
                      testdata$year %in% c("2019", "2023"), ]
      d <- droplevels(d)
      fit <- quick_fit(d, gamma = 0.005)
      value <<- popsize(fit, by = ~year)
    }
    value
  }
})

.hyp_year_sex_popsize <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      countries <- unique(testdata$country)[1:8]
      d <- testdata[testdata$country %in% countries &
                      testdata$year %in% c("2019", "2023"), ]
      d <- droplevels(d)
      fit <- quick_fit(d, gamma = 0.005)
      value <<- popsize(fit, by = ~year + sex)
    }
    value
  }
})

.hyp_boot_object <- function() {
  draws <- cbind(`2019` = 1:12, `2023` = 2 * (1:12))
  groups <- data.frame(
    .group = c("2019", "2023"),
    year = factor(c("2019", "2023")),
    group = c("2019", "2023"),
    stringsAsFactors = FALSE
  )
  popsize_full <- data.frame(
    group = c("2019", "2023"),
    plugin = c(6, 14),
    plugin_bc = NA_real_,
    boot_median = c(stats::median(draws[, 1]), stats::median(draws[, 2])),
    boot_mean = c(mean(draws[, 1]), mean(draws[, 2])),
    lower = NA_real_,
    upper = NA_real_,
    stringsAsFactors = FALSE
  )
  out <- list(
    t = draws,
    t0 = popsize_full$plugin,
    t0_bc = rep(NA_real_, 2),
    popsize = popsize_full[, c("group", "boot_median", "lower", "upper")],
    popsize_full = popsize_full,
    total = NULL,
    R = nrow(draws),
    ci_type = "perc",
    point_estimate = "median",
    level = 0.95,
    boot_params = NULL,
    n_converged = nrow(draws),
    cluster = FALSE,
    groups = groups
  )
  names(out$popsize)[2] <- "estimate"
  class(out) <- "uncounted_boot"
  out
}

test_that("hypotheses_popsize requires an explicit estimate", {
  ps <- .hyp_year_popsize()

  expect_error(hypotheses_popsize(ps, 100), "estimate")
})

test_that("numeric null tests every xi estimate", {
  ps <- .hyp_year_popsize()

  hyp <- hypotheses_popsize(ps, 100, estimate = "estimate")

  expect_s3_class(hyp, "uncounted_popsize_hypotheses")
  expect_equal(nrow(hyp), nrow(ps))
  expect_true(all(c("hypothesis", "null_hypothesis",
                    "alternative_hypothesis", "estimate", "std.error",
                    "p.value") %in% names(hyp)))
  expect_true(all(hyp$alternative == "two.sided"))
  expect_true(all(is.finite(hyp$std.error)))
})

test_that("xi filter syntax supports one-sided threshold tests", {
  ps <- .hyp_year_popsize()

  hyp <- hypotheses_popsize(
    ps,
    "xi[year == 2023] > 0",
    estimate = "estimate"
  )

  expect_equal(hyp$alternative, "greater")
  expect_equal(hyp$null_hypothesis, "xi[year == 2023] <= 0")
  expect_equal(hyp$alternative_hypothesis, "xi[year == 2023] > 0")
  expect_equal(hyp$estimate, ps$estimate[ps$group == "2023"])
  expect_true(is.finite(hyp$statistic))
  expect_true(hyp$p.value >= 0 && hyp$p.value <= 1)
})

test_that("xi filter syntax supports year-to-year increase tests", {
  ps <- .hyp_year_popsize()

  hyp_filter <- hypotheses_popsize(
    ps,
    "xi[year == 2023] - xi[year == 2019] > 0",
    estimate = "estimate"
  )
  hyp_pos <- hypotheses_popsize(ps, "b2 - b1 > 0", estimate = "estimate")

  expect_equal(hyp_filter$alternative, "greater")
  expect_equal(hyp_filter$contrast,
               ps$estimate[ps$group == "2023"] - ps$estimate[ps$group == "2019"])
  expect_equal(hyp_filter$contrast, hyp_pos$contrast)
  expect_equal(hyp_filter$std.error, hyp_pos$std.error)
})

test_that("equality tests are two-sided", {
  ps <- .hyp_year_popsize()

  hyp <- hypotheses_popsize(
    ps,
    "xi[year == 2023] = xi[year == 2019]",
    estimate = "estimate"
  )

  expect_equal(hyp$alternative, "two.sided")
  expect_equal(hyp$null_hypothesis, "xi[year == 2023] = xi[year == 2019]")
  expect_equal(hyp$alternative_hypothesis, "xi[year == 2023] != xi[year == 2019]")
  expect_true(hyp$p.value >= 0 && hyp$p.value <= 1)
})

test_that("hypothesis_side='null' reverses one-sided directional tests", {
  ps <- .hyp_year_popsize()

  as_alt <- hypotheses_popsize(ps, "xi[year == 2023] < 10000",
                               estimate = "estimate")
  as_null <- hypotheses_popsize(ps, "xi[year == 2023] < 10000",
                                estimate = "estimate",
                                hypothesis_side = "null")

  expect_equal(as_alt$alternative, "less")
  expect_equal(as_null$alternative, "greater")
  expect_equal(as_alt$null_hypothesis, "xi[year == 2023] >= 10000")
  expect_equal(as_null$null_hypothesis, "xi[year == 2023] <= 10000")
  expect_equal(as_alt$alternative_hypothesis, "xi[year == 2023] < 10000")
  expect_equal(as_null$alternative_hypothesis, "xi[year == 2023] > 10000")
  expect_equal(as_alt$p.value + as_null$p.value, 1, tolerance = 1e-12)
})

test_that("compact and label-based xi filters match explicit filters", {
  ps <- .hyp_year_popsize()

  explicit <- hypotheses_popsize(ps, "xi[year == 2023] > 0",
                                 estimate = "estimate")
  compact <- hypotheses_popsize(ps, "xi[2023] > 0", estimate = "estimate")
  label <- hypotheses_popsize(ps, "xi[group == '2023'] > 0",
                              estimate = "estimate")

  expect_equal(compact$estimate, explicit$estimate)
  expect_equal(label$estimate, explicit$estimate)
})

test_that("multi-row xi filters sum estimates and use the full covariance", {
  ps <- .hyp_year_sex_popsize()
  groups <- attr(ps, "groups")

  hyp <- hypotheses_popsize(
    ps,
    "xi[year == 2023] - xi[year == 2019] > 0",
    estimate = "estimate"
  )

  expected <- sum(ps$estimate[as.character(groups$year) == "2023"]) -
    sum(ps$estimate[as.character(groups$year) == "2019"])
  expect_equal(hyp$contrast, expected)
  expect_true(is.finite(hyp$std.error))
})

test_that("function hypotheses work with group metadata", {
  ps <- .hyp_year_popsize()

  hyp <- hypotheses_popsize(
    ps,
    function(x, groups) c(increase = unname(x["2023"] - x["2019"])),
    estimate = "estimate"
  )

  expect_equal(hyp$hypothesis, "increase")
  expect_equal(hyp$contrast,
               ps$estimate[ps$group == "2023"] - ps$estimate[ps$group == "2019"])
  expect_equal(hyp$alternative, "two.sided")
})

test_that("bootstrap hypotheses evaluate contrasts over draws", {
  boot <- .hyp_boot_object()

  hyp <- hypotheses_popsize(
    boot,
    "xi[year == 2023] - xi[year == 2019] > 0",
    estimate = "boot_mean"
  )

  expect_equal(hyp$method, "FWB Wald test")
  expect_equal(hyp$n_draws, 12)
  expect_equal(hyp$contrast, mean(boot$t[, 2]) - mean(boot$t[, 1]))
  expect_true(is.finite(hyp$std.error))
})

test_that("bootstrap hypotheses support compact xi filters and functions", {
  boot <- .hyp_boot_object()

  compact <- hypotheses_popsize(boot, "xi[2023] > 0",
                                estimate = "boot_median")
  fun <- hypotheses_popsize(
    boot,
    function(x, groups) c(diff = unname(x["2023"] - x["2019"])),
    estimate = "boot_mean"
  )

  expect_equal(compact$estimate, stats::median(boot$t[, 2]))
  expect_equal(fun$contrast, mean(boot$t[, 2]) - mean(boot$t[, 1]))
  expect_equal(fun$n_draws, 12)
})
