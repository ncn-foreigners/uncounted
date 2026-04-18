# ---- Frailty sensitivity tests (CRAN) ----

.frailty_data <- function() {
  d <- positive_data()
  keep_country <- unique(d$country)[1:12]
  d <- droplevels(d[d$country %in% keep_country, ])
  d$ukr <- as.integer(d$country %in% keep_country[1:3])
  d
}

test_that("frailty_sensitivity returns the expected structure", {
  d <- .frailty_data()
  fit <- quick_fit(
    d,
    method = "poisson",
    gamma = 0.005,
    cov_alpha = ~ year * ukr + sex,
    cov_beta = ~ year
  )

  fs <- frailty_sensitivity(
    fit,
    by = ~ year,
    r2_d = c(0, 0.05, 0.10),
    r2_y = c(0, 0.05, 0.10),
    plot = FALSE
  )

  expect_s3_class(fs, "uncounted_frailty_sensitivity")
  expect_true(all(c("baseline", "working", "surface", "robustness", "settings") %in%
                    names(fs)))
  expect_true(all(c("group", "theta_hat", "se_theta", "d_xi_dtheta", "df") %in%
                    names(fs$working)))
  expect_true(all(c("group", "r2_d", "r2_y", "bias_theta",
                    "xi_lower", "xi_upper",
                    "xi_ratio_lower", "xi_ratio_upper") %in%
                    names(fs$surface)))
  expect_true(all(c("group", "q", "rv_equal", "rv_extreme", "status") %in%
                    names(fs$robustness)))
})

test_that("frailty_sensitivity reproduces baseline Xi when r2_d = 0 or r2_y = 0", {
  d <- .frailty_data()
  fit <- quick_fit(
    d,
    method = "poisson",
    gamma = 0.005,
    cov_alpha = ~ year * ukr + sex,
    cov_beta = ~ year
  )

  fs <- frailty_sensitivity(
    fit,
    by = ~ year,
    r2_d = c(0, 0.10),
    r2_y = c(0, 0.10),
    plot = FALSE
  )

  zero_rows <- fs$surface$r2_d == 0 | fs$surface$r2_y == 0
  expect_equal(fs$surface$xi_ratio_lower[zero_rows],
               rep(1, sum(zero_rows)), tolerance = 1e-10)
  expect_equal(fs$surface$xi_ratio_upper[zero_rows],
               rep(1, sum(zero_rows)), tolerance = 1e-10)
})

test_that("frailty_sensitivity surface widens monotonically", {
  d <- .frailty_data()
  fit <- quick_fit(
    d,
    method = "poisson",
    gamma = 0.005,
    cov_alpha = ~ year * ukr + sex,
    cov_beta = ~ year
  )

  fs <- frailty_sensitivity(
    fit,
    by = ~ year,
    r2_d = c(0, 0.05, 0.10),
    r2_y = c(0.10),
    plot = FALSE
  )

  g0 <- fs$working$group[1]
  surf <- fs$surface[fs$surface$group == g0, ]

  expect_true(all(diff(surf$xi_lower) <= 1e-12))
  expect_true(all(diff(surf$xi_upper) >= -1e-12))
})

test_that("frailty_sensitivity matches popsize group labels and order", {
  d <- .frailty_data()
  fit <- quick_fit(
    d,
    method = "poisson",
    gamma = 0.005,
    cov_alpha = ~ year * ukr + sex,
    cov_beta = ~ year
  )
  ps <- popsize(fit, by = ~ year, bias_correction = FALSE, total = TRUE)
  fs <- frailty_sensitivity(fit, by = ~ year, plot = FALSE)

  expect_equal(
    fs$working$group,
    c(ps$group, "Total")
  )
})

test_that("frailty_sensitivity works for Poisson and NB vector-alpha fits", {
  d <- .frailty_data()

  for (method in c("poisson", "nb")) {
    fit <- quick_fit(
      d,
      method = method,
      gamma = 0.005,
      cov_alpha = ~ year * ukr + sex,
      cov_beta = ~ year
    )
    fs <- frailty_sensitivity(
      fit,
      by = ~ year,
      r2_d = c(0, 0.05),
      r2_y = c(0, 0.05),
      plot = FALSE
    )

    expect_true(all(is.finite(fs$working$xi_hat)), info = method)
    expect_true(all(is.finite(fs$working$se_theta)), info = method)
    expect_true(all(is.finite(fs$surface$xi_lower)), info = method)
    expect_true(all(is.finite(fs$surface$xi_upper)), info = method)
  }
})

test_that("frailty_sensitivity robustness table contains the requested targets", {
  d <- .frailty_data()
  fit <- quick_fit(
    d,
    method = "poisson",
    gamma = 0.005,
    cov_alpha = ~ year * ukr + sex,
    cov_beta = ~ year
  )
  fs <- frailty_sensitivity(
    fit,
    by = ~ year,
    q = c(0.10, 0.25),
    plot = FALSE
  )

  expect_true(all(c("rv_equal", "rv_extreme", "status") %in%
                    names(fs$robustness)))
  expect_equal(sort(unique(stats::na.omit(fs$robustness$q))), c(0.10, 0.25))
})

test_that("plot.uncounted_frailty_sensitivity runs without error", {
  d <- .frailty_data()
  fit <- quick_fit(
    d,
    method = "poisson",
    gamma = 0.005,
    cov_alpha = ~ year * ukr + sex,
    cov_beta = ~ year
  )
  fs <- frailty_sensitivity(
    fit,
    by = ~ year,
    r2_d = c(0, 0.05),
    r2_y = c(0, 0.05),
    plot = FALSE
  )

  expect_no_error(plot(fs))
})

test_that("frailty_sensitivity errors clearly for unsupported fits", {
  d <- .frailty_data()
  fit_gmm <- quick_fit(d, method = "poisson", gamma = 0.005,
                       estimator = "gmm")
  fit_ols <- quick_fit(d, method = "ols", gamma = 0.005)
  fit_logit <- quick_fit(d, method = "poisson", gamma = 0.005,
                         link_rho = "logit")
  fit_constrained <- quick_fit(d, method = "poisson", gamma = 0.005,
                               constrained = TRUE)
  fit_cov_gamma <- estimate_hidden_pop(
    data = cov_gamma_data(),
    observed = ~ m,
    auxiliary = ~ n,
    reference_pop = ~ N,
    method = "poisson",
    gamma = "estimate",
    cov_gamma = ~ z
  )

  expect_error(frailty_sensitivity(fit_gmm, plot = FALSE), "estimator = 'mle'")
  expect_error(frailty_sensitivity(fit_ols, plot = FALSE), "Poisson and NB")
  expect_error(frailty_sensitivity(fit_logit, plot = FALSE), "link_rho = 'power'")
  expect_error(frailty_sensitivity(fit_constrained, plot = FALSE), "constrained")
  expect_error(frailty_sensitivity(fit_cov_gamma, plot = FALSE), "cov_gamma")
})
