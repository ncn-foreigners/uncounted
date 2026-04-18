# ---- Dependence profile tests (CRAN) ----

.dependence_profile_data <- function() {
  d <- positive_data()
  keep_country <- unique(d$country)[1:8]
  keep_year <- levels(d$year)[1:2]
  droplevels(d[d$country %in% keep_country & d$year %in% keep_year, ])
}

test_that("profile_dependence returns the expected structure", {
  d <- .dependence_profile_data()
  fit <- quick_fit(d, method = "poisson", gamma = 0.005,
                   cov_alpha = ~ year + sex, cov_beta = ~ year)

  prof <- profile_dependence(fit, delta_grid = c(-0.2, 0, 0.2), plot = FALSE)

  expect_s3_class(prof, "data.frame")
  expect_equal(names(prof), c("delta", "kappa", "xi", "loglik"))
  expect_equal(nrow(prof), 3)
  expect_equal(prof$kappa, exp(prof$delta), tolerance = 1e-12)
})

test_that("profile_dependence reproduces baseline xi at delta = 0", {
  d <- .dependence_profile_data()
  fit <- quick_fit(d, method = "poisson", gamma = 0.005,
                   cov_alpha = ~ year + sex, cov_beta = ~ year)

  prof <- profile_dependence(fit, delta_grid = c(-0.1, 0, 0.1), plot = FALSE)
  baseline <- sum(popsize(fit)$estimate)

  expect_equal(prof$xi[prof$delta == 0], baseline, tolerance = 1e-8)
})

test_that("profile_dependence works for Poisson and NB fits", {
  d <- .dependence_profile_data()

  for (method in c("poisson", "nb")) {
    fit <- quick_fit(d, method = method, gamma = 0.005,
                     cov_alpha = ~ sex, cov_beta = ~ year)
    prof <- profile_dependence(fit, delta_grid = c(-0.1, 0, 0.1), plot = FALSE)

    expect_equal(nrow(prof), 3, info = method)
    expect_true(all(is.finite(prof$xi)), info = method)
    expect_true(all(is.finite(prof$loglik)), info = method)
  }
})

test_that("profile_dependence errors clearly for unsupported fits", {
  fit_ols <- quick_fit(small_data(), method = "ols", gamma = 0.005)
  fit_gmm <- quick_fit(small_data(), method = "poisson", gamma = 0.005,
                       estimator = "gmm")

  expect_error(profile_dependence(fit_ols, plot = FALSE), "Poisson and NB")
  expect_error(profile_dependence(fit_gmm, plot = FALSE), "estimator = 'mle'")
})

test_that("profile.uncounted dispatches to profile_dependence", {
  d <- .dependence_profile_data()
  fit <- quick_fit(d, method = "poisson", gamma = 0.005,
                   cov_alpha = ~ sex, cov_beta = ~ year)

  prof_direct <- profile_dependence(fit, delta_grid = c(-0.1, 0, 0.1), plot = FALSE)
  prof_s3 <- profile(fit, param = "dependence",
                     delta_grid = c(-0.1, 0, 0.1), plot = FALSE)

  expect_equal(prof_s3, prof_direct, tolerance = 1e-8)
})

test_that("profile_dependence with plot = TRUE runs without error", {
  fit <- quick_fit(small_data(), method = "poisson", gamma = 0.005)

  expect_no_error(
    profile_dependence(fit, delta_grid = c(-0.1, 0, 0.1), plot = TRUE)
  )
})
