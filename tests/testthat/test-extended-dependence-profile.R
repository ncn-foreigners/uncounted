# ---- Dependence profile tests (extended) ----

test_that("profile_dependence can move xi away from the baseline", {
  d <- positive_data()
  keep_country <- unique(d$country)[1:10]
  keep_year <- levels(d$year)[1:2]
  d <- droplevels(d[d$country %in% keep_country & d$year %in% keep_year, ])

  fit <- quick_fit(d, method = "poisson", gamma = 0.005,
                   cov_alpha = ~ year + sex, cov_beta = ~ year)
  prof <- profile_dependence(fit, delta_grid = c(-0.3, 0, 0.3), plot = FALSE)

  baseline <- prof$xi[prof$delta == 0]
  expect_true(any(abs(prof$xi[prof$delta != 0] - baseline) > 1e-6))
})

test_that("profile_dependence supports covariate-rich NB fits", {
  d <- positive_data()
  keep_country <- unique(d$country)[1:8]
  keep_year <- levels(d$year)[1:3]
  d <- droplevels(d[d$country %in% keep_country & d$year %in% keep_year, ])

  fit <- quick_fit(d, method = "nb", gamma = 0.005,
                   cov_alpha = ~ year + sex, cov_beta = ~ year)
  prof <- profile_dependence(fit, delta_grid = seq(-0.2, 0.2, length.out = 5),
                             plot = FALSE)

  expect_true(all(is.finite(prof$xi)))
  expect_true(all(is.finite(prof$loglik)))
})

test_that("profile_dependence works for constrained Poisson fits", {
  fit <- quick_fit(small_data(), method = "poisson", gamma = 0.005,
                   cov_alpha = ~ sex, constrained = TRUE)

  prof <- profile_dependence(fit, delta_grid = c(-0.1, 0, 0.1), plot = FALSE)

  expect_equal(nrow(prof), 3)
  expect_true(all(is.finite(prof$xi)))
})
