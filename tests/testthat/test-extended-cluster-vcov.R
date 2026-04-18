# ---- Cluster-robust vcov on real data (extended, skip on CRAN) ----

skip_on_cran()

test_that("cluster-robust SEs differ from observation-level HC3", {
  skip_on_cran()
  data(irregular_migration)
  d <- irregular_migration


  fit_no_cl <- estimate_hidden_pop(
    data = d, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
    method = "poisson", cov_alpha = ~ year + sex, cov_beta = ~ year,
    vcov = "HC3"
  )

  fit_cl <- estimate_hidden_pop(
    data = d, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
    method = "poisson", cov_alpha = ~ year + sex, cov_beta = ~ year,
    vcov = "HC1", cluster = ~ country_code
  )

  expect_equal(coef(fit_no_cl), coef(fit_cl))
  se_no_cl <- sqrt(diag(vcov(fit_no_cl)))
  se_cl <- sqrt(diag(vcov(fit_cl)))
  expect_false(all(abs(se_no_cl - se_cl) < 1e-10))
})

test_that("cluster-robust SEs work for NB", {
  skip_on_cran()
  data(irregular_migration)
  d <- irregular_migration


  fit_cl <- estimate_hidden_pop(
    data = d, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
    method = "nb", countries = ~ country_code, cluster = ~ country_code
  )

  expect_s3_class(fit_cl, "uncounted")
  expect_true(!is.null(fit_cl$cluster_var))
  expect_true(all(is.finite(sqrt(diag(vcov(fit_cl))))))
})

test_that("vcov label shows CL when clustered", {
  skip_on_cran()
  data(irregular_migration)
  d <- irregular_migration


  fit_cl <- estimate_hidden_pop(
    data = d, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
    method = "poisson", cluster = ~ country_code
  )

  out <- capture.output(print(fit_cl))
  expect_true(any(grepl("CL\\(", out)))
})

test_that("countries does NOT trigger clustering", {
  skip_on_cran()
  data(irregular_migration)
  d <- irregular_migration


  fit <- estimate_hidden_pop(
    data = d, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
    method = "poisson", countries = ~ country_code
  )

  out <- capture.output(print(fit))
  expect_false(any(grepl("CL\\(", out)))
  expect_true(any(grepl("HC3", out)))
})
