# ---- iOLS extended tests (skip on CRAN) ----

skip_on_cran()

# ---- Full data: basic fit and convergence ----

test_that("iOLS on full testdata converges", {
  skip_on_cran()
  fit <- quick_fit(testdata, method = "iols", gamma = 0.005)
  expect_equal(fit$convergence, 0L)
  expect_true(all(is.finite(coef(fit))))
})

test_that("iOLS on full testdata with covariates", {
  skip_on_cran()
  fit <- quick_fit(testdata, method = "iols", gamma = 0.005,
                   cov_alpha = ~0 + sex, cov_beta = ~0 + sex)
  expect_equal(length(coef(fit)), 4)
  expect_equal(fit$convergence, 0L)
})

test_that("iOLS popsize on full data: all positive with BC", {
  skip_on_cran()
  fit <- quick_fit(testdata, method = "iols", gamma = 0.005)
  ps <- popsize(fit, bias_correction = TRUE)
  expect_true(all(ps$estimate > 0))
  expect_true(all(ps$estimate_bc > 0))
})

test_that("iOLS popsize by sex on full data", {
  skip_on_cran()
  fit <- quick_fit(testdata, method = "iols", gamma = 0.005)
  ps <- popsize(fit, by = ~sex)
  expect_equal(nrow(ps), 2)
  expect_true(all(ps$estimate > 0))
})

# ---- Irregular migration data (the hard specification) ----

test_that("iOLS converges on year*ukr specification", {
  skip_on_cran()
  data(irregular_migration)
  d <- irregular_migration
  d$ukr <- as.integer(d$country_code == "UKR")
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "iols",
    cov_alpha = ~ year * ukr, cov_beta = ~ year,
    gamma = 0.005, countries = ~country_code)
  expect_equal(fit$convergence, 0L)
  # Normalized GPML score should be small
  mu <- fit$fitted.values
  u <- fit$m / pmax(mu, 1e-10)
  score <- crossprod(fit$model_matrix_full, u - 1) / nrow(d)
  expect_true(max(abs(score)) < 1e-4)
})

test_that("iOLS popsize on year*ukr: all BC positive", {
  skip_on_cran()
  data(irregular_migration)
  d <- irregular_migration
  d$ukr <- as.integer(d$country_code == "UKR")
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "iols",
    cov_alpha = ~ year * ukr, cov_beta = ~ year,
    gamma = 0.005, countries = ~country_code)
  ps <- popsize(fit, bias_correction = TRUE)
  expect_true(all(ps$estimate_bc > 0),
              info = paste("Negative BC:", paste(ps$group[ps$estimate_bc <= 0], collapse = ", ")))
})

test_that("iOLS with year*ukr+sex specification", {
  skip_on_cran()
  data(irregular_migration)
  d <- irregular_migration
  d$ukr <- as.integer(d$country_code == "UKR")
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "iols",
    cov_alpha = ~ year * ukr + sex, cov_beta = ~ year,
    gamma = 0.005, countries = ~country_code)
  expect_equal(fit$convergence, 0L)
  ps <- popsize(fit, bias_correction = TRUE)
  expect_true(all(ps$estimate_bc > 0))
})

# ---- Cluster-robust SEs for iOLS ----

test_that("iOLS with cluster-robust SEs", {
  skip_on_cran()
  fit_hc <- quick_fit(testdata, method = "iols", gamma = 0.005, vcov = "HC1")
  fit_cl <- quick_fit(testdata, method = "iols", gamma = 0.005,
                      vcov = "HC1", cluster = ~country)
  # Coefficients identical (clustering only affects SEs)
  expect_equal(coef(fit_hc), coef(fit_cl))
  # SEs should differ
  se_hc <- sqrt(diag(vcov(fit_hc)))
  se_cl <- sqrt(diag(vcov(fit_cl)))
  expect_false(all(abs(se_hc - se_cl) < 1e-10))
})

# ---- Bootstrap with iOLS ----

test_that("bootstrap works with iOLS", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  fit <- quick_fit(testdata, method = "iols", gamma = 0.005)
  boot <- bootstrap_popsize(fit, R = 29, seed = 42, verbose = FALSE)
  expect_s3_class(boot, "uncounted_boot")
  expect_true(boot$n_converged > 10)
  expect_true(all(is.finite(boot$popsize$estimate)))
})

# ---- LOO with iOLS ----

test_that("LOO by country works with iOLS", {
  skip_on_cran()
  fit <- quick_fit(testdata, method = "iols", gamma = 0.005,
                   countries = ~country)
  loo_res <- loo(fit, by = "country")
  expect_s3_class(loo_res, "uncounted_loo")
  expect_equal(loo_res$n_drops, length(unique(testdata$country)))
  expect_true(sum(loo_res$converged) > 10)
})

# ---- dfbeta / dfpopsize with iOLS ----

test_that("dfbeta works with iOLS on full data", {
  skip_on_cran()
  fit <- quick_fit(testdata, method = "iols", gamma = 0.005,
                   countries = ~country)
  db <- dfbeta(fit, by = "country")
  expect_equal(nrow(db), length(unique(testdata$country)))
  expect_equal(ncol(db), length(coef(fit)))
})

test_that("dfpopsize works with iOLS on full data", {
  skip_on_cran()
  fit <- quick_fit(testdata, method = "iols", gamma = 0.005,
                   countries = ~country)
  dp <- dfpopsize(fit, by = "country")
  expect_equal(length(dp), length(unique(testdata$country)))
  expect_true(sum(is.finite(dp)) > 10)
})

# ---- HC type coverage ----

test_that("iOLS supports all HC types", {
  skip_on_cran()
  d <- small_data()
  for (hc in c("HC0", "HC1", "HC2", "HC3", "HC4", "HC5")) {
    fit <- quick_fit(d, method = "iols", gamma = 0.005, vcov = hc)
    se <- sqrt(diag(vcov(fit)))
    expect_true(all(is.finite(se)), info = hc)
    expect_true(all(se > 0), info = hc)
  }
})

# ---- iOLS vs Poisson numerical comparison on full data ----

test_that("iOLS coefficients same sign as Poisson on full data", {
  skip_on_cran()
  fit_po <- quick_fit(testdata, gamma = 0.005)
  fit_io <- quick_fit(testdata, method = "iols", gamma = 0.005)
  expect_equal(sign(coef(fit_po)), sign(coef(fit_io)))
})

test_that("iOLS and Poisson fitted values positively correlated", {
  skip_on_cran()
  fit_po <- quick_fit(testdata, gamma = 0.005)
  fit_io <- quick_fit(testdata, method = "iols", gamma = 0.005)
  r <- cor(fitted(fit_po), fitted(fit_io))
  expect_true(r > 0.9,
              info = paste("Fitted correlation:", round(r, 4)))
})
