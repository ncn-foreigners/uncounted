# ---- Core estimation tests (CRAN) ----

# ---- All methods run ----

test_that("all methods run with fixed gamma", {
  d <- small_data()
  for (method in c("ols", "poisson", "nb", "nls")) {
    fit <- quick_fit(d, method = method, gamma = 0.005)
    expect_s3_class(fit, "uncounted")
    expect_true(length(coef(fit)) >= 2)
    expect_true(all(is.finite(coef(fit))))
  }
})

test_that("all methods run with estimated gamma", {
  d <- small_data()
  for (method in c("ols", "poisson", "nb")) {
    fit <- quick_fit(d, method = method, gamma = "estimate")
    expect_s3_class(fit, "uncounted")
    expect_true(!is.null(fit$gamma))
    expect_true(fit$gamma > 0)
  }
})

test_that("methods handle zeros with gamma", {
  d <- small_data()  # has zeros
  for (method in c("ols", "poisson", "nb")) {
    fit <- quick_fit(d, method = method, gamma = 0.005)
    expect_s3_class(fit, "uncounted")
    expect_true(all(is.finite(coef(fit))))
  }
})

# ---- Covariates ----

test_that("covariates in alpha work", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, cov_alpha = ~0 + sex)
  expect_true("alpha:sexF" %in% names(coef(fit)))
  expect_true("alpha:sexM" %in% names(coef(fit)))
})

test_that("covariates in both alpha and beta work", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005,
                   cov_alpha = ~0 + sex, cov_beta = ~0 + sex)
  expect_equal(length(coef(fit)), 4)
})

test_that("~0 formula treated as intercept-only", {
  d <- small_data()
  fit1 <- quick_fit(d, gamma = 0.005)
  fit2 <- quick_fit(d, gamma = 0.005, cov_alpha = ~0)
  expect_equal(as.numeric(coef(fit1)), as.numeric(coef(fit2)))
})

# ---- S3 methods ----

test_that("S3 methods return correct types", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)

  expect_true(is.numeric(coef(fit)))
  expect_true(is.matrix(vcov(fit)))
  expect_equal(nrow(vcov(fit)), length(coef(fit)))
  expect_true(is.numeric(fitted(fit)))
  expect_equal(length(fitted(fit)), nrow(d))
  expect_true(is.numeric(residuals(fit)))
  expect_equal(length(residuals(fit)), nrow(d))
})

test_that("print and summary run without error", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  expect_output(print(fit))
  expect_output(summary(fit))
})

# ---- vcov types ----

test_that("HC0 < HC1 and HC2 < HC3 SE ordering", {
  d <- small_data()
  se_list <- lapply(c("HC0", "HC1", "HC2", "HC3"), function(hc) {
    fit <- quick_fit(d, gamma = 0.005, vcov = hc)
    sqrt(diag(vcov(fit)))
  })
  expect_true(all(se_list[[1]] < se_list[[2]]))
  expect_true(all(se_list[[3]] < se_list[[4]]))
})

test_that("HC4, HC4m, HC5 produce valid SEs", {
  d <- small_data()
  for (hc in c("HC4", "HC4m", "HC5")) {
    fit <- quick_fit(d, gamma = 0.005, vcov = hc)
    se <- sqrt(diag(vcov(fit)))
    expect_true(all(is.finite(se)), info = hc)
    expect_true(all(se > 0), info = hc)
  }
})

test_that("model-based vcov is present", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, vcov = "HC3")
  expect_true(!is.null(fit$vcov_model))
  expect_true(is.matrix(fit$vcov_model))
  expect_equal(dim(fit$vcov_model), dim(fit$vcov))
})

# ---- NB extras ----

test_that("NB returns theta", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005)
  expect_true(!is.null(fit$theta))
  expect_true(fit$theta > 0)
})

test_that("NB log-likelihood >= Poisson (nested model)", {
  d <- testdata  # need more data for stable NB
  fit_p  <- quick_fit(d, method = "poisson", gamma = 0.005)
  fit_nb <- quick_fit(d, method = "nb", gamma = 0.005)
  expect_true(fit_nb$loglik >= fit_p$loglik)
})

# ---- OLS gamma profiling ----

test_that("OLS gamma='estimate' produces gamma in (0,1)", {
  d <- small_data()
  fit <- quick_fit(d, method = "ols", gamma = "estimate")
  expect_true(!is.null(fit$gamma))
  expect_true(fit$gamma > 0 & fit$gamma < 1)
  expect_true(fit$gamma_estimated)
})

# ---- Edge cases ----

test_that("invalid method gives error", {
  d <- small_data()
  expect_error(quick_fit(d, method = "invalid"))
})

test_that("no warning when data has no zeros and gamma=NULL", {
  d <- positive_data(small_data())
  expect_no_warning(quick_fit(d, gamma = NULL))
})

test_that("zeros in n with gamma are handled", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  expect_s3_class(fit, "uncounted")
  expect_true(all(is.finite(coef(fit))))
})

# ---- weights and update ----

test_that("weights method returns NULL when no weights given", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005)
  expect_null(weights(fit))
})

test_that("weights method returns weights when given", {
  d <- small_data()
  w <- rep(1, nrow(d))
  fit <- quick_fit(d, gamma = 0.005, weights = w)
  expect_equal(weights(fit), w)
})

test_that("update changes vcov type", {
  d <- small_data()
  fit_hc3 <- estimate_hidden_pop(d, ~m, ~n, ~N,
                                  method = "poisson", gamma = 0.005)
  fit_hc0 <- update(fit_hc3, vcov = "HC0")
  expect_equal(coef(fit_hc3), coef(fit_hc0))
  expect_false(all(abs(sqrt(diag(vcov(fit_hc3))) -
                       sqrt(diag(vcov(fit_hc0)))) < 1e-10))
})

test_that("update with evaluate=FALSE returns a call", {
  d <- small_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N,
                              method = "poisson", gamma = 0.005)
  cl <- update(fit, vcov = "HC0", evaluate = FALSE)
  expect_true(is.call(cl))
})

test_that("update with new weights refits the model", {
  d <- small_data()
  fit1 <- estimate_hidden_pop(d, ~m, ~n, ~N,
                               method = "poisson", gamma = 0.005)
  w <- runif(nrow(d), 0.5, 2)
  fit2 <- update(fit1, weights = w)
  expect_false(all(abs(coef(fit1) - coef(fit2)) < 1e-10))
  expect_equal(weights(fit2), w)
})

# ---- LOO (small) ----

test_that("loo runs and returns correct structure", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, countries = ~country)
  loo_res <- loo(fit, by = "obs")
  expect_s3_class(loo_res, "uncounted_loo")
  expect_equal(loo_res$n_drops, nrow(d))
  expect_true(sum(loo_res$converged) > 0)
})

test_that("loo by country works", {
  d <- small_data()
  fit <- quick_fit(d, gamma = 0.005, countries = ~country)
  loo_res <- loo(fit, by = "country")
  expect_equal(loo_res$n_drops, length(unique(d$country)))
})

# ---- LOO rank-deficiency detection ----

test_that("loo by country warns and skips when dropping removes covariate level", {
  d <- small_data()
  # Add a binary indicator for a single country
  d$is_special <- as.integer(d$country == d$country[1])
  fit <- quick_fit(d, gamma = 0.005, countries = ~country,
                   cov_alpha = ~is_special)
  expect_warning(
    loo_res <- loo(fit, by = "country"),
    "rank-deficient"
  )
  # The special country should be marked as not converged
  special_idx <- which(loo_res$dropped == d$country[1])
  expect_false(loo_res$converged[special_idx])
  # Other countries should still converge
  expect_true(sum(loo_res$converged) > 0)
})
