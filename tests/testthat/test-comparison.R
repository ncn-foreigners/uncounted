make_test_data <- function(n = 50, seed = 42) {
  set.seed(seed)
  N <- round(exp(rnorm(n, 5, 1.5)))
  n_aux <- rpois(n, exp(-3 + 0.8 * log(N)))
  ratio <- n_aux / N
  mu <- N^0.7 * (0.005 + ratio)^0.5
  m <- rpois(n, mu)
  data.frame(m = m, n = n_aux, N = N)
}

test_that("compare_models works with 2 models", {
  dt <- make_test_data()
  fit_po <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  fit_nb <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  comp <- compare_models(fit_po, fit_nb)
  expect_s3_class(comp, "uncounted_comparison")
  expect_equal(nrow(comp$table), 2)
  expect_true(all(c("AIC", "BIC", "logLik") %in% names(comp$table)))
})

test_that("compare_models includes pseudo R^2 columns", {
  dt <- make_test_data()
  fit_po <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  fit_nb <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  comp <- compare_models(fit_po, fit_nb)
  expect_true(all(c("R2_cor", "R2_D", "R2_CW") %in% names(comp$table)))
  # R^2 should be between 0 and 1
  expect_true(all(comp$table$R2_cor >= 0 & comp$table$R2_cor <= 1))
})

test_that("pseudo R^2_D is NA for OLS models", {
  dt <- make_test_data()
  dt_pos <- dt[dt$m > 0 & dt$n > 0, ]
  fit_ols <- estimate_hidden_pop(dt_pos, ~m, ~n, ~N, method = "ols", gamma = NULL)
  fit_po <- estimate_hidden_pop(dt_pos, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  comp <- suppressWarnings(compare_models(OLS = fit_ols, Poisson = fit_po))
  ols_row <- comp$table[comp$table$Model == "OLS", ]
  expect_true(is.na(ols_row$R2_D))
  expect_true(is.na(ols_row$R2_CW))
})

test_that("R^2_D uses null model without covariates", {
  dt <- make_test_data(n = 60)
  dt$grp <- factor(rep(1:2, each = 30))
  fit_cov <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson",
                                  gamma = 0.005, cov_alpha = ~grp)
  fit_null <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson",
                                   gamma = 0.005)
  comp <- compare_models(WithCov = fit_cov, NoCov = fit_null)
  # Model with covariates should have R^2_D >= 0 (or very close)
  expect_true(comp$table$R2_D[comp$table$Model == "WithCov"] >= -0.01)
  # Null model should have R^2_D = 0 (it IS the null)
  expect_true(abs(comp$table$R2_D[comp$table$Model == "NoCov"]) < 0.01)
})

test_that("compare_models prints without error", {
  dt <- make_test_data()
  fit_po <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  fit_nb <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  expect_output(print(compare_models(fit_po, fit_nb)), "Model comparison")
})

test_that("compare_models warns with mixed types", {
  dt <- make_test_data()
  dt_pos <- dt[dt$m > 0 & dt$n > 0, ]
  fit_ols <- estimate_hidden_pop(dt_pos, ~m, ~n, ~N, method = "ols", gamma = NULL)
  fit_po <- estimate_hidden_pop(dt_pos, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  expect_warning(compare_models(fit_ols, fit_po), "pseudo-loglik")
})

test_that("lrtest works for Poisson vs NB", {
  dt <- make_test_data()
  fit_po <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  fit_nb <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  lr <- lrtest(fit_po, fit_nb)
  expect_s3_class(lr, "uncounted_lrtest")
  expect_true(lr$statistic >= 0)
  expect_true(lr$p_value >= 0 && lr$p_value <= 1)
  expect_true(lr$boundary)  # Poisson vs NB is boundary test
})

test_that("lrtest prints without error", {
  dt <- make_test_data()
  fit_po <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  fit_nb <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  expect_output(print(lrtest(fit_po, fit_nb)), "Likelihood ratio test")
})

test_that("lrtest errors with OLS models", {
  dt <- make_test_data()
  dt_pos <- dt[dt$m > 0 & dt$n > 0, ]
  fit1 <- estimate_hidden_pop(dt_pos, ~m, ~n, ~N, method = "ols", gamma = NULL)
  fit2 <- estimate_hidden_pop(dt_pos, ~m, ~n, ~N, method = "ols", gamma = NULL)
  expect_error(lrtest(fit1, fit2), "log-likelihood")
})


# ---- compare_loo tests ----

test_that("compare_loo returns correct class", {
  dt <- make_test_data(n = 30, seed = 42)
  fit_po <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  fit_nb <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  loo_po <- loo(fit_po, by = "obs")
  loo_nb <- loo(fit_nb, by = "obs")
  comp <- compare_loo(loo_po, loo_nb, labels = c("Poisson", "NB"))
  expect_s3_class(comp, "uncounted_loo_compare")
  expect_true("table" %in% names(comp))
  expect_equal(nrow(comp$table), 30)
})

test_that("compare_loo prints without error", {
  dt <- make_test_data(n = 30, seed = 42)
  fit_po <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  fit_nb <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  loo_po <- loo(fit_po, by = "obs")
  loo_nb <- loo(fit_nb, by = "obs")
  comp <- compare_loo(loo_po, loo_nb, labels = c("Poisson", "NB"))
  expect_output(print(comp), "LOO comparison")
})

test_that("compare_loo scatter plot runs without error", {
  dt <- make_test_data(n = 30, seed = 42)
  fit_po <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  fit_nb <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  loo_po <- loo(fit_po, by = "obs")
  loo_nb <- loo(fit_nb, by = "obs")
  comp <- compare_loo(loo_po, loo_nb)
  expect_no_error(plot(comp, type = "scatter"))
})

test_that("compare_loo bar plot runs without error", {
  dt <- make_test_data(n = 30, seed = 42)
  fit_po <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  fit_nb <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  loo_po <- loo(fit_po, by = "obs")
  loo_nb <- loo(fit_nb, by = "obs")
  comp <- compare_loo(loo_po, loo_nb)
  expect_no_error(plot(comp, type = "bar"))
})

test_that("compare_loo errors with mismatched by", {
  dt <- make_test_data(n = 30, seed = 42)
  dt$country <- rep(letters[1:6], each = 5)
  fit_po <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson",
                                 gamma = 0.005, countries = ~country)
  fit_nb <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb",
                                 gamma = 0.005, countries = ~country)
  loo_po <- loo(fit_po, by = "obs")
  loo_nb <- loo(fit_nb, by = "country")
  expect_error(compare_loo(loo_po, loo_nb), "same 'by'")
})

test_that("compare_loo table is sorted by max influence", {
  dt <- make_test_data(n = 30, seed = 42)
  fit_po <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "poisson", gamma = 0.005)
  fit_nb <- estimate_hidden_pop(dt, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  loo_po <- loo(fit_po, by = "obs")
  loo_nb <- loo(fit_nb, by = "obs")
  comp <- compare_loo(loo_po, loo_nb)
  max_abs <- comp$table$max_abs_pct
  expect_true(all(diff(max_abs) <= 0))  # sorted descending
})
