# ---- Bayesian brms support (CRAN-safe tests) ----

.fake_bayes_fit <- function(constrained = FALSE) {
  d <- data.frame(
    year = factor(c("2019", "2019", "2024", "2024")),
    sex = factor(c("F", "M", "F", "M")),
    m = c(1L, 2L, 3L, 4L),
    n = c(1, 1, 2, 2),
    N = c(100, 120, 140, 160)
  )
  X_alpha <- model.matrix(~ year, data = d)
  colnames(X_alpha) <- paste0("alpha:", colnames(X_alpha))
  X_beta <- model.matrix(~ 1, data = d)
  colnames(X_beta) <- "beta"

  alpha_draws <- matrix(
    c(
      0.50, 0.50, 0.65, 0.65,
      0.55, 0.55, 0.70, 0.70,
      0.60, 0.60, 0.75, 0.75,
      0.65, 0.65, 0.80, 0.80
    ),
    nrow = 4, byrow = TRUE
  )

  out <- list(
    fit = NULL,
    engine = "test",
    backend = "test",
    method = "poisson",
    estimator = "bayes",
    gamma = 0.005,
    gamma_estimated = FALSE,
    constrained = constrained,
    link_rho = "power",
    data = d,
    m = d$m,
    N = d$N,
    n_aux = d$n,
    X_alpha = X_alpha,
    X_beta = X_beta,
    cov_alpha_vars = d["year"],
    countries_var = NULL,
    obs_weights = NULL,
    n_obs = nrow(d),
    alpha_draws = alpha_draws,
    coef_draws = cbind(
      b_alpha_Intercept = c(0.50, 0.55, 0.60, 0.65),
      b_alpha_year2024 = c(0.15, 0.15, 0.15, 0.15),
      b_beta_Intercept = c(0.80, 0.82, 0.84, 0.86)
    ),
    fitted_draws = matrix(
      c(
        1.0, 2.0, 3.0, 4.0,
        1.5, 2.5, 3.5, 4.5,
        2.0, 3.0, 4.0, 5.0,
        2.5, 3.5, 4.5, 5.5
      ),
      nrow = 4, byrow = TRUE
    )
  )
  class(out) <- "uncounted_bayes"
  out
}

test_that("Bayesian summary and coefficient methods use posterior draws", {
  skip_if_not_installed("posterior")
  fit <- .fake_bayes_fit()

  coefs <- coef(fit)
  expect_named(coefs, c("alpha:(Intercept)", "alpha:year2024",
                        "beta:(Intercept)"))
  expect_equal(unname(coefs["alpha:(Intercept)"]), 0.575)

  coef_draws <- coef(fit, summary = "draws")
  expect_equal(dim(coef_draws), c(4, 3))

  out <- capture.output(s <- summary(fit, total = TRUE))
  expect_true(any(grepl("Bayesian unauthorized population estimation", out)))
  expect_s3_class(s, "summary_uncounted_bayes")
  expect_true(all(c("term", "median", "sd", "rhat", "ess_bulk") %in%
                    names(s$coefficients)))
  expect_s3_class(s$popsize, "uncounted_popsize_bayes")
})

test_that("Bayesian fitted, predict, and residual methods summarize draws", {
  fit <- .fake_bayes_fit()

  expect_equal(fitted(fit, summary = "mean"), c(1.75, 2.75, 3.75, 4.75))
  expect_equal(fitted(fit, summary = "median"), c(1.75, 2.75, 3.75, 4.75))
  expect_equal(predict(fit, type = "response"), c(1.75, 2.75, 3.75, 4.75))
  expect_equal(predict(fit, type = "draws"), fit$fitted_draws)

  res <- residuals(fit)
  expect_equal(res, fit$m - c(1.75, 2.75, 3.75, 4.75))
  pearson <- residuals(fit, type = "pearson")
  expect_true(all(is.finite(pearson)))
})

test_that("Bayesian tidy and glance methods return data frames", {
  skip_if_not_installed("posterior")
  fit <- .fake_bayes_fit()

  td <- tidy(fit, conf.int = TRUE)
  expect_s3_class(td, "data.frame")
  expect_true(all(c("term", "estimate", "std.error", "conf.low",
                    "conf.high") %in% names(td)))
  expect_true(all(is.na(td$p.value)))

  gl <- glance(fit)
  expect_s3_class(gl, "data.frame")
  expect_equal(gl$estimator, "BAYES")
  expect_equal(gl$nobs, fit$n_obs)
  expect_true("max_rhat" %in% names(gl))
})

test_that("popsize() summarizes deterministic Bayesian posterior draws", {
  fit <- .fake_bayes_fit()
  ps <- popsize(fit, by = ~ year, total = TRUE)

  expect_s3_class(ps, "uncounted_popsize_bayes")
  expect_equal(nrow(ps), 2)
  expect_true(all(c("estimate", "median", "mean", "sd", "lower", "upper") %in%
                    names(ps)))

  draws <- posterior_popsize_draws(ps, include_total = TRUE)
  expect_equal(dim(draws), c(4, 3))
  expect_equal(colnames(draws), c("2019", "2024", "Total"))
  expect_equal(unname(ps$median),
               unname(apply(draws[, 1:2, drop = FALSE], 2, stats::median)))
  expect_equal(attr(ps, "total")$median, stats::median(draws[, "Total"]))
})

test_that("Bayesian hypotheses report posterior proposition probabilities", {
  fit <- .fake_bayes_fit()
  ps <- popsize(fit, by = ~ year)

  hyp <- hypotheses_popsize(
    ps,
    "xi[year == 2024] - xi[year == 2019] > 0"
  )

  expect_equal(hyp$method, "Bayesian posterior probability")
  expect_equal(hyp$alternative, "greater")
  expect_equal(hyp$p_h1, 1)
  expect_equal(hyp$p_h0, 0)
  expect_true(is.infinite(hyp$posterior_odds))
  expect_true(is.na(hyp$p.value))
})

test_that("hypothesis_side='null' reverses Bayesian posterior propositions", {
  fit <- .fake_bayes_fit()
  ps <- popsize(fit, by = ~ year)

  hyp <- hypotheses_popsize(
    ps,
    "xi[year == 2024] - xi[year == 2019] > 0",
    hypothesis_side = "null"
  )

  expect_equal(hyp$alternative, "less")
  expect_equal(hyp$p_h1, 0)
  expect_equal(hyp$p_h0, 1)
  expect_equal(hyp$null_hypothesis,
               "xi[year == 2024] - xi[year == 2019] >= 0")
  expect_equal(hyp$alternative_hypothesis,
               "xi[year == 2024] - xi[year == 2019] < 0")
})

test_that("Bayesian hypotheses support functions, numeric nulls, and ROPE", {
  fit <- .fake_bayes_fit()
  ps <- popsize(fit, by = ~ year)

  numeric_hyp <- hypotheses_popsize(ps, 10)
  expect_equal(nrow(numeric_hyp), 2)
  expect_true(all(is.na(numeric_hyp$p_h1)))

  fun_hyp <- hypotheses_popsize(
    ps,
    function(x, groups) c(diff = unname(x["2024"] - x["2019"])),
    rope = c(-1, 1)
  )
  expect_equal(fun_hyp$hypothesis, "diff")
  expect_true(is.finite(fun_hyp$std.error))
  expect_equal(fun_hyp$rope.low, -1)
  expect_equal(fun_hyp$rope.high, 1)
  expect_true(fun_hyp$p_rope >= 0 && fun_hyp$p_rope <= 1)
})

test_that("Bayesian brms specs generate expected Stan likelihood code", {
  skip_if_not_installed("brms")

  countries <- unique(testdata$country)[1:4]
  d <- testdata[testdata$country %in% countries &
                  testdata$year %in% c("2019", "2023"), ]
  d <- droplevels(d)
  spec <- uncounted:::.bayes_make_brms_spec(
    data = d,
    observed = ~ m,
    auxiliary = ~ n,
    reference_pop = ~ N,
    method = "poisson",
    cov_alpha = ~ year + sex,
    cov_beta = ~ year,
    gamma = 0.005,
    link_rho = "power"
  )
  code <- brms::make_stancode(spec$formula, data = spec$data,
                              family = spec$family, prior = spec$prior)
  expect_match(code, "poisson_log", fixed = TRUE)

  spec_nb <- uncounted:::.bayes_make_brms_spec(
    data = d,
    observed = ~ m,
    auxiliary = ~ n,
    reference_pop = ~ N,
    method = "nb",
    cov_alpha = ~ year + sex,
    cov_beta = ~ year,
    gamma = 0.005,
    link_rho = "power"
  )
  code_nb <- brms::make_stancode(spec_nb$formula, data = spec_nb$data,
                                 family = spec_nb$family,
                                 prior = spec_nb$prior)
  expect_match(code_nb, "neg_binomial_2_log", fixed = TRUE)
})

test_that("estimated gamma and bounded links are represented in brms code", {
  skip_if_not_installed("brms")

  countries <- unique(testdata$country)[1:4]
  d <- testdata[testdata$country %in% countries &
                  testdata$year %in% c("2019", "2023"), ]
  d <- droplevels(d)
  spec <- uncounted:::.bayes_make_brms_spec(
    data = d,
    observed = ~ m,
    auxiliary = ~ n,
    reference_pop = ~ N,
    method = "poisson",
    cov_alpha = ~ year + sex,
    cov_beta = ~ year,
    gamma = "estimate",
    link_rho = "cloglog",
    constrained = TRUE
  )
  code <- brms::make_stancode(spec$formula, data = spec$data,
                              family = spec$family, prior = spec$prior,
                              stanvars = spec$stanvars)

  expect_match(code, "uncounted_log_cloglog", fixed = TRUE)
  expect_match(code, "inv_logit", fixed = TRUE)
  expect_match(code, "K_gamma", fixed = TRUE)
})

test_that("Bayesian sampling smoke test can be enabled explicitly", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("rstan")
  skip_if(Sys.getenv("UNCOUNTED_RUN_STAN_TESTS") != "true",
          "Set UNCOUNTED_RUN_STAN_TESTS=true to run Stan smoke tests.")

  d <- small_data()[1:8, ]
  fit <- estimate_hidden_pop_bayes(
    data = d,
    observed = ~ m,
    auxiliary = ~ n,
    reference_pop = ~ N,
    method = "poisson",
    gamma = 0.005,
    chains = 1,
    iter = 200,
    warmup = 100,
    seed = 1,
    refresh = 0
  )
  ps <- popsize(fit)
  expect_s3_class(ps, "uncounted_popsize_bayes")
  expect_true(all(is.finite(ps$estimate)))
})
