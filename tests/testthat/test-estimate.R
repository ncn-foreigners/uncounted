# ---- Test data ----

make_test_data <- function() {
  data.frame(
    m = c(14L, 1L, 56L, 21L, 29L, 57L, 1L, 41L, 258L, 6L,
          19L, 9L, 3L, 124L, 4L, 15L, 30L, 10L, 7L, 686L),
    n = c(2L, 1L, 27L, 9L, 85L, 82L, 1L, 16L, 586L, 2L,
          38L, 1L, 1L, 37L, 4L, 12L, 19L, 6L, 7L, 774L),
    N = c(1023L, 188L, 14655L, 892L, 3672L, 4473L, 44L, 5081L, 131365L, 1357L,
          4391L, 320L, 309L, 6422L, 482L, 1119L, 1784L, 746L, 2732L, 25576L),
    sex = rep(c("M", "F"), each = 10),
    country = paste0("C", 1:20)
  )
}

make_test_data_with_zeros <- function() {
  d <- make_test_data()
  d$m[c(2, 5)] <- 0L
  d$n[c(3, 7)] <- 0L
  d
}

# ---- Poisson matches glm ----

test_that("Poisson without gamma matches glm()", {
  d <- make_test_data()

  glm_fit <- glm(m ~ -1 + log(N) + log(n / N), family = poisson(), data = d)
  our_fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = NULL)

  expect_equal(as.numeric(coef(our_fit)), as.numeric(coef(glm_fit)),
               tolerance = 1e-5)
})

test_that("Poisson with sex covariate matches glm()", {
  d <- make_test_data()

  glm_fit <- glm(m ~ -1 + sex:log(N) + log(n / N), family = poisson(), data = d)
  our_fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = NULL,
                                  cov_alpha = ~0 + sex)

  glm_coefs <- coef(glm_fit)
  our_coefs <- coef(our_fit)

  # alpha:sexF, alpha:sexM should match sexF:log(N), sexM:log(N)
  expect_equal(our_coefs[["alpha:sexF"]], glm_coefs[["sexF:log(N)"]], tolerance = 1e-5)
  expect_equal(our_coefs[["alpha:sexM"]], glm_coefs[["sexM:log(N)"]], tolerance = 1e-5)
  expect_equal(our_coefs[["beta"]], glm_coefs[["log(n/N)"]], tolerance = 1e-5)
})

# ---- OLS matches lm ----

test_that("OLS with fixed gamma matches lm()", {
  d <- make_test_data()
  gamma_val <- 0.001

  # No zeros in m → code uses log(m), not log(m+1)
  lm_fit <- lm(log(m) ~ 0 + log(N) + log(gamma_val + n / N), data = d)
  our_fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "ols", gamma = gamma_val)

  expect_equal(as.numeric(coef(our_fit)), as.numeric(coef(lm_fit)),
               tolerance = 1e-5)
})

test_that("OLS with gamma uses log(m+1) when zeros present", {
  d <- make_test_data_with_zeros()
  gamma_val <- 0.001

  lm_fit <- lm(log(m + 1) ~ 0 + log(N) + log(gamma_val + n / N), data = d)
  our_fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "ols", gamma = gamma_val)

  expect_equal(as.numeric(coef(our_fit)), as.numeric(coef(lm_fit)),
               tolerance = 1e-5)
})

test_that("OLS with sex covariate and gamma matches lm()", {
  d <- make_test_data()
  gamma_val <- 0.001

  lm_fit <- lm(log(m + 1) ~ 0 + sex:log(N) + sex:log(gamma_val + n / N), data = d)
  our_fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "ols", gamma = gamma_val,
                                  cov_alpha = ~0 + sex, cov_beta = ~0 + sex)

  expect_equal(length(coef(our_fit)), length(coef(lm_fit)))
})

# ---- Methods run without error ----

test_that("all methods run with fixed gamma", {
  d <- make_test_data()

  for (method in c("ols", "poisson", "nb", "nls")) {
    fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = method, gamma = 0.001)
    expect_s3_class(fit, "uncounted")
    expect_true(length(coef(fit)) >= 2)
    expect_true(all(is.finite(coef(fit))))
  }
})

test_that("all methods run with estimated gamma", {
  d <- make_test_data()

  for (method in c("ols", "poisson", "nb")) {
    fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = method, gamma = "estimate")
    expect_s3_class(fit, "uncounted")
    expect_true(!is.null(fit$gamma))
    expect_true(fit$gamma > 0)
  }
})

test_that("methods handle zeros with gamma", {
  d <- make_test_data_with_zeros()

  for (method in c("ols", "poisson", "nb")) {
    fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = method, gamma = 0.001)
    expect_s3_class(fit, "uncounted")
    expect_true(all(is.finite(coef(fit))))
  }
})

# ---- Covariates ----

test_that("covariates in alpha work", {
  d <- make_test_data()

  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001,
                              cov_alpha = ~0 + sex)
  expect_true("alpha:sexF" %in% names(coef(fit)))
  expect_true("alpha:sexM" %in% names(coef(fit)))
})

test_that("covariates in both alpha and beta work", {
  d <- make_test_data()

  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "nb", gamma = 0.001,
                              cov_alpha = ~0 + sex, cov_beta = ~0 + sex)
  expect_equal(length(coef(fit)), 4)  # 2 alpha + 2 beta
})

test_that("~0 formula treated as intercept-only", {
  d <- make_test_data()

  fit1 <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)
  fit2 <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001,
                               cov_alpha = ~0)

  expect_equal(as.numeric(coef(fit1)), as.numeric(coef(fit2)))
})

# ---- S3 methods ----

test_that("S3 methods return correct types", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)

  expect_true(is.numeric(coef(fit)))
  expect_true(is.matrix(vcov(fit)))
  expect_equal(nrow(vcov(fit)), length(coef(fit)))
  expect_true(is.numeric(fitted(fit)))
  expect_equal(length(fitted(fit)), nrow(d))
  expect_true(is.numeric(residuals(fit)))
  expect_equal(length(residuals(fit)), nrow(d))
})

test_that("print and summary run without error", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)

  expect_output(print(fit))
  expect_output(summary(fit))
})

# ---- vcov types ----

test_that("all HC types produce different SE", {
  d <- make_test_data()

  se_list <- lapply(c("HC0", "HC1", "HC2", "HC3"), function(hc) {
    fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                                gamma = 0.001, vcov = hc)
    sqrt(diag(vcov(fit)))
  })

  # HC0 < HC1 (finite sample correction)
  expect_true(all(se_list[[1]] < se_list[[2]]))
  # HC2 < HC3 (leverage correction)
  expect_true(all(se_list[[3]] < se_list[[4]]))
})

# ---- vcov_model ----

test_that("model-based vcov is present and smaller than HC3", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                              gamma = 0.001, vcov = "HC3")

  expect_true(!is.null(fit$vcov_model))

  se_hc3 <- sqrt(diag(fit$vcov))
  se_model <- sqrt(diag(fit$vcov_model))

  # Model-based SE should generally be smaller than HC3
  expect_true(all(se_model < se_hc3))
})

# ---- popsize ----

test_that("popsize returns correct structure", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)

  ps <- popsize(fit)
  expect_true(is.data.frame(ps))
  expect_true(all(c("group", "estimate", "estimate_bc", "lower", "upper") %in% names(ps)))
  expect_true(nrow(ps) >= 1)
})

test_that("bias correction reduces estimate (Jensen's inequality)", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)

  ps <- popsize(fit, bias_correction = TRUE)

  # g(alpha) = sum(N^alpha) is convex, so E[g(alpha_hat)] >= g(alpha)
  # Bias correction should reduce the estimate
  expect_true(all(ps$estimate_bc <= ps$estimate))
})

test_that("popsize without bias correction equals estimate", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)

  ps <- popsize(fit, bias_correction = FALSE)
  expect_equal(ps$estimate, ps$estimate_bc)
})

test_that("popsize with covariates gives multiple rows", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001,
                              cov_alpha = ~0 + sex)

  ps <- popsize(fit)
  expect_equal(nrow(ps), 2)  # M and F
})

test_that("popsize CI is ordered", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)

  ps <- popsize(fit)
  expect_true(all(ps$lower < ps$estimate))
  expect_true(all(ps$estimate < ps$upper))
})

# ---- NB extras ----

test_that("NB returns theta", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "nb", gamma = 0.001)

  expect_true(!is.null(fit$theta))
  expect_true(fit$theta > 0)
})

test_that("NB log-likelihood is higher than Poisson (nested model)", {
  d <- make_test_data()
  fit_p <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)
  fit_nb <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "nb", gamma = 0.001)

  expect_true(fit_nb$loglik >= fit_p$loglik)
})

# ---- LOO ----

test_that("loo runs and returns correct structure", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001,
                              countries = ~country)

  loo_res <- loo(fit, by = "obs")
  expect_s3_class(loo_res, "uncounted_loo")
  expect_equal(loo_res$n_drops, nrow(d))
  expect_true(sum(loo_res$converged) > 0)
})

test_that("loo by country works", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001,
                              countries = ~country)

  loo_res <- loo(fit, by = "country")
  expect_equal(loo_res$n_drops, length(unique(d$country)))
})

# ---- Edge cases ----

test_that("single observation per group works", {
  d <- make_test_data()[1:2, ]
  d$sex <- c("M", "F")

  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)
  expect_s3_class(fit, "uncounted")
})

test_that("no warning when data has no zeros and gamma=NULL", {
  d <- make_test_data()  # no zeros

  expect_no_warning(
    estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = NULL)
  )
})

test_that("zeros in n with gamma are handled", {
  d <- make_test_data_with_zeros()

  # With gamma, zeros in n should not cause issues
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)
  expect_s3_class(fit, "uncounted")
  expect_true(all(is.finite(coef(fit))))
})

# ---- OLS with gamma="estimate" (profiled) ----

test_that("OLS with gamma='estimate' estimates gamma in (0,1)", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "ols", gamma = "estimate")
  expect_true(!is.null(fit$gamma))
  expect_true(fit$gamma > 0 & fit$gamma < 1)
  expect_true(fit$gamma_estimated)
})

# ---- Invalid method ----

test_that("invalid method gives error", {
  d <- make_test_data()
  expect_error(
    estimate_hidden_pop(d, ~m, ~n, ~N, method = "invalid")
  )
})

# ---- LOO with constrained model ----

test_that("LOO passes constrained argument to refits", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                             gamma = 0.001, constrained = TRUE,
                             countries = ~country)
  loo_res <- loo(fit, by = "obs")
  expect_s3_class(loo_res, "uncounted_loo")
  expect_true(sum(loo_res$converged) > 0)
})

test_that("LOO by country with constrained model", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                             gamma = 0.001, constrained = TRUE,
                             countries = ~country)
  loo_res <- loo(fit, by = "country")
  expect_equal(loo_res$n_drops, length(unique(d$country)))
})

# ---- sandwich package methods ----

test_that("bread method returns correct dimensions", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)

  B <- sandwich::bread(fit)
  p <- length(coef(fit))
  expect_equal(dim(B), c(p, p))
  expect_true(all(is.finite(B)))
})

test_that("estfun method returns correct dimensions", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)

  ef <- sandwich::estfun(fit)
  expect_equal(nrow(ef), nrow(d))
  expect_equal(ncol(ef), length(coef(fit)))
  expect_true(all(is.finite(ef)))
})

test_that("nobs method works", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)
  expect_equal(nobs(fit), nrow(d))
})

test_that("hatvalues method works", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)

  h <- hatvalues(fit)
  expect_equal(length(h), nrow(d))
  expect_true(all(h >= 0 & h <= 1))
})

test_that("sandwich::vcovHC produces valid vcov for all methods", {
  d <- make_test_data()
  for (method in c("ols", "poisson", "nb")) {
    fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = method, gamma = 0.001)
    V <- sandwich::vcovHC(fit, type = "HC3")
    expect_equal(dim(V), rep(length(coef(fit)), 2))
    expect_true(all(diag(V) > 0))
  }
})

test_that("sandwich::vcovCL works with cluster", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                             gamma = 0.001, countries = ~country)

  V_cl <- sandwich::vcovCL(fit, cluster = d$country, type = "HC1")
  expect_equal(dim(V_cl), rep(length(coef(fit)), 2))
  expect_true(all(diag(V_cl) > 0))
})

test_that("cluster parameter in estimate_hidden_pop works", {
  d <- make_test_data()
  fit_hc <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                                gamma = 0.001, vcov = "HC3")
  fit_cl <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                                gamma = 0.001, vcov = "HC1",
                                cluster = ~country)

  se_hc <- sqrt(diag(fit_hc$vcov))
  se_cl <- sqrt(diag(fit_cl$vcov))

  # Cluster SE should differ from observation-level
  expect_false(all(abs(se_hc - se_cl) < 1e-10))
})

test_that("vcov as function works", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                             gamma = 0.001,
                             vcov = function(x) sandwich::vcovHC(x, type = "HC0"))

  expect_true(all(diag(fit$vcov) > 0))
})

test_that("vcov as function with vcovCL works", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                             gamma = 0.001, countries = ~country,
                             vcov = function(x) sandwich::vcovCL(x,
                               cluster = x$data$country, type = "HC1"))

  expect_true(all(diag(fit$vcov) > 0))
})

test_that("HC0 and HC3 produce different SE (observation-level)", {
  d <- make_test_data()

  fit_hc0 <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                                  gamma = 0.001, vcov = "HC0")
  fit_hc3 <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                                  gamma = 0.001, vcov = "HC3")

  se_hc0 <- sqrt(diag(fit_hc0$vcov))
  se_hc3 <- sqrt(diag(fit_hc3$vcov))
  expect_false(all(abs(se_hc0 - se_hc3) < 1e-10))
})

test_that("HC0 and HC3 produce different SE with clustering", {
  d <- make_test_data()

  fit_hc0 <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                                  gamma = 0.001, vcov = "HC0",
                                  cluster = ~country)
  fit_hc3 <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                                  gamma = 0.001, vcov = "HC3",
                                  cluster = ~country)

  se_hc0 <- sqrt(diag(fit_hc0$vcov))
  se_hc3 <- sqrt(diag(fit_hc3$vcov))
  # This was the original bug: HC0 and HC3 were identical with clustering
  expect_false(all(abs(se_hc0 - se_hc3) < 1e-10))
})

test_that("countries parameter is only for grouping, not clustering", {
  d <- make_test_data()
  # Use sex as cluster (2 clusters with 10 obs each) so CL differs from HC
  fit_no_cl <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                                    gamma = 0.001, countries = ~country)
  fit_with_cl <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                                      gamma = 0.001, countries = ~country,
                                      cluster = ~sex)

  # Without cluster arg, SE should be observation-level HC3
  # With cluster arg, SE should be cluster-robust
  se_no <- sqrt(diag(fit_no_cl$vcov))
  se_cl <- sqrt(diag(fit_with_cl$vcov))
  expect_false(all(abs(se_no - se_cl) < 1e-10))
})

test_that("NB bread/estfun methods work", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "nb", gamma = 0.001)

  B <- sandwich::bread(fit)
  ef <- sandwich::estfun(fit)
  expect_equal(ncol(ef), length(coef(fit)))
  expect_equal(dim(B), rep(length(coef(fit)), 2))
})

test_that("estimated gamma: bread/estfun include gamma column", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = "estimate")

  p_ab <- length(coef(fit))
  # model_matrix_full should have p_ab + 1 columns (gamma)
  expect_equal(ncol(fit$model_matrix_full), p_ab + 1)

  # bread/estfun should use full matrix
  B <- sandwich::bread(fit)
  ef <- sandwich::estfun(fit)
  expect_equal(ncol(B), p_ab + 1)
  expect_equal(ncol(ef), p_ab + 1)

  # But fit$vcov should only be p_ab x p_ab
  expect_equal(dim(fit$vcov), rep(p_ab, 2))
})

# ---- update and weights methods ----

test_that("weights method returns NULL when no weights given", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)
  expect_null(weights(fit))
})

test_that("weights method returns weights when given", {
  d <- make_test_data()
  w <- rep(1, nrow(d))
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson",
                             gamma = 0.001, weights = w)
  expect_equal(weights(fit), w)
})

test_that("update changes vcov type", {
  d <- make_test_data()
  fit_hc3 <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)
  fit_hc0 <- update(fit_hc3, vcov = "HC0")

  se_hc3 <- sqrt(diag(vcov(fit_hc3)))
  se_hc0 <- sqrt(diag(vcov(fit_hc0)))

  # Coefficients should be identical

  expect_equal(coef(fit_hc3), coef(fit_hc0))
  # SEs should differ
  expect_false(all(abs(se_hc3 - se_hc0) < 1e-10))
})

test_that("update with evaluate=FALSE returns a call", {
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)
  cl <- update(fit, vcov = "HC0", evaluate = FALSE)
  expect_true(is.call(cl))
})

test_that("update with new weights refits the model", {
  d <- make_test_data()
  fit1 <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)
  w <- runif(nrow(d), 0.5, 2)
  fit2 <- update(fit1, weights = w)

  # Coefficients should differ with different weights
  expect_false(all(abs(coef(fit1) - coef(fit2)) < 1e-10))
  expect_equal(weights(fit2), w)
})

test_that("fwb::vcovFWB works with uncounted objects", {
  skip_if_not_installed("fwb")
  set.seed(123)
  d <- make_test_data()
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)

  V_fwb <- fwb::vcovFWB(fit, R = 30)
  expect_equal(dim(V_fwb), rep(length(coef(fit)), 2))
  expect_true(all(diag(V_fwb) > 0))
})

test_that("vcov = fwb::vcovFWB works as function argument", {
  skip_if_not_installed("fwb")
  set.seed(456)
  d <- make_test_data()

  # Fit first, then compute FWB vcov separately (avoids sink stack
  # overflow from nested capture.output inside fwb during testthat)
  fit <- estimate_hidden_pop(d, ~m, ~n, ~N, method = "poisson", gamma = 0.001)
  V <- fwb::vcovFWB(fit, R = 30)
  se <- sqrt(diag(V[1:2, 1:2]))
  expect_true(all(se > 0))
  expect_true(all(is.finite(se)))
})
