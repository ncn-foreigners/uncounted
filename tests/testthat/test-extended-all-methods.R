# ---- All methods on full data (extended, skip on CRAN) ----

test_that("all methods with estimated gamma on full data", {
  skip_on_cran()
  for (method in c("ols", "poisson", "nb")) {
    fit <- quick_fit(testdata, method = method, gamma = "estimate")
    expect_s3_class(fit, "uncounted")
    expect_true(!is.null(fit$gamma))
    expect_true(fit$gamma > 0)
  }
})

test_that("all methods with covariates on full data", {
  skip_on_cran()
  for (method in c("ols", "poisson", "nb")) {
    fit <- quick_fit(testdata, method = method, gamma = 0.005,
                     cov_alpha = ~0 + sex, cov_beta = ~0 + sex)
    expect_equal(length(coef(fit)), 4)
    expect_true(all(is.finite(coef(fit))))
  }
})

test_that("DGP recovery: alpha estimate close to truth", {
  skip_on_cran()
  fit <- quick_fit(testdata, gamma = dgp$gamma)
  alpha_hat <- coef(fit)[["alpha"]]
  expect_true(abs(alpha_hat - dgp$alpha) < 0.2,
              info = paste("alpha_hat =", round(alpha_hat, 3),
                           "vs truth =", dgp$alpha))
})

test_that("DGP recovery: estimated gamma in reasonable range", {
  skip_on_cran()
  fit <- quick_fit(testdata, gamma = "estimate")
  expect_true(fit$gamma > 0.0001 & fit$gamma < 0.1,
              info = paste("gamma_hat =", fit$gamma,
                           "vs truth =", dgp$gamma))
})

test_that("Poisson on full data: base R comparison", {
  skip_on_cran()
  d <- positive_data()
  g <- 0.005

  glm_fit <- glm(m ~ -1 + log(N) + log(g + n / N), family = poisson(), data = d)
  our_fit <- quick_fit(d, gamma = g)

  expect_equal(as.numeric(coef(our_fit)), as.numeric(coef(glm_fit)),
               tolerance = 1e-5)

  V_glm <- sandwich::vcovHC(glm_fit, type = "HC3")
  V_our <- vcov(our_fit)
  expect_equal(as.numeric(V_our), as.numeric(V_glm), tolerance = 1e-4)
})

test_that("OLS on full positive data: base R comparison", {
  skip_on_cran()
  d <- positive_data()
  g <- 0.005

  lm_fit <- lm(log(m) ~ 0 + log(N) + log(g + n / N), data = d)
  our_fit <- quick_fit(d, method = "ols", gamma = g)

  expect_equal(as.numeric(coef(our_fit)), as.numeric(coef(lm_fit)),
               tolerance = 1e-5)
})

test_that("Poisson with covariates on full data: base R comparison", {
  skip_on_cran()
  d <- positive_data()

  glm_fit <- glm(m ~ -1 + sex:log(N) + log(n / N), family = poisson(), data = d)
  our_fit <- quick_fit(d, gamma = NULL, cov_alpha = ~0 + sex)

  expect_equal(our_fit |> coef() |> (\(x) x[["alpha:sexF"]])(),
               coef(glm_fit)[["sexF:log(N)"]], tolerance = 1e-5)
  expect_equal(our_fit |> coef() |> (\(x) x[["alpha:sexM"]])(),
               coef(glm_fit)[["sexM:log(N)"]], tolerance = 1e-5)
})

test_that("popsize on full data: all methods give positive estimates", {
  skip_on_cran()
  for (method in c("poisson", "nb")) {
    fit <- quick_fit(testdata, method = method, gamma = 0.005)
    ps <- popsize(fit, bias_correction = TRUE)
    expect_true(all(ps$estimate > 0), info = paste("Method:", method))
    expect_true(all(ps$estimate_bc > 0), info = paste("Method:", method))
    expect_true(all(ps$lower < ps$estimate), info = paste("Method:", method))
    expect_true(all(ps$estimate < ps$upper), info = paste("Method:", method))
  }
})

# ---- iOLS simulation study ----

test_that("iOLS DGP recovery: alpha close to truth", {
  skip_on_cran()
  fit <- quick_fit(testdata, method = "iols", gamma = dgp$gamma)
  alpha_hat <- coef(fit)[["alpha"]]
  expect_true(abs(alpha_hat - dgp$alpha) < 0.2,
              info = paste("iOLS alpha_hat =", round(alpha_hat, 3),
                           "vs truth =", dgp$alpha))
})

test_that("iOLS DGP recovery: beta close to truth", {
  skip_on_cran()
  fit <- quick_fit(testdata, method = "iols", gamma = dgp$gamma)
  beta_hat <- coef(fit)[["beta"]]
  expect_true(abs(beta_hat - dgp$beta) < 0.3,
              info = paste("iOLS beta_hat =", round(beta_hat, 3),
                           "vs truth =", dgp$beta))
})

test_that("iOLS converges and satisfies GPML score on full data", {
  skip_on_cran()
  fit <- quick_fit(testdata, method = "iols", gamma = dgp$gamma)
  expect_equal(fit$convergence, 0L)
  mu <- fit$fitted.values
  score <- crossprod(fit$model_matrix_full, fit$m / pmax(mu, 1e-10) - 1) / nrow(testdata)
  expect_true(max(abs(score)) < 1e-4)
})

test_that("iOLS vs Poisson: coefficients in same direction on full data", {
  skip_on_cran()
  fit_po <- quick_fit(testdata, gamma = dgp$gamma)
  fit_io <- quick_fit(testdata, method = "iols", gamma = dgp$gamma)
  # Same sign for both coefficients
  expect_equal(sign(coef(fit_po)), sign(coef(fit_io)))
  # Within 50% of each other
  rel_diff <- abs(coef(fit_io) - coef(fit_po)) / (abs(coef(fit_po)) + 0.01)
  expect_true(all(rel_diff < 0.5),
              info = paste("Poisson:", paste(round(coef(fit_po), 4), collapse = ", "),
                           "| iOLS:", paste(round(coef(fit_io), 4), collapse = ", ")))
})

test_that("iOLS with covariates on full data", {
  skip_on_cran()
  fit <- quick_fit(testdata, method = "iols", gamma = 0.005,
                   cov_alpha = ~0 + sex, cov_beta = ~0 + sex)
  expect_equal(length(coef(fit)), 4)
  expect_true(all(is.finite(coef(fit))))
  expect_equal(fit$convergence, 0L)
})

test_that("iOLS popsize on full data: positive estimates", {
  skip_on_cran()
  fit <- quick_fit(testdata, method = "iols", gamma = 0.005)
  ps <- popsize(fit, bias_correction = TRUE)
  expect_true(all(ps$estimate > 0))
  expect_true(all(ps$estimate_bc > 0))
})

test_that("All methods comparison: DGP recovery summary", {
  skip_on_cran()
  # Fit all methods with true gamma
  methods <- c("ols", "poisson", "nb", "iols")
  results <- lapply(methods, function(m) {
    suppressWarnings(quick_fit(testdata, method = m, gamma = dgp$gamma))
  })
  names(results) <- methods

  # Compare alpha estimates
  alphas <- sapply(results, function(f) coef(f)[["alpha"]])
  cat("\n  DGP alpha =", dgp$alpha, "\n")
  for (m in methods) cat("  ", m, "alpha =", round(alphas[m], 4), "\n")

  # All alpha estimates should be within 0.3 of truth
  expect_true(all(abs(alphas - dgp$alpha) < 0.3),
              info = paste("Alphas:", paste(round(alphas, 4), collapse = ", ")))

  # Compare beta estimates
  betas <- sapply(results, function(f) coef(f)[["beta"]])
  cat("  DGP beta =", dgp$beta, "\n")
  for (m in methods) cat("  ", m, "beta =", round(betas[m], 4), "\n")

  expect_true(all(abs(betas - dgp$beta) < 0.3),
              info = paste("Betas:", paste(round(betas, 4), collapse = ", ")))

  # Compare popsize
  popsizes <- sapply(results, function(f) sum(popsize(f)$estimate))
  cat("  Population sizes:\n")
  for (m in methods) cat("  ", m, "xi =", round(popsizes[m]), "\n")
})

# ---- Monte Carlo simulation: bias correction ----

test_that("Simulation: multiplicative BC reduces bias for xi", {
  skip_on_cran()
  set.seed(2025)
  R <- 200  # enough for bias estimation, not too slow
  N_vec <- testdata$N
  n_vec <- testdata$n
  ratio_vec <- n_vec / N_vec
  mu_true <- N_vec^dgp$alpha * (dgp$gamma + ratio_vec)^dgp$beta
  xi_true <- sum(N_vec^dgp$alpha)

  xi_plugin <- xi_bc <- numeric(R)
  n_ok <- 0
  for (r in seq_len(R)) {
    m_sim <- rpois(length(mu_true), mu_true)
    d_sim <- data.frame(m = m_sim, n = n_vec, N = N_vec)
    fit <- tryCatch(
      estimate_hidden_pop(d_sim, ~m, ~n, ~N, method = "poisson", gamma = dgp$gamma),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    n_ok <- n_ok + 1
    xi_plugin[n_ok] <- sum(popsize(fit, bias_correction = FALSE)$estimate)
    xi_bc[n_ok] <- sum(popsize(fit, bias_correction = TRUE)$estimate_bc)
  }

  xi_plugin <- xi_plugin[seq_len(n_ok)]
  xi_bc <- xi_bc[seq_len(n_ok)]

  cat("\n  xi_true =", round(xi_true), "\n")
  cat("  mean(xi_plugin) =", round(mean(xi_plugin)),
      " bias =", round(mean(xi_plugin) - xi_true), "\n")
  cat("  mean(xi_bc) =", round(mean(xi_bc)),
      " bias =", round(mean(xi_bc) - xi_true), "\n")
  cat("  n_ok =", n_ok, "/", R, "\n")

  # Plugin should be upward biased (Jensen's inequality)
  expect_true(mean(xi_plugin) > xi_true,
              info = "Plugin should overestimate xi (Jensen)")

  # Multiplicative BC should reduce bias
  expect_true(abs(mean(xi_bc) - xi_true) < abs(mean(xi_plugin) - xi_true),
              info = "BC should reduce bias relative to plugin")

  # BC should always be positive
  expect_true(all(xi_bc > 0), info = "BC must always be positive")
})
