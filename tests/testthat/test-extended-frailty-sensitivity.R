# ---- Frailty sensitivity tests (extended) ----

test_that("frailty_sensitivity auto degrees of freedom use clusters when available", {
  d <- positive_data()
  keep_country <- unique(d$country)[1:12]
  d <- droplevels(d[d$country %in% keep_country, ])

  fit <- estimate_hidden_pop(
    data = d,
    observed = ~ m,
    auxiliary = ~ n,
    reference_pop = ~ N,
    method = "poisson",
    gamma = 0.005,
    cov_alpha = ~ sex,
    cluster = ~ country,
    vcov = "HC1"
  )

  fs <- frailty_sensitivity(fit, by = ~ year, df_method = "auto", plot = FALSE)
  x_work <- uncounted:::.frailty_working_model(fit)$x_work
  expected_df <- length(unique(fit$cluster_var)) - qr(x_work)$rank

  expect_true(expected_df > 0)
  expect_true(all(fs$working$df == expected_df))
  expect_identical(fs$settings$resolved_df_method, "cluster")
})

test_that("frailty target derivative matches a numerical perturbation", {
  d <- positive_data()
  keep_country <- unique(d$country)[1:10]
  d <- droplevels(d[d$country %in% keep_country, ])
  d$ukr <- as.integer(d$country %in% keep_country[1:2])

  fit <- quick_fit(
    d,
    method = "poisson",
    gamma = 0.005,
    cov_alpha = ~ year * ukr + sex,
    cov_beta = ~ year
  )

  group_idx <- uncounted:::.popsize_group_index(fit, by = ~ year)
  v_alpha <- fit$vcov[seq_len(fit$p_alpha), seq_len(fit$p_alpha), drop = FALSE]
  target <- uncounted:::.frailty_target_stats(
    fit,
    idx = group_idx[[1]],
    label = names(group_idx)[1],
    v_alpha = v_alpha,
    df = 10
  )

  eps <- 1e-6
  alpha_eps <- fit$alpha_coefs + eps * target$contrast
  alpha_new <- as.numeric(fit$X_alpha[group_idx[[1]], , drop = FALSE] %*% alpha_eps)
  xi_new <- sum(fit$N[group_idx[[1]]]^alpha_new)
  numeric_deriv <- (xi_new - target$summary$xi_hat) / eps

  expect_equal(numeric_deriv, target$summary$d_xi_dtheta, tolerance = 1e-4)
})

test_that("frailty_sensitivity threshold handling reports already-crossed targets", {
  d <- positive_data()
  keep_country <- unique(d$country)[1:8]
  d <- droplevels(d[d$country %in% keep_country, ])
  d$ukr <- as.integer(d$country %in% keep_country[1:2])

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
    threshold = 1e12,
    direction = "decrease",
    plot = FALSE
  )

  expect_true(all(fs$robustness$status == "already_crossed"))
  expect_true(all(fs$robustness$rv_equal == 0))
  expect_true(all(fs$robustness$rv_extreme == 0))
})
