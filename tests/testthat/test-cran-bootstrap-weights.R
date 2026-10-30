# ---- bootstrap_popsize(): refits use the original model (CRAN) ----

fit_w <- function(w) {
  estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005, weights = w)
}

boot_quiet <- function(fit, ...) {
  bootstrap_popsize(fit, seed = 1, verbose = FALSE, ...)
}

test_that("FWB weights multiply the user's weights", {
  w1 <- ifelse(testdata$N > median(testdata$N), 1, 0.01)
  w2 <- rep(c(1, 3, 0.5, 2), length.out = nrow(testdata))
  b1 <- boot_quiet(fit_w(w1), R = 20)
  b2 <- boot_quiet(fit_w(w2), R = 20)
  ## Used to be identical: the replicates ignored the weights
  expect_false(isTRUE(all.equal(b1$t, b2$t)))

  ## Replicates agree with a manual refit under weights = user * FWB
  set.seed(1)
  fwb_out <- fwb::fwb(testdata, function(data, w) {
    fit <- estimate_hidden_pop(data, ~m, ~n, ~N, gamma = 0.005, weights = w1 * w,
                               vcov = "HC0")
    sum(popsize(fit, bias_correction = FALSE)$estimate)
  }, R = 20, verbose = FALSE)
  expect_equal(as.numeric(b1$t), as.numeric(fwb_out$t), tolerance = 1e-5)
})

test_that("constant weights give the unweighted bootstrap", {
  un <- boot_quiet(estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005), R = 20)
  w2 <- boot_quiet(fit_w(rep(2, nrow(testdata))), R = 20)
  expect_equal(un$t, w2$t, tolerance = 1e-5)
})

test_that("the bootstrap keeps the rest of the model specification", {
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", cov_alpha = ~sex,
                             gamma = "estimate", gamma_bounds = c(1e-6, 0.1),
                             constrained = TRUE, link_rho = "cloglog")
  b <- boot_quiet(fit, R = 5)
  expect_gt(b$n_converged, 0L)
  gam <- b$boot_params[, "gamma"]
  expect_true(all(gam[is.finite(gam)] <= 0.1))
})

test_that("cluster bootstrap with one cluster gives a clear error", {
  d <- single_cluster_data()
  fit <- quick_fit(d)
  expect_error(bootstrap_popsize(fit, R = 5, cluster = ~one_cluster, verbose = FALSE),
               "at least 2 clusters")
})

test_that("cheap CI matches fwb and works with very few replications", {
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005)
  b <- boot_quiet(fit, R = 5, ci_type = "cheap")
  ## Same formula as fwb.ci(type = "cheap") with h = log
  xi0 <- b$t0
  s <- sqrt(mean((log(b$t[, 1]) - log(xi0))^2))
  expect_equal(b$popsize$lower, exp(log(xi0) - qt(0.975, 5) * s))
  expect_equal(b$popsize$upper, exp(log(xi0) + qt(0.975, 5) * s))

  set.seed(1)
  ## fwb.ci() allows only h = identity, so compare on the log statistic
  f <- fwb::fwb(testdata, function(data, w) log(stats::weighted.mean(data$m, w)),
                R = 5, verbose = FALSE)
  ours <- .bootstrap_ci(exp(f$t[, 1]), exp(f$t0), 0.95, "cheap")
  theirs <- fwb::fwb.ci(f, type = "cheap")$cheap
  expect_equal(ours, exp(as.numeric(theirs[, 2:3])), tolerance = 1e-8)

  ## Point estimate defaults to the plug-in; R = 1 is allowed
  expect_equal(b$popsize$estimate, xi0)
  b1 <- boot_quiet(fit, R = 1, ci_type = "cheap")
  expect_true(is.finite(b1$popsize$lower) && b1$popsize$lower < b1$popsize$upper)
})

# ---- boot_params lines up with the replicates ----

test_that("boot_params rows are the replicates, not the full-sample fit", {
  ## fwb calls the statistic once on the full sample before the replicates, so
  ## a counter-based store had the full-sample fit in row 1 and lost the last
  ## replicate.
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005)
  b <- boot_quiet(fit, R = 10)
  expect_equal(nrow(b$boot_params), 10L)
  xi_from_params <- vapply(b$boot_params[, "alpha"],
                           function(a) sum(testdata$N^a), numeric(1))
  expect_equal(unname(xi_from_params), unname(b$t[, 1]))
  expect_false(isTRUE(all.equal(unname(b$boot_params[1, ]), unname(coef(fit)))))
})

test_that("boot_params names cover gamma and theta", {
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb")
  b <- boot_quiet(fit, R = 3)
  expect_identical(colnames(b$boot_params), c("alpha", "beta", "gamma", "theta"))
  expect_true(all(is.finite(b$boot_params)))
  expect_identical(colnames(b$t), b$popsize$group)
})
