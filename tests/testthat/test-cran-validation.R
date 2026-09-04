## Input validation for estimate_hidden_pop().
##
## Each of these inputs previously either failed deep inside optim()/sandwich
## with a message that named none of the user's arguments, or was accepted
## silently. The assertions below pin the diagnostics to the offending
## argument.

base_fit <- function(...) {
  estimate_hidden_pop(
    data = testdata, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
    ...
  )
}

test_that("gamma must be 'estimate', a single non-negative number, or NULL", {
  expect_error(base_fit(method = "poisson", gamma = "foo"), "must be")
  expect_error(base_fit(method = "poisson", gamma = -0.01),
               "single finite non-negative")
  expect_error(base_fit(method = "poisson", gamma = c(0.1, 0.2)),
               "single finite non-negative")
  expect_error(base_fit(method = "poisson", gamma = NA_real_),
               "single finite non-negative")
  expect_error(base_fit(method = "poisson", gamma = TRUE), "must be")

  ## The three documented forms still work.
  expect_s3_class(base_fit(method = "poisson", gamma = 0.005), "uncounted")
  expect_s3_class(base_fit(method = "poisson", gamma = "estimate"), "uncounted")
})

test_that("gamma_bounds must be a finite increasing pair", {
  expect_error(
    base_fit(method = "poisson", gamma = "estimate",
             gamma_bounds = c(0.5, 1e-10)),
    "0 <= lower < upper"
  )
  expect_error(
    base_fit(method = "poisson", gamma = "estimate", gamma_bounds = 0.5),
    "length 2"
  )
  expect_error(
    base_fit(method = "poisson", gamma = "estimate",
             gamma_bounds = c(1e-10, Inf)),
    "finite"
  )
  ## Bounds are only meaningful when gamma is estimated, so a fixed gamma
  ## must not be rejected because of them.
  expect_s3_class(
    base_fit(method = "poisson", gamma = 0.005, gamma_bounds = c(0.5, 1e-10)),
    "uncounted"
  )
})

test_that("weights are checked against the data rather than recycled", {
  n_obs <- nrow(testdata)
  ## Previously recycled silently, so 5 weights fitted 200 observations.
  expect_error(base_fit(method = "poisson", gamma = 0.005, weights = 1:5),
               "one value per observation")
  expect_error(
    base_fit(method = "poisson", gamma = 0.005, weights = rep(-1, n_obs)),
    "non-negative"
  )
  expect_error(
    base_fit(method = "poisson", gamma = 0.005,
             weights = c(NA, rep(1, n_obs - 1))),
    "missing or non-finite"
  )
  expect_error(
    base_fit(method = "poisson", gamma = 0.005, weights = rep(0, n_obs)),
    "all zero"
  )
  expect_error(
    base_fit(method = "poisson", gamma = 0.005,
             weights = rep("a", n_obs)),
    "must be numeric"
  )
})

test_that("observed, auxiliary and reference_pop are checked", {
  d_na <- testdata; d_na$m[1] <- NA
  expect_error(
    estimate_hidden_pop(d_na, ~ m, ~ n, ~ N, method = "poisson", gamma = 0.005),
    "missing value"
  )

  d_neg <- testdata; d_neg$m[1] <- -5
  expect_error(
    estimate_hidden_pop(d_neg, ~ m, ~ n, ~ N, method = "poisson", gamma = 0.005),
    "negative"
  )

  d_negn <- testdata; d_negn$n[1] <- -1
  expect_error(
    estimate_hidden_pop(d_negn, ~ m, ~ n, ~ N, method = "poisson", gamma = 0.005),
    "negative"
  )

  d_zero <- testdata; d_zero$N[1] <- 0
  expect_error(
    estimate_hidden_pop(d_zero, ~ m, ~ n, ~ N, method = "poisson", gamma = 0.005),
    "strictly positive"
  )

  ## N < 1 is fitted but flagged: log(N) < 0 flips the sign of that
  ## observation's contribution to xi = sum(N^alpha).
  d_small <- testdata; d_small$N[1] <- 0.5
  expect_warning(
    estimate_hidden_pop(d_small, ~ m, ~ n, ~ N, method = "poisson",
                        gamma = 0.005),
    "below 1"
  )
})

test_that("gamma = NULL is rejected when the auxiliary count has zeros", {
  ## log(n/N) is -Inf for those rows; this used to warn and then die with
  ## "NA/NaN/Inf in 'x'" from the optimizer.
  expect_true(any(testdata$n == 0))
  expect_error(base_fit(method = "poisson", gamma = NULL), "undefined")

  ## With all-positive auxiliary counts it is a valid specification.
  d_pos <- testdata[testdata$n > 0, ]
  expect_s3_class(
    estimate_hidden_pop(d_pos, ~ m, ~ n, ~ N, method = "poisson", gamma = NULL),
    "uncounted"
  )
})

test_that("theta_start must be a single positive value for NB", {
  expect_error(base_fit(method = "nb", gamma = 0.005, theta_start = -1),
               "positive")
  expect_error(base_fit(method = "nb", gamma = 0.005, theta_start = c(1, 2)),
               "positive")
  ## Irrelevant for non-NB methods, so it must not be validated there.
  expect_s3_class(
    base_fit(method = "poisson", gamma = 0.005, theta_start = -1),
    "uncounted"
  )
})
