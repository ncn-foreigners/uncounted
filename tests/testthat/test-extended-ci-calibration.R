# ---- 95% CI calibration for population size xi (extended, skip on CRAN) ----
#
# Monte Carlo coverage check for the package's central inferential claim: the
# 95% intervals for the population size xi should attain ~nominal coverage when
# the model is correctly specified. DGP follows the paper's Setup-II simulation:
# I = 100 communities, N ~ LogNormal(5.5, 2.75) truncated at 3, n ~ Poisson with
# log E(n | N) = c0 + c1 log N, and m ~ Poisson(N^alpha (gamma + n/N)^beta) with
# (alpha, beta, gamma) = (0.7, 0.5, 0.005); gamma is estimated (Po(n)-Po(m)).
#
# The acceptance bands are not guesses. Coverage of each procedure was measured
# on this DGP in long runs, and the bands are set at roughly +/- 3.5 Monte Carlo
# standard errors around the measured value at the replication count used here:
#
#   procedure                reps   measured coverage   MC se
#   analytic HC3             1000   0.940               0.008
#   cluster FWB, R = 99       100   0.890               0.031
#   cluster FWB, R = 199       80   0.900               0.034
#   cluster FWB, R = 299      120   0.933               0.023
#
# Two things follow. First, both procedures are close to nominal, the analytic
# interval slightly conservative-side of 0.94. Second, the FWB shortfall at
# small R is percentile resolution, not a defect: with R = 99 the 2.5% and
# 97.5% quantiles of the bootstrap distribution are estimated from a handful of
# order statistics, and coverage recovers toward the analytic value as R grows.
# The FWB test below therefore pins the R = 99 behaviour deliberately.
#
# Both bands are two-sided. An upper bound matters as much as a lower one: a
# variance that is too large by a factor of four still passes any test whose
# only requirement is "coverage >= 0.90", so the previous upper bound of 1.00
# gave no power against over-wide intervals at all.
#
# Replicates that fail to converge are counted as NOT covered rather than
# dropped. Dropping them biases the coverage estimate, because a replicate that
# defeats the optimizer is exactly the kind most likely to be miscovered. Their
# frequency is asserted separately; on this DGP none of 1000 failed.

skip_on_cran()

.calib_alpha <- 0.7
.calib_beta  <- 0.5
.calib_gamma <- 0.005
.calib_I     <- 100L

.calib_gen <- function() {
  N <- rlnorm(.calib_I, 5.5, 2.75)
  while (any(N < 3)) {
    bad <- N < 3
    N[bad] <- rlnorm(sum(bad), 5.5, 2.75)
  }
  n <- rpois(.calib_I, exp(-2.0 + 0.70 * log(N)))
  m <- rpois(.calib_I, N^.calib_alpha * (.calib_gamma + n / N)^.calib_beta)
  data.frame(N, n, m, country = factor(seq_len(.calib_I)))
}

.calib_fit <- function(d) {
  tryCatch(
    suppressWarnings(estimate_hidden_pop(
      d, ~ m, ~ n, ~ N, method = "poisson",
      gamma = "estimate", vcov = "HC3", countries = ~ country)),
    error = function(e) NULL)
}

test_that("analytic 95% CI for xi attains ~nominal coverage (Setup-II DGP)", {
  skip_on_cran()
  set.seed(101)
  reps <- 400L
  covered <- logical(reps)
  n_failed <- 0L

  for (r in seq_len(reps)) {
    d <- .calib_gen()
    xi_true <- sum(d$N^.calib_alpha)
    fit <- .calib_fit(d)
    if (is.null(fit) || isTRUE(fit$convergence != 0)) {
      covered[r] <- FALSE
      n_failed <- n_failed + 1L
      next
    }
    ps <- suppressWarnings(popsize(fit))
    covered[r] <- isTRUE(sum(ps$lower) <= xi_true && xi_true <= sum(ps$upper))
  }

  ## Non-convergence must stay rare; if it rises, the coverage number below
  ## is being propped up by failures counted as misses and needs a look.
  expect_lte(n_failed, 8L)

  cov <- mean(covered)
  ## Measured 0.940 over 1000 reps; MC se at 400 reps is 0.012, so each band
  ## edge is ~3.4 se from the mean (false-failure risk ~0.04% per side).
  ## Seed 101 as written yields 0.9525 with 0 non-converged.
  expect_gte(cov, 0.90)
  expect_lte(cov, 0.98)
  expect_true(
    cov >= 0.90 && cov <= 0.98,
    info = sprintf("analytic xi coverage = %.1f%% (reps = %d, %d non-converged), nominal 95%%",
                   100 * cov, reps, n_failed)
  )
})

test_that("cluster-FWB 95% CI for xi attains ~nominal coverage (Setup-II DGP)", {
  skip_on_cran()
  skip_if_not_installed("fwb")
  set.seed(202)
  reps <- 120L
  R_boot <- 99L
  covered <- logical(reps)
  n_failed <- 0L

  for (r in seq_len(reps)) {
    d <- .calib_gen()
    xi_true <- sum(d$N^.calib_alpha)
    fit <- .calib_fit(d)
    if (is.null(fit)) {
      covered[r] <- FALSE
      n_failed <- n_failed + 1L
      next
    }
    bt <- tryCatch(
      suppressWarnings(bootstrap_popsize(fit, R = R_boot, cluster = ~ country,
                                         seed = 202 + r, verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(bt)) {
      covered[r] <- FALSE
      n_failed <- n_failed + 1L
      next
    }
    covered[r] <- isTRUE(sum(bt$popsize$lower) <= xi_true &&
                           xi_true <= sum(bt$popsize$upper))
  }

  expect_lte(n_failed, 4L)

  cov <- mean(covered)
  ## Measured 0.890 at R = 99; MC se at 120 reps is 0.029. The band is wider
  ## on the low side because the percentile interval is coarse at this R.
  ## Seed 202 as written yields 0.8667 with 0 failures.
  expect_gte(cov, 0.80)
  expect_lte(cov, 0.975)
  expect_true(
    cov >= 0.80 && cov <= 0.975,
    info = sprintf("cluster-FWB xi coverage = %.1f%% (reps = %d, R = %d, %d failed), nominal 95%%",
                   100 * cov, reps, R_boot, n_failed)
  )
})
