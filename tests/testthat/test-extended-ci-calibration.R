# ---- 95% CI calibration for population size xi (extended, skip on CRAN) ----
#
# Monte Carlo coverage check for the package's central inferential claim: the
# 95% intervals for the population size xi should attain ~nominal coverage when
# the model is correctly specified. DGP follows the paper's Setup-II simulation:
# I = 100 communities, N ~ LogNormal(5.5, 2.75) truncated at 3, n ~ Poisson with
# log E(n | N) = c0 + c1 log N, and m ~ Poisson(N^alpha (gamma + n/N)^beta) with
# (alpha, beta, gamma) = (0.7, 0.5, 0.005); gamma is estimated (Po(n)-Po(m)).
#
# Seeds are fixed for reproducibility, but the assertions are bands around the
# nominal 0.95 -- they hold because the estimator is calibrated, not because of
# a specific stream, so they remain valid across platforms/RNG versions. The
# lower bounds are the meaningful guard: they catch gross under-coverage such as
# the ~25-40% seen when the mean structure is misspecified.

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
  reps <- 200L
  covered <- logical(0)
  for (r in seq_len(reps)) {
    d <- .calib_gen()
    xi_true <- sum(d$N^.calib_alpha)
    fit <- .calib_fit(d)
    if (is.null(fit) || isTRUE(fit$convergence != 0)) next
    ps <- suppressWarnings(popsize(fit))
    covered <- c(covered, isTRUE(sum(ps$lower) <= xi_true && xi_true <= sum(ps$upper)))
  }
  expect_gte(length(covered), 190L)  # nearly all reps converge
  cov <- mean(covered)
  # nominal 0.95; seed 101 yields ~0.965. Lower bound flags gross under-coverage.
  expect_true(cov >= 0.90 && cov <= 1.00,
              info = sprintf("analytic xi coverage = %.1f%% (n = %d), nominal 95%%",
                             100 * cov, length(covered)))
})

test_that("cluster-FWB 95% CI for xi attains ~nominal coverage (Setup-II DGP)", {
  skip_on_cran()
  set.seed(202)
  reps <- 30L
  R_boot <- 49L
  covered <- logical(0)
  for (r in seq_len(reps)) {
    d <- .calib_gen()
    xi_true <- sum(d$N^.calib_alpha)
    fit <- .calib_fit(d)
    if (is.null(fit) || isTRUE(fit$convergence != 0)) next
    bt <- tryCatch(
      suppressWarnings(bootstrap_popsize(fit, R = R_boot, cluster = ~ country,
                                         seed = 202 + r, verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(bt)) next
    covered <- c(covered, isTRUE(sum(bt$popsize$lower) <= xi_true &&
                                   xi_true <= sum(bt$popsize$upper)))
  }
  expect_gte(length(covered), 28L)
  cov <- mean(covered)
  # nominal 0.95; seed 202 yields ~0.967 over 30 reps. Wider floor (fewer reps).
  expect_true(cov >= 0.80 && cov <= 1.00,
              info = sprintf("cluster-FWB xi coverage = %.1f%% (n = %d), nominal 95%%",
                             100 * cov, length(covered)))
})
