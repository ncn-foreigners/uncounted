# ---- lrtest(): conditions of the likelihood ratio test (CRAN) ----

test_that("lrtest() requires the same data and weights", {
  ## Only n_obs used to be compared: weights of 2 gave LR = 1376 for a covariate
  ## that does nothing
  f0 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005)
  f1 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~sex, gamma = 0.005)
  fw <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005,
                            weights = rep(2, nrow(testdata)))
  expect_error(lrtest(fw, f1), "differ in: weights")

  d <- testdata
  d$m <- rev(d$m)
  fy <- suppressWarnings(estimate_hidden_pop(d, ~m, ~n, ~N, gamma = 0.005))
  expect_error(lrtest(fy, f1), "differ in: observed")

  ## Unit weights are the same as no weights; equal weights pass
  fu <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005,
                            weights = rep(1, nrow(testdata)))
  expect_no_error(lrtest(fu, f1))
  fw1 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~sex, gamma = 0.005,
                             weights = rep(2, nrow(testdata)))
  expect_no_error(lrtest(fw, fw1))
})

test_that("an NB model is never treated as nested in a Poisson model", {
  ## Parameter counts used to decide, so NB(intercept) vs Poisson(~country) ran as
  ## "NB nested in Poisson" with no warning
  nb0 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  po_c <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~country, gamma = 0.005)
  expect_error(lrtest(nb0, po_c), "not nested")
  expect_error(lrtest(po_c, nb0), "not nested")
  ## Same parameter count, Poisson with the extra covariate
  po_s <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~sex, gamma = 0.005)
  expect_error(lrtest(po_s, nb0), "not nested")
})

test_that("Poisson nested in NB is accepted in either argument order", {
  po <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005)
  nb <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  nb_s <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", cov_alpha = ~sex,
                              gamma = 0.005)
  r1 <- lrtest(po, nb)
  r2 <- lrtest(nb, po)
  expect_equal(r1$statistic, r2$statistic)
  expect_equal(r1$df, 1)
  expect_match(r2$model1, "POISSON")
  expect_equal(lrtest(nb_s, po)$df, 2)
})

test_that("the boundary p-value is the chi-squared mixture", {
  po <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005)
  nb <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  r <- lrtest(po, nb)
  expect_equal(r$p_value, 0.5 * pchisq(r$statistic, 1, lower.tail = FALSE))

  ## Extra covariates in the NB model: used to be a plain chi2 on df degrees
  nb_s <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", cov_alpha = ~sex,
                              gamma = 0.005)
  r2 <- lrtest(po, nb_s)
  expect_equal(r2$df, 2)
  expect_equal(r2$p_value, 0.5 * pchisq(r2$statistic, 1, lower.tail = FALSE) +
                 0.5 * pchisq(r2$statistic, 2, lower.tail = FALSE))
  expect_output(print(r2), "0.5 \\* chi2_1 \\+ 0.5 \\* chi2_2")
})

test_that("LR = 0 at the boundary gives p = 1", {
  ## The old formula, 0.5 * P(chi2_1 > 0), gave 0.5
  expect_equal(.lr_boundary_p(0, 1), 1)
  expect_equal(.lr_boundary_p(0, 2), 1)
  expect_equal(.lr_boundary_p(3.84, 1), 0.5 * pchisq(3.84, 1, lower.tail = FALSE))

  ## Underdispersed counts: the NB fit stays at the Poisson limit
  d <- testdata
  d$m <- round(d$N^0.7 * (0.005 + d$n / d$N)^0.55)
  po <- estimate_hidden_pop(d, ~m, ~n, ~N, gamma = 0.005)
  nb <- suppressWarnings(estimate_hidden_pop(d, ~m, ~n, ~N, method = "nb", gamma = 0.005))
  r <- lrtest(po, nb)
  expect_lt(r$statistic, 1e-3)
  expect_gt(r$p_value, 0.45)
})

test_that("a same-family comparison keeps the plain chi-squared", {
  f0 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005)
  f1 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~sex, gamma = 0.005)
  r <- lrtest(f0, f1)
  expect_false(r$boundary)
  expect_equal(r$p_value, pchisq(r$statistic, 1, lower.tail = FALSE))
  expect_false(any(grepl("Boundary", capture.output(print(r)))))
})

test_that("lrtest() reports overdispersion for same-family comparisons", {
  mu <- testdata$N^0.7 * (0.005 + testdata$n / testdata$N)^0.55
  set.seed(1)
  d <- testdata
  d$m <- rnbinom(nrow(d), mu = mu, size = 2)
  f0 <- estimate_hidden_pop(d, ~m, ~n, ~N, gamma = 0.005)
  f1 <- estimate_hidden_pop(d, ~m, ~n, ~N, cov_alpha = ~sex, gamma = 0.005)
  r <- lrtest(f0, f1)
  p <- .count_params(f1)
  expect_equal(r$dispersion,
               sum((d$m - fitted(f1))^2 / fitted(f1)) / (nrow(d) - p))
  expect_gt(r$dispersion, 1.5)
  expect_output(print(r), "Pearson dispersion")

  set.seed(1)
  d$m <- rpois(nrow(d), mu)
  f0 <- estimate_hidden_pop(d, ~m, ~n, ~N, gamma = 0.005)
  f1 <- estimate_hidden_pop(d, ~m, ~n, ~N, cov_alpha = ~sex, gamma = 0.005)
  r <- lrtest(f0, f1)
  expect_lt(r$dispersion, 1.5)
  expect_false(any(grepl("Pearson", capture.output(print(r)))))
})

test_that("Poisson against NB does not get the dispersion note", {
  po <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005)
  nb <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  r <- lrtest(po, nb)
  expect_true(is.na(r$dispersion))
  expect_false(any(grepl("Pearson", capture.output(print(r)))))
})

# ---- test = "Wald" ----

test_that("the Wald test of one added dummy is the squared z statistic", {
  f0 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005)
  f1 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~sex, gamma = 0.005)
  r <- lrtest(f0, f1, test = "Wald")
  j <- "alpha:sexM"
  expect_equal(r$test, "Wald")
  expect_equal(r$df, 1L)
  expect_equal(r$statistic, unname(coef(f1)[j]^2 / vcov(f1)[j, j]))
  expect_equal(r$p_value, pchisq(r$statistic, 1, lower.tail = FALSE))
  expect_output(print(r), "Wald statistic")
  expect_equal(lrtest(f1, f0, test = "Wald")$statistic, r$statistic)
})

test_that("the Wald statistic does not depend on the parameterisation", {
  f0 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005)
  f1 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~year, cov_beta = ~sex,
                            gamma = 0.005)
  r <- lrtest(f0, f1, test = "Wald")
  n_year <- nlevels(testdata$year)
  expect_equal(r$df, n_year)

  ## Direct: all contrasts are zero
  idx <- grep("year|sex", names(coef(f1)))
  b <- coef(f1)[idx]
  expect_equal(r$statistic, as.numeric(crossprod(b, solve(vcov(f1)[idx, idx], b))),
               tolerance = 1e-6)

  ## Nested but not by dropping columns: ~sex within ~sex + year, without intercepts
  a <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~0 + sex, gamma = 0.005)
  b0 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~sex, gamma = 0.005)
  c1 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~sex + year, gamma = 0.005)
  expect_equal(lrtest(a, c1, test = "Wald")$statistic,
               lrtest(b0, c1, test = "Wald")$statistic, tolerance = 1e-6)
})

test_that("test = 'Wald' refuses what it cannot test", {
  po <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005)
  nb <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  expect_error(lrtest(po, nb, test = "Wald"), "same method")

  sex <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~sex, gamma = 0.005)
  year <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~year, gamma = 0.005)
  expect_warning(expect_error(lrtest(sex, year, test = "Wald"), "nested"), "not be nested")

  est <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~sex)
  expect_error(lrtest(po, est, test = "Wald"), "same gamma")
  expect_error(lrtest(po, po, test = "Wald"), "nothing")
})

test_that("the Wald test works for NB and with an estimated gamma", {
  n0 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", gamma = 0.005)
  n1 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, method = "nb", cov_alpha = ~sex, gamma = 0.005)
  j <- "alpha:sexM"
  expect_equal(lrtest(n0, n1, test = "Wald")$statistic,
               unname(coef(n1)[j]^2 / vcov(n1)[j, j]))
  g0 <- estimate_hidden_pop(testdata, ~m, ~n, ~N)
  g1 <- estimate_hidden_pop(testdata, ~m, ~n, ~N, cov_alpha = ~sex)
  expect_equal(lrtest(g0, g1, test = "Wald")$statistic,
               unname(coef(g1)[j]^2 / vcov(g1)[j, j]))
})
