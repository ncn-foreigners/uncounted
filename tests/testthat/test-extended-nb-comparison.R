# ---- NB vs MASS::glm.nb comparison (extended, skip on CRAN) ----

test_that("NB coefficients close to MASS::glm.nb", {
  skip_on_cran()
  d <- positive_data()
  g <- 0.005

  mass_fit <- MASS::glm.nb(m ~ -1 + log(N) + log(g + n / N), data = d)
  our_fit  <- quick_fit(d, method = "nb", gamma = g)

  mass_coef <- as.numeric(coef(mass_fit))
  our_coef  <- as.numeric(coef(our_fit))

  # Coefficients within 10% relative difference
  rel_diff <- abs(our_coef - mass_coef) / (abs(mass_coef) + 1e-6)
  expect_true(all(rel_diff < 0.10),
              info = paste("Relative diff:", paste(round(rel_diff, 4), collapse = ", ")))
})

test_that("NB theta close to MASS theta", {
  skip_on_cran()
  d <- positive_data()
  g <- 0.005

  mass_fit <- MASS::glm.nb(m ~ -1 + log(N) + log(g + n / N), data = d)
  our_fit  <- quick_fit(d, method = "nb", gamma = g)

  # Theta within 50% (different optimization paths)
  rel_diff <- abs(our_fit$theta - mass_fit$theta) / mass_fit$theta
  expect_true(rel_diff < 0.50,
              info = paste("our theta =", round(our_fit$theta, 3),
                           "MASS theta =", round(mass_fit$theta, 3)))
})

test_that("NB log-likelihood close to MASS", {
  skip_on_cran()
  d <- positive_data()
  g <- 0.005

  mass_fit <- MASS::glm.nb(m ~ -1 + log(N) + log(g + n / N), data = d)
  our_fit  <- quick_fit(d, method = "nb", gamma = g)

  # Log-likelihoods within 2 units
  expect_true(abs(our_fit$loglik - as.numeric(logLik(mass_fit))) < 2,
              info = paste("our ll =", round(our_fit$loglik, 2),
                           "MASS ll =", round(as.numeric(logLik(mass_fit)), 2)))
})

test_that("NB with covariates: coefficients close to MASS", {
  skip_on_cran()
  d <- positive_data()
  g <- 0.005

  mass_fit <- MASS::glm.nb(m ~ -1 + sex:log(N) + log(g + n / N), data = d)
  our_fit  <- quick_fit(d, method = "nb", gamma = g, cov_alpha = ~0 + sex)

  # Check beta coefficient (shared between models)
  our_beta  <- coef(our_fit)[["beta"]]
  mass_beta <- coef(mass_fit)[["log(g + n/N)"]]

  rel_diff <- abs(our_beta - mass_beta) / (abs(mass_beta) + 1e-6)
  expect_true(rel_diff < 0.10,
              info = paste("our beta =", round(our_beta, 4),
                           "MASS beta =", round(mass_beta, 4)))
})
