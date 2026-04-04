# ---- NB sandwich with theta: analytical verification ----

test_that("NB per-obs score sums to near-zero at MLE", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005)
  expect_true(!is.null(fit$score_full))
  score_sums <- colSums(fit$score_full)
  expect_true(all(abs(score_sums) < 0.1),
              info = paste("Score sums:", paste(round(score_sums, 4), collapse = ", ")))
})

test_that("NB Hessian is positive definite (valid MLE)", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005)
  expect_true(!is.null(fit$hessian_nll))
  eigenvalues <- eigen(fit$hessian_nll, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0),
              info = paste("Eigenvalues:", paste(round(eigenvalues, 3), collapse = ", ")))
})

test_that("NB meat matrix is positive semi-definite", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005)
  meat <- crossprod(fit$score_full)
  eigenvalues <- eigen(meat, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues >= -1e-10))
})

test_that("NB full vcov is positive definite", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005)
  V_full <- fit$vcov_full
  expect_true(nrow(V_full) > length(coef(fit)))
  eigenvalues <- eigen(V_full, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0))
})

test_that("NB coefficient vcov equals submatrix of full vcov", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005)
  p_ab <- length(coef(fit))
  V_sub <- fit$vcov_full[seq_len(p_ab), seq_len(p_ab)]
  expect_equal(as.numeric(fit$vcov), as.numeric(V_sub), tolerance = 1e-10)
})

test_that("NB theta_se is finite and positive", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005)
  expect_true(!is.null(fit$theta_se))
  expect_true(is.finite(fit$theta_se))
  expect_true(fit$theta_se > 0)
})

test_that("NB full-sandwich SEs are finite and positive", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005)
  se <- sqrt(diag(fit$vcov))
  expect_true(all(is.finite(se)))
  expect_true(all(se > 0))
})

test_that("NB with estimated gamma: vcov_full includes theta and gamma", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = "estimate")
  p_ab <- length(coef(fit))
  # Full vcov: p_ab + theta + gamma = p_ab + 2
  expect_equal(nrow(fit$vcov_full), p_ab + 2)
  expect_equal(ncol(fit$vcov_full), p_ab + 2)
  # Coefficient vcov still p_ab x p_ab
  expect_equal(nrow(fit$vcov), p_ab)
})

test_that("NB cluster-robust sandwich with theta works", {
  d <- small_data()
  fit <- quick_fit(d, method = "nb", gamma = 0.005,
                   cluster = ~country, vcov = "HC1")
  expect_true(all(diag(fit$vcov) > 0))
  expect_true(!is.null(fit$theta_se))
  expect_true(fit$theta_se > 0)
})

test_that("Poisson does NOT have theta_se or score_full", {
  d <- small_data()
  fit <- quick_fit(d, method = "poisson", gamma = 0.005)
  expect_null(fit$theta_se)
  expect_null(fit$score_full)
})
