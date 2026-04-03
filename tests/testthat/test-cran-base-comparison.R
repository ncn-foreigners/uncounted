# ---- Oracle tests: uncounted vs base R (glm / lm / sandwich) ----

# ---- Poisson without gamma matches glm ----

test_that("Poisson gamma=NULL: coefficients match glm()", {
  d <- positive_data(small_data())

  glm_fit <- glm(m ~ -1 + log(N) + log(n / N), family = poisson(), data = d)
  our_fit <- quick_fit(d, method = "poisson", gamma = NULL)

  expect_equal(as.numeric(coef(our_fit)), as.numeric(coef(glm_fit)),
               tolerance = 1e-5)
})

test_that("Poisson fixed gamma: coefficients match glm()", {
  d <- small_data()
  g <- 0.005

  glm_fit <- glm(m ~ -1 + log(N) + log(g + n / N), family = poisson(), data = d)
  our_fit <- quick_fit(d, method = "poisson", gamma = g)

  expect_equal(as.numeric(coef(our_fit)), as.numeric(coef(glm_fit)),
               tolerance = 1e-5)
})

# ---- OLS matches lm ----

test_that("OLS fixed gamma: coefficients match lm()", {
  d <- positive_data(small_data())
  g <- 0.005

  lm_fit <- lm(log(m) ~ 0 + log(N) + log(g + n / N), data = d)
  our_fit <- quick_fit(d, method = "ols", gamma = g)

  expect_equal(as.numeric(coef(our_fit)), as.numeric(coef(lm_fit)),
               tolerance = 1e-5)
})

test_that("OLS with zeros uses log(m+1) and matches lm()", {
  # Use full testdata which has zeros in m
  d <- testdata
  g <- 0.005

  lm_fit <- lm(log(m + 1) ~ 0 + log(N) + log(g + n / N), data = d)
  our_fit <- quick_fit(d, method = "ols", gamma = g)

  expect_equal(as.numeric(coef(our_fit)), as.numeric(coef(lm_fit)),
               tolerance = 1e-5)
})

# ---- Covariates ----

test_that("Poisson cov_alpha=~sex: coefficients match glm interaction", {
  d <- positive_data(small_data())

  glm_fit <- glm(m ~ -1 + sex:log(N) + log(n / N), family = poisson(), data = d)
  our_fit <- quick_fit(d, method = "poisson", gamma = NULL, cov_alpha = ~0 + sex)

  glm_coefs <- coef(glm_fit)
  our_coefs <- coef(our_fit)

  expect_equal(our_coefs[["alpha:sexF"]], glm_coefs[["sexF:log(N)"]], tolerance = 1e-5)
  expect_equal(our_coefs[["alpha:sexM"]], glm_coefs[["sexM:log(N)"]], tolerance = 1e-5)
  expect_equal(our_coefs[["beta"]], glm_coefs[["log(n/N)"]], tolerance = 1e-5)
})

test_that("OLS cov_alpha+cov_beta=~sex: same structure as lm", {
  d <- positive_data(small_data())
  g <- 0.005

  lm_fit <- lm(log(m) ~ 0 + sex:log(N) + sex:log(g + n / N), data = d)
  our_fit <- quick_fit(d, method = "ols", gamma = g,
                       cov_alpha = ~0 + sex, cov_beta = ~0 + sex)

  expect_equal(length(coef(our_fit)), length(coef(lm_fit)))
})

# ---- Sandwich vcov matches ----

test_that("Poisson HC3 vcov matches sandwich::vcovHC(glm, 'HC3')", {
  d <- positive_data(small_data())

  glm_fit <- glm(m ~ -1 + log(N) + log(n / N), family = poisson(), data = d)
  our_fit <- quick_fit(d, method = "poisson", gamma = NULL, vcov = "HC3")

  V_glm <- sandwich::vcovHC(glm_fit, type = "HC3")
  V_our <- vcov(our_fit)

  expect_equal(as.numeric(V_our), as.numeric(V_glm), tolerance = 1e-4)
})

test_that("Poisson HC0 vcov matches sandwich::vcovHC(glm, 'HC0')", {
  d <- positive_data(small_data())

  glm_fit <- glm(m ~ -1 + log(N) + log(n / N), family = poisson(), data = d)
  our_fit <- quick_fit(d, method = "poisson", gamma = NULL, vcov = "HC0")

  V_glm <- sandwich::vcovHC(glm_fit, type = "HC0")
  V_our <- vcov(our_fit)

  expect_equal(as.numeric(V_our), as.numeric(V_glm), tolerance = 1e-4)
})

test_that("Poisson fixed gamma HC3: vcov matches sandwich on glm", {
  d <- small_data()
  g <- 0.005

  glm_fit <- glm(m ~ -1 + log(N) + log(g + n / N), family = poisson(), data = d)
  our_fit <- quick_fit(d, method = "poisson", gamma = g, vcov = "HC3")

  V_glm <- sandwich::vcovHC(glm_fit, type = "HC3")
  V_our <- vcov(our_fit)

  expect_equal(as.numeric(V_our), as.numeric(V_glm), tolerance = 1e-4)
})
