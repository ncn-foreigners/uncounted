# Auto-sourced by testthat before every test file.

# ---- Load shared test fixture ----
testdata <- readRDS(test_path("testdata.rds"))
testdata$year <- as.factor(testdata$year)
dgp <- attr(testdata, "dgp")  # list(alpha, beta, gamma)

# ---- Subset helpers ----

small_data <- function() {
  # Deterministic 20-row subset: one M and one F per 10 countries
  countries_10 <- unique(testdata$country)[1:10]
  d <- testdata[testdata$country %in% countries_10 & testdata$year == 2019, ]
  d[order(d$country, d$sex), ]
}

positive_data <- function(data = testdata) {
  data[data$m > 0 & data$n > 0, ]
}

single_cluster_data <- function(data = small_data()) {
  data$one_cluster <- "A"
  data
}

rank_deficient_country_data <- function() {
  d <- small_data()
  d$is_special <- as.integer(d$country == d$country[1])
  d
}

cov_gamma_data <- function() {
  d <- small_data()
  d$z <- rep(c(0, 1), length.out = nrow(d))
  d
}

# ---- Convenience wrapper ----

quick_fit <- function(data = small_data(), method = "poisson", gamma = 0.005, ...) {
  estimate_hidden_pop(
    data = data, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
    method = method, gamma = gamma, ...
  )
}
