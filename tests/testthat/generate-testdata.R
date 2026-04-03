#!/usr/bin/env Rscript
# Run once to create testdata.rds:
#   Rscript tests/testthat/generate-testdata.R
#
# This script is NOT auto-sourced by testthat. The resulting .rds file
# is committed to the repository so tests never depend on cross-platform
# RNG reproducibility.

set.seed(20250403)

# ---- DGP parameters ----
alpha_true <- 0.70
beta_true  <- 0.55
gamma_true <- 0.005

# ---- Grid: 20 countries x 5 years x 2 sexes = 200 rows ----
grid <- expand.grid(
  country = sprintf("C%02d", 1:20),
  year    = 2019:2023,
  sex     = c("M", "F"),
  stringsAsFactors = FALSE
)

# Country-level base log-population (stable across years/sex)
country_base <- setNames(rnorm(20, mean = 9, sd = 1.5), sprintf("C%02d", 1:20))

# Reference population: lognormal with country base + small year trend + noise
grid$N <- as.integer(pmax(100, round(exp(
  country_base[grid$country] +
    0.02 * (grid$year - 2019) +
    rnorm(nrow(grid), 0, 0.3)
))))

# Auxiliary count
grid$n <- rpois(nrow(grid), lambda = exp(-2 + 0.7 * log(grid$N)))

# Expected count under the DGP
ratio <- grid$n / grid$N
mu <- grid$N^alpha_true * (gamma_true + ratio)^beta_true

# Observed count
grid$m <- rpois(nrow(grid), lambda = mu)

# Ensure at least some zeros in m and n (realistic for migration data)
zero_m <- which(grid$m == 0)
if (length(zero_m) < 8) {
  force_idx <- sample(which(grid$m > 0 & grid$m <= 3), 8 - length(zero_m))
  grid$m[force_idx] <- 0L
}
zero_n <- which(grid$n == 0)
if (length(zero_n) < 5) {
  candidates <- which(grid$n > 0 & grid$n <= 3)
  force_idx <- sample(candidates, min(5 - length(zero_n), length(candidates)))
  grid$n[force_idx] <- 0L
}

# Reorder columns for clarity
grid <- grid[, c("country", "year", "sex", "N", "n", "m")]

# Attach DGP parameters as attribute
attr(grid, "dgp") <- list(
  alpha = alpha_true,
  beta  = beta_true,
  gamma = gamma_true
)

# Save
out_path <- file.path(dirname(sys.frame(1)$ofile %||% "."), "testdata.rds")
if (!file.exists(out_path)) {
  out_path <- "tests/testthat/testdata.rds"
}
saveRDS(grid, out_path)
cat("Saved", nrow(grid), "rows to", out_path, "\n")
cat("DGP: alpha =", alpha_true, ", beta =", beta_true, ", gamma =", gamma_true, "\n")
cat("Zeros in m:", sum(grid$m == 0), ", zeros in n:", sum(grid$n == 0), "\n")
