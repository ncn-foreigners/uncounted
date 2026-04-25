# ---- Bootstrap exceedance tests (CRAN) ----

.make_boot_object <- function(t, groups = NULL, include_total = FALSE) {
  t <- as.matrix(t)
  if (is.null(groups)) {
    groups <- paste0("G", seq_len(ncol(t)))
  }

  pop_groups <- if (include_total) c(groups, "Total") else groups
  out <- list(
    t = t,
    t0 = rep(NA_real_, ncol(t)),
    t0_bc = rep(NA_real_, ncol(t)),
    popsize = data.frame(
      group = pop_groups,
      estimate = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      stringsAsFactors = FALSE
    ),
    popsize_full = data.frame(
      group = pop_groups,
      plugin = NA_real_,
      plugin_bc = NA_real_,
      boot_median = NA_real_,
      boot_mean = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      stringsAsFactors = FALSE
    ),
    total = if (include_total) list(plugin = NA_real_) else NULL,
    R = nrow(t),
    ci_type = "perc",
    point_estimate = "median",
    level = 0.95,
    boot_params = NULL,
    n_converged = sum(apply(t, 1, function(x) all(is.finite(x)))),
    cluster = FALSE
  )
  class(out) <- "uncounted_boot"
  out
}

test_that("exceedance_popsize computes the total exceedance probability", {
  boot <- .make_boot_object(
    t = cbind(seq(10, 100, by = 10), seq(1, 10)),
    groups = c("A", "B"),
    include_total = TRUE
  )

  exc <- exceedance_popsize(boot, threshold = 55)

  expect_s3_class(exc, "uncounted_popsize_exceedance")
  expect_identical(exc$group, "Total")
  expect_equal(exc$count, 5)
  expect_equal(exc$n_finite, 10)
  expect_equal(exc$estimate, 0.5)
})

test_that("exceedance_popsize supports named-group selection", {
  boot <- .make_boot_object(
    t = matrix(1:20, ncol = 2),
    groups = c("A", "B"),
    include_total = TRUE
  )

  exc <- exceedance_popsize(boot, threshold = 5, group = "A", direction = "below")

  expect_identical(exc$group, "A")
  expect_equal(exc$count, sum(boot$t[, 1] < 5))
  expect_equal(exc$estimate, exc$count / exc$n_finite)
})

test_that("exceedance_popsize uses the only group when no total is available", {
  boot <- .make_boot_object(t = matrix(1:10, ncol = 1), groups = "Only")

  exc <- exceedance_popsize(boot, threshold = 5)

  expect_identical(exc$group, "Only")
  expect_equal(exc$count, 5)
  expect_equal(exc$estimate, 0.5)
})

test_that("exceedance_popsize validates threshold, group, and direction", {
  boot <- .make_boot_object(
    t = matrix(1:20, ncol = 2),
    groups = c("A", "B"),
    include_total = TRUE
  )

  expect_error(exceedance_popsize(boot, threshold = NA_real_), "finite numeric")
  expect_error(exceedance_popsize(boot, threshold = c(1, 2)), "finite numeric")
  expect_error(exceedance_popsize(boot, threshold = 5, group = "Total"),
               "non-total bootstrap groups")
  expect_error(exceedance_popsize(boot, threshold = 5, group = "C"),
               "non-total bootstrap groups")
})

test_that("exceedance_popsize errors when too few finite draws remain", {
  t <- matrix(c(rep(NA_real_, 4), 1:5), ncol = 1)
  boot <- .make_boot_object(t = t, groups = "Only")

  expect_error(exceedance_popsize(boot, threshold = 2),
               "At least 10 finite bootstrap draws")
})

test_that("print.uncounted_popsize_exceedance shows a stable header", {
  boot <- .make_boot_object(t = matrix(1:10, ncol = 1), groups = "Only")
  exc <- exceedance_popsize(boot, threshold = 5)

  expect_output(print(exc), "Bootstrap exceedance probability")
})
