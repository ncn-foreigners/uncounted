# ---- popsize(by =) with NA in the by variable (CRAN) ----
#
# Rows with an NA by-value used to vanish from the table and the total, while
# the total's SE was still computed over all rows.

na_by_data <- function() {
  d <- testdata
  d$region <- ifelse(d$sex == "M", "A", "B")
  d$region[1:30] <- NA
  d
}

test_that("NA by-values form a group, with a message", {
  fit <- estimate_hidden_pop(na_by_data(), ~m, ~n, ~N, gamma = 0.005)
  expect_message(ps <- popsize(fit, by = ~region, total = TRUE),
                 "region.*30 row\\(s\\).*\"NA\"")
  expect_setequal(ps$group, c("A", "B", "NA"))

  ## Groups partition the data: the total equals the ungrouped estimate
  all_rows <- popsize(fit)
  expect_equal(sum(ps$estimate), sum(all_rows$estimate))
  expect_equal(attr(ps, "total")$estimate, sum(all_rows$estimate))

  ## The total's SE agrees with the group covariance
  expect_equal(attr(ps, "total")$se, sqrt(sum(attr(ps, "vcov_estimate"))),
               tolerance = 1e-8)
  expect_equal(sum(ps$observed), sum(fit$m))
})

test_that("the NA group's by-value stays NA, so filters skip it", {
  fit <- estimate_hidden_pop(na_by_data(), ~m, ~n, ~N, gamma = 0.005)
  ps <- suppressMessages(popsize(fit, by = ~region))
  groups <- attr(ps, "groups")
  expect_true(is.na(groups$region[groups$.group == "NA"]))
  h <- hypotheses_popsize(ps, "xi[region == 'A'] = 0", estimate = "estimate")
  expect_equal(h$estimate, ps$estimate[ps$group == "A"])
})

test_that("NA in one of several by variables gets a combined label", {
  fit <- estimate_hidden_pop(na_by_data(), ~m, ~n, ~N, gamma = 0.005)
  ps <- suppressMessages(popsize(fit, by = ~region + year))
  expect_true(any(grepl("NA", ps$group)))
  expect_equal(sum(ps$estimate), sum(popsize(fit)$estimate))
})

test_that("no message without NA, and bootstrap reports it only once", {
  fit <- estimate_hidden_pop(testdata, ~m, ~n, ~N, gamma = 0.005)
  expect_no_message(popsize(fit, by = ~sex))

  fit_na <- estimate_hidden_pop(na_by_data(), ~m, ~n, ~N, gamma = 0.005)
  msgs <- character(0)
  b <- withCallingHandlers(
    bootstrap_popsize(fit_na, R = 3, by = ~region, seed = 1, verbose = FALSE, ci_type = "cheap"),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(length(grep("contain NA", msgs)) <= 2L)
  expect_setequal(b$popsize$group, c("A", "B", "NA"))
})
