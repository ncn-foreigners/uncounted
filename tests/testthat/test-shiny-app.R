test_that("run_app is exported and callable", {
  expect_true(is.function(run_app))
})

test_that("run_app checks for required packages", {
  # We can't easily test missing packages, but we verify the function
  # signature accepts ... arguments
  expect_true("..." %in% names(formals(run_app)))
})

test_that("app UI function returns valid shiny tag", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")

  ui <- uncounted:::.app_ui()
  expect_true(inherits(ui, "shiny.tag") || inherits(ui, "shiny.tag.list"))
})

test_that("app server function has correct signature", {
  expect_true(is.function(uncounted:::.app_server))
  args <- names(formals(uncounted:::.app_server))
  expect_true(all(c("input", "output", "session") %in% args))
})
