# uncounted package — development notes

## Next update: coverage improvement plan

Current overall coverage: ~75% (1785/2388 lines). Target: 85%+.

### Priority 1: helpers.R (54.81% → ~90%)

Missing tests for:
- **HC4, HC4m, HC5** branches in `.compute_sandwich_vcov()` (lines 139-150). Add tests that call `estimate_hidden_pop(..., vcov = "HC4")` etc. and verify SEs differ from HC3.
- **`.compute_nb_sandwich()` HC1 vs HC0** — add explicit test comparing HC0 and HC1 NB results.

### Priority 2: loo.R (69.89% → ~85%)

Missing tests for:
- **`summary.uncounted_loo()`** — no test calls `summary()` on a LOO result. Add `expect_output(summary(loo_res))`.
- **`plot(loo_res, type = "coef")`** — DFBETA coefficient panel plots untested. Add `expect_no_error(plot(loo_res, type = "coef"))`.
- **`compare_loo()` with `label_vars`** — the data/label_vars arguments never tested. Add test with `compare_loo(loo1, loo2, data = d, label_vars = ~country)`.

### Priority 3: bootstrap.R (78.64% → ~90%)

Missing tests for:
- **`point_estimate = "plugin"`** and **`"mean"`** — only default `"median"` tested. Add explicit test for each.
- **`total = TRUE`** in `bootstrap_popsize()` — per-replicate total summing never tested. Add test verifying total row appears in `boot$popsize`.
- **`by` parameter** — stratified bootstrap never tested. Add `bootstrap_popsize(fit, R = 29, by = ~sex)`.

### Priority 4: estimate_hidden_pop.R (86.21% → ~92%)

Missing tests for:
- **`.print_response_summary()`** for constrained models (lines 645-715). Add `expect_output(summary(fit_constrained), "Response-scale")`.
- **Summary with `total = TRUE` for constrained** — verify printed output.

### Priority 5: plot_uncounted.R (86.84% → ~95%)

Missing tests for:
- **`profile_gamma()` with `plot = TRUE`** — only tested with `plot = FALSE`. Add `expect_no_error(profile_gamma(fit, plot = TRUE))`.
- **`plot_explore()` with constrained model** — add test.

### Excluded: shiny_app.R (33.19%)

Would need `shinytest2` infrastructure. Low priority — Shiny server logic is reactive and hard to unit-test. Consider adding basic integration tests later if needed.

### Estimated new tests: ~15-20 tests, targeting 85%+ overall coverage.

## Architecture notes

- Test fixture: `tests/testthat/testdata.rds` (200 rows, DGP: alpha=0.70, beta=0.55, gamma=0.005)
- Helper: `tests/testthat/helper-data.R` — `testdata`, `small_data()`, `positive_data()`, `quick_fit()`
- CRAN tests: `test-cran-*.R` (fast, no skip)
- Extended tests: `test-extended-*.R` (skip_on_cran())
- NB sandwich: dedicated path in `.compute_nb_sandwich()`, bypasses `sandwich::vcovHC()`
- lrtest nesting: column-space QR projection in `.is_col_subspace()`
