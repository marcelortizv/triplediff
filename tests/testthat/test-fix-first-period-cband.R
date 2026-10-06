# Item 10 (0.3.0): first-period drop without never-treated units (b); cband without bootstrap (f).

test_that("units treated in period 1 are dropped without losing the comparison-only cohort (item 10b)", {
  d <- fx_dgp_nonever(1600, 99)
  f <- function(dd) ddd(yname = "y", tname = "time", idname = "id", gname = "S", pname = "Q", xformla = ~x1,
                        data = dd, control_group = "notyettreated", base_period = "varying", est_method = "dr")
  expect_warning(o <- f(d), "already treated in the first period")
  r <- f(d[S != 1])
  expect_identical(o$ATT, r$ATT); expect_identical(o$se, r$se)
  expect_identical(o$groups, r$groups); expect_identical(o$periods, r$periods); expect_identical(o$n, r$n)
  expect_equal(length(o$ATT), length(o$se)); expect_equal(length(o$ATT), length(o$groups))
})

test_that("cband = TRUE without the bootstrap is handled before estimation (item 10f)", {
  m <- fx_mp_base()
  expect_warning(o <- fx_fit_mp(m, cband = TRUE, boot = FALSE), "cband = TRUE requires boot = TRUE")
  expect_false(o$argu$cband)
  expect_no_error(agg_ddd(o, type = "eventstudy"))
})
