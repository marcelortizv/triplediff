# N19 (0.3.0): two-period repeated cross-sections without an id never modify the user's data.table.

test_that("two-period ddd() without an id leaves the user's data.table unchanged (N19)", {
  dt <- data.table::as.data.table(fx_2p_base())
  dt[, id := NULL]
  saved <- data.table::copy(dt)
  f <- function() ddd(yname = "y", tname = "time", idname = NULL, gname = "state", pname = "partition",
                      xformla = ~cov1 + cov2, data = dt, est_method = "dr", panel = FALSE, inffunc = TRUE)
  o1 <- f()
  expect_identical(dt, saved)
  expect_false(".row_id" %in% names(dt))
  o2 <- f()
  expect_identical(dt, saved)
  expect_identical(o2$ATT, o1$ATT); expect_identical(o2$se, o1$se)
  expect_identical(o2$att_inf_func, o1$att_inf_func)
  # same fit as with a data.frame input, which ddd() always copies
  o3 <- ddd(yname = "y", tname = "time", idname = NULL, gname = "state", pname = "partition",
            xformla = ~cov1 + cov2, data = as.data.frame(saved), est_method = "dr", panel = FALSE, inffunc = TRUE)
  expect_identical(o1$ATT, o3$ATT); expect_identical(o1$se, o3$se)
})
