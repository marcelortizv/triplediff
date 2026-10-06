# Item 9 (0.3.0): multi-period repeated cross-sections without an id work with a cluster variable.

test_that("multi-period RCS without an id runs with cluster and equals the unique-row-id fit (item 9)", {
  d <- fx_mp_base(); d$id <- NULL
  f <- function(dd, idname, ...) ddd(yname = "y", tname = "time", idname = idname, gname = "state", pname = "partition",
                                     xformla = ~cov1 + cov2, data = dd, control_group = "notyettreated",
                                     base_period = "varying", est_method = "dr", panel = FALSE, ...)
  o <- f(d, NULL, cluster = "cluster")
  d2 <- d; d2$rid <- seq_len(nrow(d2))
  r <- f(d2, "rid", cluster = "cluster")
  expect_equal(o$ATT, r$ATT, tolerance = 1e-12)
  expect_equal(o$se, r$se, tolerance = 1e-12)
  dt <- data.table::as.data.table(d)
  f(dt, NULL, cluster = "cluster")
  expect_false(".row_id" %in% names(dt))
  expect_no_error(f(dt, NULL))
})
