# Item 6 (0.3.0): input validation.

test_that("duplicated (id, time) rows stop with a clear error (item 6a)", {
  d <- fx_2p_base(); dd <- rbind(d, d[3, ])
  expect_error(fx_fit_2p(dd), "must be unique")
  expect_error(fx_fit_2p(dd, panel = FALSE), "must be unique")
  m <- fx_mp_base()
  expect_error(fx_fit_mp(rbind(m, m[5, ])), "must be unique")
  m2 <- m; m2$id <- NULL                    # RCS without an id: rows are the units
  expect_no_error(ddd(yname = "y", tname = "time", idname = NULL, gname = "state", pname = "partition",
                      xformla = ~cov1 + cov2, data = m2, control_group = "notyettreated",
                      base_period = "varying", est_method = "dr", panel = FALSE))
})

test_that("pname must be coded 0 and 1 (item 6b)", {
  d <- fx_2p_base(); d$partition <- d$partition + 1
  expect_error(ddd(yname = "y", tname = "time", idname = "id", gname = "state", pname = "partition",
                   xformla = ~1, data = d, est_method = "ipw"), "only two values \\(0 and 1\\)")
  m <- fx_mp_base(); m$partition <- m$partition + 1
  expect_error(fx_fit_mp(m), "only two values \\(0 and 1\\)")
})

test_that("two-period gname: Inf means never treated; other codes stop (item 6c)", {
  d <- fx_2p_base(); ref <- fx_fit_2p(d)
  dinf <- d; dinf$state[dinf$state == 0] <- Inf
  oi <- fx_fit_2p(dinf)
  expect_identical(oi$ATT, ref$ATT); expect_identical(oi$se, ref$se)
  expect_equal(fx_fit_2p(dinf, skip_data_checks = TRUE)$ATT, fx_fit_2p(d, skip_data_checks = TRUE)$ATT, tolerance = 1e-12)
  d21 <- d; d21$state <- ifelse(d21$state == 0, 2, 1)
  expect_error(fx_fit_2p(d21), "never-treated units must have gname = 0")
  dti <- data.table::as.data.table(dinf); fx_fit_2p(dti)
  expect_true(any(is.infinite(dti$state)))    # the user's data are not modified
})
