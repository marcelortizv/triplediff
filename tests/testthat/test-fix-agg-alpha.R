# Item 4 (0.3.0): agg_ddd() uses the ddd() alpha unless alpha is given.

test_that("agg_ddd() inherits alpha from the ddd object by default (item 4)", {
  set.seed(1234)
  dat <- gen_dgp_mult_periods(size = 500, dgp_type = 1)[["data"]]
  o <- ddd(yname = "y", tname = "time", idname = "id", gname = "state", pname = "partition",
           xformla = ~cov1 + cov2 + cov3 + cov4, data = dat, control_group = "nevertreated",
           base_period = "varying", est_method = "dr", alpha = 0.10)
  a <- agg_ddd(o, type = "eventstudy")$aggte_ddd
  expect_equal(a$argu$alpha, 0.10)
  expect_equal(unname(a$crit.val.egt), qnorm(0.95))
  a2 <- agg_ddd(o, type = "eventstudy", alpha = 0.05)$aggte_ddd
  expect_equal(a2$argu$alpha, 0.05)
  expect_equal(unname(a2$crit.val.egt), qnorm(0.975))
})
