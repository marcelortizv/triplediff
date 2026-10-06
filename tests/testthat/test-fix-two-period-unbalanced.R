# Item 5 (0.3.0): two-period panels are sorted and balanced (or routed) before the covariate check.

test_that("a balanced two-period panel gives the same result in any row order (item 5)", {
  d <- fx_2p_base(); o1 <- fx_fit_2p(d)
  set.seed(1); o2 <- fx_fit_2p(d[sample(nrow(d)), ])
  expect_identical(o2$ATT, o1$ATT); expect_identical(o2$se, o1$se)
})

test_that("an unbalanced two-period panel is balanced before the covariate check (item 5)", {
  d <- fx_2p_base(); u <- fx_2p_unbalanced(d)
  expect_warning(oa <- fx_fit_2p(u$da), "Dropped observations from 10 units")
  ra <- fx_fit_2p(d[!(d$id %in% u$A), ])
  expect_identical(oa$ATT, ra$ATT); expect_identical(oa$se, ra$se)
  expect_warning(ob <- fx_fit_2p(u$db), "Dropped observations from 20 units")
  rb <- fx_fit_2p(d[!(d$id %in% c(u$A, u$B)), ])
  expect_identical(ob$ATT, rb$ATT); expect_identical(ob$se, rb$se)
})

test_that("allow_unbalanced_panel = TRUE uses the repeated cross-section estimator (item 5)", {
  d <- fx_2p_base(); u <- fx_2p_unbalanced(d)
  expect_warning(ou <- fx_fit_2p(u$db, allow_unbalanced_panel = TRUE), "repeated cross-section estimator")
  orc <- fx_fit_2p(u$db, panel = FALSE)
  expect_equal(ou$ATT, orc$ATT, tolerance = 1e-12)
  # Unit-level standard error: the influence function is summed by unit
  dd <- data.table::as.data.table(u$db); data.table::setorderv(dd, c("id", "time"))
  nid <- data.table::uniqueN(dd$id)
  psi <- (nid / nrow(dd)) * rowsum(orc$att_inf_func, dd$id)[, 1]
  expect_equal(length(ou$att_inf_func), nid)
  expect_equal(ou$se, stats::sd(psi) / sqrt(nid), tolerance = 1e-10)
})
