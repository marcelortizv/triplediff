# Item 8 (0.3.0): cluster = idname; weights that vary within a unit.

test_that("cluster = idname is unit-level clustering (item 8b)", {
  m <- fx_mp_base()
  o0 <- fx_fit_mp(m); o1 <- fx_fit_mp(m, cluster = "id")
  expect_identical(o1$ATT, o0$ATT); expect_identical(o1$se, o0$se); expect_null(o1$argu$cluster)
  expect_no_warning(a1 <- agg_ddd(o1, type = "eventstudy"))
  expect_identical(a1$aggte_ddd$overall.se, agg_ddd(o0, type = "eventstudy")$aggte_ddd$overall.se)
  oc <- fx_fit_mp(m, cluster = c("id", "cluster")); occ <- fx_fit_mp(m, cluster = "cluster")
  expect_identical(oc$se, occ$se); expect_identical(oc$argu$cluster, "cluster")
  d <- fx_2p_base()
  expect_identical(fx_fit_2p(d, cluster = "id")$se, fx_fit_2p(d)$se)
})

test_that("multi-period weights must be constant within a unit (item 8c)", {
  m <- fx_mp_base(); set.seed(83)
  m$w <- runif(nrow(m), .5, 1.5)
  expect_error(fx_fit_mp(m, weightsname = "w"), "must be the same across all periods")
  uw <- runif(length(unique(m$id)), .5, 1.5); m$w2 <- uw[match(m$id, sort(unique(m$id)))]
  expect_no_error(fx_fit_mp(m, weightsname = "w2"))
})
