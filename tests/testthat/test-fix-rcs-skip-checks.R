# N20 (0.3.0): two-period repeated cross-sections without an id work with skip_data_checks = TRUE.

fx_n20_fit <- function(dd, skip, method = "dr", ...) {
  ddd(yname = "y", tname = "time", idname = NULL, gname = "state", pname = "partition",
      xformla = ~cov1 + cov2, data = dd, est_method = method, panel = FALSE, inffunc = TRUE,
      skip_data_checks = skip, ...)
}
fx_n20_same <- function(a, b) {
  for (k in c("ATT", "se", "lci", "uci", "att_inf_func", "subgroup_counts", "argu")) {
    expect_identical(a[[k]], b[[k]], label = paste0("skip_data_checks = TRUE: ", k))
  }
}

test_that("two-period RCS without an id: skip_data_checks = TRUE equals the checked fit (N20)", {
  d <- fx_2p_base(); d$id <- NULL
  # rows in a non-sorted order: the internal row id must follow the user's row order in both paths
  set.seed(20); d <- d[sample(nrow(d)), ]; rownames(d) <- NULL
  for (m in c("dr", "reg", "ipw")) {
    fx_n20_same(fx_n20_fit(d, TRUE, m), fx_n20_fit(d, FALSE, m))
  }
  # sampling weights
  set.seed(21); d$w <- stats::runif(nrow(d), 0.5, 1.5)
  fx_n20_same(fx_n20_fit(d, TRUE, weightsname = "w"), fx_n20_fit(d, FALSE, weightsname = "w"))
})

test_that("two-period RCS without an id and skip_data_checks = TRUE leaves the user's data.table unchanged (N20)", {
  dt <- data.table::as.data.table(fx_2p_base())
  dt[, id := NULL]
  saved <- data.table::copy(dt)
  o1 <- fx_n20_fit(dt, TRUE)
  expect_identical(dt, saved)
  expect_false(".row_id" %in% names(dt))
  o2 <- fx_n20_fit(dt, TRUE)
  expect_identical(o2$ATT, o1$ATT); expect_identical(o2$se, o1$se)
  expect_identical(o2$att_inf_func, o1$att_inf_func)
  fx_n20_same(o1, fx_n20_fit(dt, FALSE))
})
