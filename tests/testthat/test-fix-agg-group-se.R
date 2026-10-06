# Item 2 (0.3.0): overall standard error of agg_ddd(type = "group").

test_that("group overall influence function is centered and matches a numerical derivative (item 2)", {
  d <- fx_dgp_a(1200, 20261008)
  ids <- sort(unique(d$id)); u <- d[time == 1][order(id)]
  agg_at <- function(wv) {
    dd <- data.table::copy(d); dd[, wt := wv[match(id, ids)]]
    o <- ddd(yname = "y", tname = "time", idname = "id", gname = "S", pname = "Q", xformla = ~x1 + x2, data = dd,
             control_group = "nevertreated", base_period = "varying", est_method = "dr", weightsname = "wt")
    agg_ddd(o, type = "group")$aggte_ddd
  }
  w0 <- rep(1, length(ids))
  psi <- as.numeric(agg_at(w0)$inf.function$selective.inf.func)
  expect_lt(abs(mean(psi)), 1e-8 * sqrt(mean(psi^2)))
  H <- cbind(coh3 = (u$S == 3) - mean(u$S == 3), x1 = as.numeric(scale(u$x1)))
  for (k in colnames(H)) {
    cd <- fx_num_deriv(function(wv) agg_at(wv)$overall.att, w0, H[, k])
    hif <- mean(H[, k] * psi)
    expect_lt(abs(cd - hif), 1e-6 * max(1, abs(hif)))
  }
})
