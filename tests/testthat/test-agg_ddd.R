# agg_ddd(): cohorts are weighted by P(S = g, Q = 1) in the aggregations.

fx_agg_fit <- function() {
  d <- fx_dgp_a(1500, 20261009)
  ddd(yname = "y", tname = "time", idname = "id", gname = "S", pname = "Q", xformla = ~x1 + x2, data = d,
      control_group = "nevertreated", base_period = "varying", est_method = "dr")
}

test_that("aggregate estimates use the share of treated units as cohort weights", {
  o <- fx_agg_fit(); fp <- o$first_period_dta
  pg <- sapply(o$glist, function(g) mean(fp$weights * (fp$first_treat == g & fp$partition == 1)))
  w_of <- function(gs) pg[match(gs, o$glist)]
  G <- o$groups; Tt <- o$periods; A <- o$ATT
  es <- agg_ddd(o, type = "eventstudy")$aggte_ddd
  for (i in seq_along(es$egt)) {
    k <- which(Tt - G == es$egt[i])
    expect_equal(es$att.egt[i], sum(w_of(G[k]) * A[k]) / sum(w_of(G[k])), tolerance = 1e-10)
  }
  k <- which(Tt >= G)
  expect_equal(agg_ddd(o, type = "simple")$aggte_ddd$overall.att, sum(w_of(G[k]) * A[k]) / sum(w_of(G[k])), tolerance = 1e-10)
  th <- sapply(o$glist, function(g) mean(A[G == g & Tt >= g]))
  expect_equal(agg_ddd(o, type = "group")$aggte_ddd$overall.att, sum(pg * th) / sum(pg), tolerance = 1e-10)
  ca <- agg_ddd(o, type = "calendar")$aggte_ddd
  for (i in seq_along(ca$egt)) {
    k <- which(Tt == ca$egt[i] & G <= Tt)
    expect_equal(ca$att.egt[i], sum(w_of(G[k]) * A[k]) / sum(w_of(G[k])), tolerance = 1e-10)
  }
})

test_that("event-study SE uses treated-unit indicators in the weight influence function", {
  o <- fx_agg_fit(); fp <- o$first_period_dta; n <- o$n
  pg <- sapply(o$glist, function(g) mean(fp$weights * (fp$first_treat == g & fp$partition == 1)))
  G <- o$groups; Tt <- o$periods; A <- o$ATT; IF <- as.matrix(o$inf_func_mat)
  Gu <- fp$first_treat * fp$partition
  k <- which(Tt - G == 0); p <- pg[match(G[k], o$glist)]; S <- sum(p)
  ind <- sapply(G[k], function(g) fp$weights * (Gu == g)) - matrix(p, n, length(k), byrow = TRUE)
  wif <- ind / S - rowSums(ind) %o% (p / S^2)
  psi <- IF[, k] %*% (p / S) + wif %*% A[k]
  es <- agg_ddd(o, type = "eventstudy")$aggte_ddd
  expect_equal(es$se.egt[es$egt == 0], sqrt(mean(psi^2) / n), tolerance = 1e-8)
})

# agg_ddd(type = "group"): overall standard error.

test_that("group overall influence function is centered and matches a numerical derivative", {
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

# agg_ddd(): alpha is taken from the ddd() object unless alpha is given.

test_that("agg_ddd() inherits alpha from the ddd object by default", {
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
