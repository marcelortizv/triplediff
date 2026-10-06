# Item 3 (0.3.0): aggregations weight cohorts by P(S = g, Q = 1).

fx_item3_fit <- function() {
  d <- fx_dgp_a(1500, 20261009)
  ddd(yname = "y", tname = "time", idname = "id", gname = "S", pname = "Q", xformla = ~x1 + x2, data = d,
      control_group = "nevertreated", base_period = "varying", est_method = "dr")
}

test_that("aggregate estimates use the share of treated units as cohort weights (item 3)", {
  o <- fx_item3_fit(); fp <- o$first_period_dta
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

test_that("event-study SE uses treated-unit indicators in the weight influence function (item 3)", {
  o <- fx_item3_fit(); fp <- o$first_period_dta; n <- o$n
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
