# ddd(est_method = "dr") with covariates, repeated cross-sections and unbalanced panels:
# the influence function is the derivative of the estimate with respect to the unit weights.
# Designs (fx_dgp_rcs, fx_dgp_panel in helper.R): misspecified propensity score, error SD 4 in the pre-period comparison cells.
# Every derivative check uses fx_rc_deriv(): the test passes when the gap is below 1e-6 * max(1, |mean(h * psi)|), ratio < 1.

# compute_did_rc() through att_dr_rc(): two periods

test_that("two-period repeated cross-section: ATT influence function is the derivative of the estimate", {
  d <- fx_dgp_rcs(2000, 20261023)
  r <- fx_rc_deriv(d, function(x) fx_rc_fit(x, "dr", "rcs2", w = "wt"))
  expect_lt(max(r$ratio), 1)
})

test_that("two-period unbalanced panel (routed to the cross-section code): influence function is the derivative in unit weights", {
  d <- fx_dgp_panel(2000, 20261024, drop = 0.15)
  expect_length(fx_rc_fit(d, "dr", "unb2")$att_inf_func, length(unique(d$id)))
  r <- fx_rc_deriv(d, function(x) fx_rc_fit(x, "dr", "unb2", w = "wt"))
  expect_lt(max(r$ratio), 1)
})

test_that("single-observation derivative of the two-period repeated cross-section fit equals the influence function", {
  d <- fx_dgp_rcs(2000, 20261023); n <- nrow(d); sub <- fx_rc_sub(d); h <- 0.01
  att_w <- function(w) { d$wt <- w; fx_rc_fit(d, "dr", "rcs2", w = "wt")$ATT }
  psi <- as.numeric(fx_rc_fit(transform(d, wt = 1), "dr", "rcs2", w = "wt")$att_inf_func)
  for (s in 1:4) for (tt in 1:2) for (i in head(which(sub == s & d$time == tt), 2)) {
    w <- rep(1, n); w[i] <- 1 + h; ap <- att_w(w); w[i] <- 1 - h; am <- att_w(w)
    num <- (ap - am) / (h / (n + h) + h / (n - h))
    expect_lt(abs(num - psi[i]), 1e-5 * max(1, abs(psi[i])))
  }
})

test_that("two-period trimming: influence function is the derivative with the trimming indicator held fixed", {
  d <- fx_rc_trim_design(fx_dgp_rcs(2000, 20261031))
  set.seed(99); w <- 0.5 + runif(nrow(d))
  for (cnt in list(fx_rc_trim_count(d), fx_rc_trim_count(d, w))) {
    expect_true(all(cnt[, "trimmed"] >= 3))      # comparison rows are trimmed in all three propensity score models
    expect_true(all(cnt[, "near"] == 0))         # no row is close enough to 0.995 to flip in a derivative step
  }
  r <- fx_rc_deriv(d, function(x) fx_rc_fit(x, "dr", "rcs2", w = "wt"))
  expect_lt(max(r$ratio), 1)
})

# att_gt(): several periods, never-treated comparison group

test_that("multiple-period repeated cross-section, never-treated: cell influence functions are derivatives", {
  skip_on_cran()
  d <- fx_dgp_rcs(3000, 20261026, T = 4, ft = c(0, 3))
  for (bp in c("varying", "universal")) {
    cells <- data.frame(g = 3, t = c(3, 4, if (bp == "varying") 2 else 1))      # two post cells and one pre cell
    r <- fx_rc_deriv(d, function(x) fx_rc_fit(x, "dr", "mp_rcs", w = "wt", bp = bp), fx_rc_est(cells), fx_rc_psi(cells))
    expect_lt(max(r$ratio), 1)
  }
})

test_that("multiple-period unbalanced panel, never-treated: cell influence functions are derivatives", {
  skip_on_cran()
  d <- fx_dgp_panel(2000, 20261027, T = 4, ft = c(0, 3), drop = 0.2)
  for (bp in c("varying", "universal")) {
    cells <- data.frame(g = 3, t = c(3, 4, if (bp == "varying") 2 else 1))
    r <- fx_rc_deriv(d, function(x) fx_rc_fit(x, "dr", "mp_unb", w = "wt", bp = bp), fx_rc_est(cells), fx_rc_psi(cells))
    expect_lt(max(r$ratio), 1)
  }
})

# att_gt(): not-yet-treated (GMM) cells. The cell influence function holds the estimated GMM weights fixed by design,
# so only its components are checked against derivatives, and the cell is rebuilt by hand from the components.

test_that("not-yet-treated cells: components are exact derivatives and the cell is their GMM combination", {
  d <- fx_dgp_rcs(4000, 20261029, T = 5, ft = c(0, 3, 5))
  for (tt in 3:4) for (cc in c(0, 5)) {
    sub <- d[d$state %in% c(3, cc) & d$time %in% c(2, tt), ]
    sub$state <- as.numeric(sub$state == 3); sub$time <- ifelse(sub$time == 2, 1, 2)
    r <- fx_rc_deriv(sub, function(x) fx_rc_fit(x, "dr", "rcs2", w = "wt"))
    expect_lt(max(r$ratio), 1)
  }
  o <- fx_rc_fit(d, "dr", "mp_rcs", cg = "notyettreated")
  IF <- as.matrix(o$inf_func_mat)
  dS <- data.table::as.data.table(d); data.table::setnames(dS, c("state", "partition"), c("S", "Q"))
  for (tt in 3:4) {
    j <- which(o$groups == 3 & o$periods == tt)
    h <- fx_hand_gmm(dS, 3, tt, 2, c(0, 5), ~x1 + x2, "dr", panel = FALSE)
    expect_equal(o$ATT[j], h$att, tolerance = 1e-8)
    expect_equal(o$se[j], h$se, tolerance = 1e-8)
    expect_equal(as.numeric(IF[, j]), h$IF, tolerance = 1e-8)
  }
})

# agg_ddd(): overall ATT

test_that("agg_ddd overall ATT: simple and event-study influence functions are derivatives", {
  skip_on_cran()
  d <- fx_dgp_rcs(3000, 20261026, T = 4, ft = c(0, 3))
  for (type in c("simple", "eventstudy")) {
    inf_name <- if (type == "simple") "simple.att" else "dynamic.inf.func"
    r <- fx_rc_deriv(d, function(x) agg_ddd(fx_rc_fit(x, "dr", "mp_rcs", w = "wt", bp = "universal"), type = type),
                     est = function(a) a$aggte_ddd$overall.att,
                     psi = function(a) matrix(as.numeric(a$aggte_ddd$inf.function[[inf_name]]), ncol = 1))
    expect_lt(max(r$ratio), 1)
  }
})

# Controls: these pass with and without the fix.

test_that("controls: reg, ipw, dr without covariates and balanced-panel dr are exact derivatives", {
  skip_on_cran()
  d2r <- fx_dgp_rcs(2000, 20261023); d2u <- fx_dgp_panel(2000, 20261024, drop = 0.15); d2p <- fx_dgp_panel(2000, 20261025)
  dmr <- fx_dgp_rcs(3000, 20261026, T = 4, ft = c(0, 3)); dmu <- fx_dgp_panel(2000, 20261027, T = 4, ft = c(0, 3), drop = 0.2)
  dmp <- fx_dgp_panel(2000, 20261028, T = 4, ft = c(0, 3))
  cells <- data.frame(g = 3, t = c(3, 4, 2))
  chk <- function(d, est, path, xf = ~x1 + x2) {
    multi <- path %in% c("mp_rcs", "mp_unb", "mp_pan")
    r <- fx_rc_deriv(d, function(x) fx_rc_fit(x, est, path, w = "wt", xf = xf),
                     if (multi) fx_rc_est(cells) else fx_rc_est(), if (multi) fx_rc_psi(cells) else fx_rc_psi(), w0s = "rand")
    expect_true(max(r$ratio) < 1, info = paste(path, est, deparse(xf)))
  }
  for (e in c("reg", "ipw")) chk(d2r, e, "rcs2")
  chk(d2r, "dr", "rcs2", ~1)
  chk(d2u, "reg", "unb2"); chk(d2u, "dr", "unb2", ~1)
  chk(d2p, "dr", "pan2")
  for (e in c("reg", "ipw")) chk(dmr, e, "mp_rcs")
  chk(dmr, "dr", "mp_rcs", ~1)
  chk(dmu, "reg", "mp_unb")
  chk(dmp, "dr", "mp_pan")
})
