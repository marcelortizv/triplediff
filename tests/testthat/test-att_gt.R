# ddd(control_group = "notyettreated"): GMM cells whose comparison set leaves out a treated cohort.

test_that("not-yet-treated GMM cells equal the GMM of two-period fits", {
  d <- fx_dgp_a(1500, 20261006)
  o <- ddd(yname = "y", tname = "time", idname = "id", gname = "S", pname = "Q", xformla = ~x1 + x2,
           data = d, control_group = "notyettreated", base_period = "varying", est_method = "dr", inffunc = TRUE)
  IF <- as.matrix(o$inf_func_mat)
  S_unit <- d[time == 1][order(id)]$S
  # (g, t, base, comparisons): (2,3) post, base 1, comparisons 0 and 4; (4,2) pre, base 1, comparisons 0 and 3
  for (cell in list(list(2, 3, 1, c(0, 4)), list(4, 2, 1, c(0, 3)))) {
    g <- cell[[1]]; t <- cell[[2]]; b <- cell[[3]]; comps <- cell[[4]]
    j <- which(o$groups == g & o$periods == t)
    h <- fx_hand_gmm(d, g, t, b, comps, ~x1 + x2, "dr")
    expect_equal(o$ATT[j], h$att, tolerance = 1e-10)
    expect_equal(o$se[j], h$se, tolerance = 1e-8)
    expect_equal(as.numeric(IF[, j]), h$IF, tolerance = 1e-8)
    expect_true(all(IF[!(S_unit %in% c(g, comps)), j] == 0))
  }
})

test_that("not-yet-treated results do not depend on how units are numbered", {
  d <- fx_dgp_a(1500, 20261007)
  f <- function(dd) ddd(yname = "y", tname = "time", idname = "id", gname = "S", pname = "Q", xformla = ~x1 + x2,
                        data = dd, control_group = "notyettreated", base_period = "varying", est_method = "dr")
  o1 <- f(d)
  set.seed(7)
  map <- sample(100000:999999, data.table::uniqueN(d$id))
  d2 <- data.table::copy(d)
  d2[, id := map[id]]
  d2 <- d2[sample(.N)]
  o2 <- f(d2)
  expect_equal(o2$ATT, o1$ATT, tolerance = 1e-10)
  expect_equal(o2$se, o1$se, tolerance = 1e-10)
})

# ddd(base_period = "universal"): pre-treatment cells use the period-t covariates for both nuisance models.

fx_c_fit <- function(d, m, base = "universal", w = NULL, cg = "nevertreated") {
  ddd(yname = "y", tname = "time", idname = "id", gname = "S", pname = "Q", xformla = ~x, data = d,
      control_group = cg, base_period = base, est_method = m, weightsname = w, inffunc = TRUE)
}

test_that("pre cell t < g - 1 equals the two-period fit with period-t covariates", {
  d <- fx_dgp_c(1000, 20261010)
  for (m in c("reg", "dr")) {
    o <- fx_c_fit(d, m); j <- which(o$groups == 3 & o$periods == 1)
    h <- fx_hand_gmm(d, g = 3, t = 1, b = 2, comps = 0, ~x, m)
    expect_equal(o$ATT[j], h$att, tolerance = 1e-10)
    expect_equal(o$se[j], h$se, tolerance = 1e-8)
    expect_equal(as.numeric(as.matrix(o$inf_func_mat)[, j]), h$IF, tolerance = 1e-8)
    d2 <- data.table::copy(d); d2[time > 1, x := W]          # later-period covariates do not matter
    o2 <- fx_c_fit(d2, m)
    expect_identical(o2$ATT[j], o$ATT[j]); expect_identical(o2$se[j], o$se[j])
  }
})

test_that("influence function of a pre cell t < g - 1 matches a numerical derivative", {
  d <- fx_dgp_c(1200, 20261011)
  ids <- sort(unique(d$id)); u <- d[time == 1][order(id)]
  set.seed(20261012)
  w0 <- 0.5 + runif(length(ids))
  H <- cbind(W = as.numeric(scale(u$W)), Z = as.numeric(scale(u$Z)), N = rnorm(length(ids)))
  for (m in c("reg", "dr")) {
    fit_w <- function(wv) { dd <- data.table::copy(d); dd[, wt := wv[match(id, ids)]]; fx_c_fit(dd, m, w = "wt") }
    o <- fit_w(w0); j <- which(o$groups == 3 & o$periods == 1)
    psi <- as.numeric(as.matrix(o$inf_func_mat)[, j])
    for (k in colnames(H)) {
      cd <- fx_num_deriv(function(wv) { f <- fit_w(wv); f$ATT[f$groups == 3 & f$periods == 1] }, w0, H[, k])
      hif <- mean(H[, k] * psi)
      expect_lt(abs(cd - hif), 1e-6 * max(1, abs(hif)))
    }
  }
})

test_that("not-yet-treated GMM pre cell uses period-t covariates", {
  d <- fx_dgp_d2(1500, 20261013)
  for (m in c("reg", "dr")) {
    o <- fx_c_fit(d, m, cg = "notyettreated"); j <- which(o$groups == 3 & o$periods == 1)
    h <- fx_hand_gmm(d, g = 3, t = 1, b = 2, comps = c(0, 4), ~x, m)
    expect_equal(o$ATT[j], h$att, tolerance = 1e-10)
    expect_equal(o$se[j], h$se, tolerance = 1e-8)
    expect_equal(as.numeric(as.matrix(o$inf_func_mat)[, j]), h$IF, tolerance = 1e-8)
  }
})

test_that("post-treatment cells are the same under both base periods", {
  d <- fx_dgp_c(1000, 20261010)
  for (m in c("reg", "dr")) {
    ou <- fx_c_fit(d, m, "universal"); ov <- fx_c_fit(d, m, "varying")
    expect_equal(ou$ATT[ou$groups == 3 & ou$periods == 3], ov$ATT[ov$groups == 3 & ov$periods == 3], tolerance = 1e-12)
  }
})
