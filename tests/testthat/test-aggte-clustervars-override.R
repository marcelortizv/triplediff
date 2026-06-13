# =============================================================================
# Tests for analytic cluster-robust SEs in agg_ddd()/compute_aggregation() and
# the warn + i.i.d.-fallback behaviour when clustering is requested on an object
# that was NOT clustered on the requested variable.
#
# Note on reachability:
#   An earlier revision of run_preprocess_multPeriods() forced boot = TRUE (and
#   cband = TRUE) for ddd(..., cluster = "cl", boot = FALSE); that override has
#   since been REMOVED. The aggregation-level analytic SEs are exercised by
#   calling agg_ddd(..., boot = FALSE, cband = FALSE), which recomputes SEs from
#   the stored influence function independently of the ddd object's boot flag.
# =============================================================================

make_mult_data_with_clusters <- function(seed = 123, size = 500, n_clusters = 25,
                                         add_other = FALSE) {
  set.seed(seed)
  d <- gen_dgp_mult_periods(size = size, dgp_type = 1)[["data"]]
  ids <- sort(unique(d$id))
  clmap <- data.frame(id = ids, cl = rep(seq_len(n_clusters), length.out = length(ids)))
  if (add_other) clmap$other <- rep(seq_len(5), length.out = length(ids))
  merge(d, clmap, by = "id")
}

fit_ddd <- function(d, cluster = NULL) {
  suppressWarnings(ddd(
    yname = "y", tname = "time", idname = "id", gname = "state",
    pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
    data = d, control_group = "nevertreated", base_period = "varying",
    est_method = "dr", cluster = cluster, boot = FALSE))
}

# -----------------------------------------------------------------------------
# CHECK 3-of-aggregation flavour: analytic clustered aggregation SE equals the
# manual cluster-sum CRVE of the aggregated influence function, for every
# aggregation type. (This is the aggregation-level analogue of (g,t) check 2.)
# -----------------------------------------------------------------------------
test_that("agg_ddd analytic clustered SE == manual cluster-sum CRVE (all types)", {
  d <- make_mult_data_with_clusters(seed = 123, size = 500, n_clusters = 25)
  obj <- fit_ddd(d, cluster = "cl")
  cv  <- obj$cluster_vector
  expect_equal(obj$cluster_var, "cl")
  expect_false(is.null(cv))

  manual_crve <- function(ifu) {
    ifu <- as.numeric(ifu)
    n <- length(ifu)
    sqrt(sum(rowsum(ifu, cv)^2)) / n
  }

  # simple
  a <- agg_ddd(obj, type = "simple", boot = FALSE, cband = FALSE)$aggte_ddd
  expect_equal(a$overall.se, manual_crve(a$inf.function$simple.att), tolerance = 1e-8)

  # eventstudy (overall + per-event)
  es <- agg_ddd(obj, type = "eventstudy", boot = FALSE, cband = FALSE)$aggte_ddd
  expect_equal(es$overall.se, manual_crve(es$inf.function$dynamic.inf.func), tolerance = 1e-8)
  ife <- es$inf.function$dynamic.inf.func.e
  man_e <- apply(as.matrix(ife), 2, manual_crve)
  ok_e <- !is.na(es$se.egt) & !is.na(man_e)
  expect_equal(as.numeric(es$se.egt[ok_e]), as.numeric(man_e[ok_e]), tolerance = 1e-8)

  # group
  gr <- agg_ddd(obj, type = "group", boot = FALSE, cband = FALSE)$aggte_ddd
  expect_equal(gr$overall.se, manual_crve(gr$inf.function$selective.inf.func), tolerance = 1e-8)

  # calendar
  ca <- agg_ddd(obj, type = "calendar", boot = FALSE, cband = FALSE)$aggte_ddd
  expect_equal(ca$overall.se, manual_crve(ca$inf.function$calendar.inf.func), tolerance = 1e-8)
})

# -----------------------------------------------------------------------------
# Companion: analytic i.i.d. aggregation SE (no clustering) equals the i.i.d.
# formula sqrt(mean(IF^2)/n) from the aggregated influence function.
# -----------------------------------------------------------------------------
test_that("agg_ddd analytic i.i.d. SE == sqrt(mean(IF^2)/n) (no clustering)", {
  d <- make_mult_data_with_clusters(seed = 321, size = 500, n_clusters = 25)
  obj <- fit_ddd(d, cluster = NULL)
  expect_true(is.null(obj$cluster_vector))

  manual_iid <- function(ifu) {
    ifu <- as.numeric(ifu); n <- length(ifu)
    sqrt(mean(ifu^2) / n)
  }

  a <- agg_ddd(obj, type = "simple", boot = FALSE, cband = FALSE)$aggte_ddd
  expect_equal(a$overall.se, manual_iid(a$inf.function$simple.att), tolerance = 1e-9)

  es <- agg_ddd(obj, type = "eventstudy", boot = FALSE, cband = FALSE)$aggte_ddd
  expect_equal(es$overall.se, manual_iid(es$inf.function$dynamic.inf.func), tolerance = 1e-9)
})

# -----------------------------------------------------------------------------
# CHECK 6(a): clustering requested on an object clustered on a DIFFERENT variable
# must warn and fall back to i.i.d. SEs (no error, no silent mis-report).
# -----------------------------------------------------------------------------
test_that("agg_ddd warns and falls back to i.i.d. when cluster var differs", {
  d <- make_mult_data_with_clusters(seed = 7, size = 500, n_clusters = 25, add_other = TRUE)
  obj <- fit_ddd(d, cluster = "cl")          # clustered on 'cl'
  expect_equal(obj$cluster_var, "cl")

  # requesting a DIFFERENT cluster variable must warn
  expect_warning(
    agg_ddd(obj, type = "simple", cluster = "other", boot = FALSE, cband = FALSE),
    regexp = "cannot switch|do NOT account|not available"
  )

  ovr <- suppressWarnings(
    agg_ddd(obj, type = "simple", cluster = "other", boot = FALSE, cband = FALSE)$aggte_ddd)

  # fallback SEs must equal the i.i.d. analytic formula (NOT a cluster-robust value)
  ifu <- as.numeric(ovr$inf.function$simple.att); n <- length(ifu)
  expect_equal(ovr$overall.se, sqrt(mean(ifu^2) / n), tolerance = 1e-9)

  # and the returned argu reflects the fallback to NULL clustering
  expect_true(is.null(ovr$argu$cluster))
})

# -----------------------------------------------------------------------------
# CHECK 6(b): clustering requested on an object built with cluster = NULL must
# warn and fall back to i.i.d. SEs; the fallback SEs equal the no-cluster SEs.
# -----------------------------------------------------------------------------
test_that("agg_ddd warns and falls back when ddd object was not clustered", {
  d <- make_mult_data_with_clusters(seed = 11, size = 500, n_clusters = 25)
  obj <- fit_ddd(d, cluster = NULL)          # NOT clustered
  expect_true(is.null(obj$cluster_var))

  expect_warning(
    agg_ddd(obj, type = "eventstudy", cluster = "cl", boot = FALSE, cband = FALSE),
    regexp = "not supplied|do NOT account|not available"
  )

  # i.i.d. baseline: no cluster requested at all
  base <- agg_ddd(obj, type = "eventstudy", boot = FALSE, cband = FALSE)$aggte_ddd
  ovr  <- suppressWarnings(
    agg_ddd(obj, type = "eventstudy", cluster = "cl", boot = FALSE, cband = FALSE)$aggte_ddd)

  expect_equal(ovr$overall.se, base$overall.se, tolerance = 1e-10)
  ok <- !is.na(base$se.egt) & !is.na(ovr$se.egt)
  expect_equal(as.numeric(ovr$se.egt[ok]), as.numeric(base$se.egt[ok]), tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# CHECK 1 (aggregation flavour): analytic clustered aggregation SE approx the
# clustered bootstrap aggregation SE (large nboot). Monte-Carlo tolerance.
# This path is reachable (aggregation honours its own boot flag), so it is run.
# -----------------------------------------------------------------------------
test_that("agg_ddd analytic clustered SE approx bootstrap clustered SE", {
  skip_on_cran()
  d <- make_mult_data_with_clusters(seed = 2025, size = 800, n_clusters = 40)
  obj <- fit_ddd(d, cluster = "cl")

  for (ty in c("simple", "eventstudy")) {
    a <- agg_ddd(obj, type = ty, boot = FALSE, cband = FALSE)$aggte_ddd
    set.seed(99)
    b <- suppressWarnings(
      agg_ddd(obj, type = ty, boot = TRUE, nboot = 5000, cband = FALSE)$aggte_ddd)
    rel <- abs(a$overall.se - b$overall.se) / a$overall.se
    expect_lt(rel, 0.10)
  }
})
