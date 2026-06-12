# =============================================================================
# Tests for the clustered multiplier bootstrap mboot() (CS 2021, Remark 10).
#
# Feature under test:
#   mboot() now aggregates the influence function to cluster SUMS and scales the
#   SE by bSigma * sqrt(n_clusters) / n. For EQUAL-SIZED clusters of size m this
#   must coincide with the pre-change cluster-MEANS result (bSigma_mean / sqrt(nc)),
#   because cluster_sum = m * cluster_mean and the multiplier bootstrap is linear
#   in the influence-function columns:
#
#       se_sum  = bSigma_sum  * sqrt(nc) / n
#               = (m * bSigma_mean) * sqrt(nc) / (m * nc)
#               = bSigma_mean / sqrt(nc)
#               = se_mean.
#
#   These tests exercise mboot() directly on a synthetic influence function so
#   the identity is checked deterministically (fixed multiplier seed), free of
#   the multi-period preprocessing path.
# =============================================================================

# small helper to build the minimal did_preprocessed list mboot() consumes
make_dp_for_mboot <- function(cl, nboot) {
  n  <- length(cl)
  dt <- data.table::data.table(id = seq_len(n), cluster = cl)
  dt2 <- rbind(data.table::copy(dt)[, period := 1L],
               data.table::copy(dt)[, period := 2L])
  list(preprocessed_data = dt2,
       cluster = "cluster",
       cluster_vector = cl,
       nboot = nboot,
       alpha = 0.05,
       panel = TRUE,
       allow_unbalanced_panel = FALSE)
}

# IQR-based robust sd, identical to the one used inside mboot()
iqr_sd <- function(bres) {
  apply(bres, 2, function(b)
    (stats::quantile(b, .75, type = 1, na.rm = TRUE) -
       stats::quantile(b, .25, type = 1, na.rm = TRUE)) /
      (stats::qnorm(.75) - stats::qnorm(.25)))
}

# -----------------------------------------------------------------------------
# CHECK 3 (scaling identity): for equal-sized clusters of size m with nc clusters
# and n = m * nc units, the SE scaling sqrt(nc)/n equals 1/(m*sqrt(nc)). This is
# the algebraic backbone of "sum-based == mean-based" for balanced clusters.
# -----------------------------------------------------------------------------
test_that("equal-sized clusters: SE scaling sqrt(nc)/n == 1/(m*sqrt(nc))", {
  for (m in c(2L, 5L, 10L)) {
    nc <- 30L
    n  <- m * nc
    expect_equal(sqrt(nc) / n, 1 / (m * sqrt(nc)), tolerance = 1e-12)
  }
})

# -----------------------------------------------------------------------------
# CHECK 3 (main): the Remark-10 cluster-SUM mboot SE equals the cluster-MEANS SE
# for equal-sized clusters, to machine precision, using the same multiplier seed.
# -----------------------------------------------------------------------------
test_that("equal-sized clusters: cluster-SUM mboot SE == cluster-MEANS SE", {
  skip_on_cran()
  set.seed(1)
  n  <- 240L; k <- 4L; m <- 6L
  nc <- n / m
  IF <- matrix(stats::rnorm(n * k), n, k)
  cl <- rep(seq_len(nc), each = m)          # exactly equal-sized clusters

  dp <- make_dp_for_mboot(cl, nboot = 4000)

  # mboot() performs the cluster-SUM (Remark 10) computation
  set.seed(123)
  r_sum <- mboot(IF, dp)

  # manual cluster-MEANS computation with the SAME multiplier draws
  cmean <- rowsum(IF, cl) / m
  set.seed(123)
  bres_mean <- sqrt(nc) * BMisc::multiplier_bootstrap(cmean, 4000)
  se_mean   <- as.numeric(iqr_sd(bres_mean)) / sqrt(nc)

  expect_equal(as.numeric(r_sum$se), se_mean, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# CHECK 3 (cross-check): the cluster-SUM mboot SE also equals the manual
# cluster-SUM SE bSigma_sum * sqrt(nc) / n (definition of the implemented path).
# -----------------------------------------------------------------------------
test_that("cluster-SUM mboot SE matches its own definition bSigma*sqrt(nc)/n", {
  skip_on_cran()
  set.seed(2)
  n  <- 300L; k <- 3L; m <- 5L
  nc <- n / m
  IF <- matrix(stats::rnorm(n * k), n, k)
  cl <- rep(seq_len(nc), each = m)

  dp <- make_dp_for_mboot(cl, nboot = 3000)

  set.seed(321)
  r_sum <- mboot(IF, dp)

  csum <- rowsum(IF, cl)
  set.seed(321)
  bres_sum <- sqrt(nc) * BMisc::multiplier_bootstrap(csum, 3000)
  se_sum   <- as.numeric(iqr_sd(bres_sum)) * sqrt(nc) / n

  expect_equal(as.numeric(r_sum$se), se_sum, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Sanity: with one unit per cluster (n_clusters == n), the cluster bootstrap
# reduces to the i.i.d. multiplier bootstrap (sqrt(n) * MB(IF)), so the SE equals
# bSigma / sqrt(n).
# -----------------------------------------------------------------------------
test_that("one-unit-per-cluster reduces to the i.i.d. multiplier bootstrap SE", {
  skip_on_cran()
  set.seed(3)
  n  <- 250L; k <- 3L
  IF <- matrix(stats::rnorm(n * k), n, k)
  cl <- seq_len(n)                          # each unit is its own cluster

  dp <- make_dp_for_mboot(cl, nboot = 3000)

  set.seed(555)
  r_cluster <- mboot(IF, dp)

  # i.i.d. bootstrap with same seed: bres = sqrt(n) * MB(IF), se = bSigma / sqrt(n)
  set.seed(555)
  bres_iid <- sqrt(n) * BMisc::multiplier_bootstrap(IF, 3000)
  se_iid   <- as.numeric(iqr_sd(bres_iid)) / sqrt(n)

  expect_equal(as.numeric(r_cluster$se), se_iid, tolerance = 1e-10)
})
