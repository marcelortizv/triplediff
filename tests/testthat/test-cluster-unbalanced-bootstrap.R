# =============================================================================
# Tests that the clustered multiplier bootstrap actually implements the
# Callaway & Sant'Anna (2021, Remark 10) cluster-SUM aggregation on the regime
# where it MATTERS: UNBALANCED clusters.
#
# Why this file exists (adversarial gap found in the original suite):
#   test-mboot-cluster.R only exercises EQUAL-SIZED clusters, where cluster
#   means and cluster sums coincide. On equal-sized clusters the new (sum-based)
#   and old (mean-based) implementations are numerically identical, so reverting
#   Feature 2 entirely is a SURVIVING MUTANT against that file. These tests fail
#   if mboot() is reverted to cluster MEANS or to the old bSigma/sqrt(n_clusters)
#   scaling, because on unbalanced clusters sums != means.
# =============================================================================

# self-contained helpers (do NOT depend on test-mboot-cluster.R internals) -----
.make_dp_unbal <- function(cl, nboot) {
  n  <- length(cl)
  dt <- data.table::data.table(id = seq_len(n), cluster = cl)
  dt2 <- rbind(data.table::copy(dt)[, period := 1L],
               data.table::copy(dt)[, period := 2L])
  list(preprocessed_data = dt2, cluster = "cluster", cluster_vector = cl,
       nboot = nboot, alpha = 0.05, panel = TRUE, allow_unbalanced_panel = FALSE)
}
.iqr_sd <- function(bres) apply(bres, 2, function(b)
  (stats::quantile(b, .75, type = 1, na.rm = TRUE) -
     stats::quantile(b, .25, type = 1, na.rm = TRUE)) /
    (stats::qnorm(.75) - stats::qnorm(.25)))

# unbalanced cluster structure shared by the tests below
.unbalanced_clusters <- function() {
  nc <- 20L
  sizes <- c(rep(2L, 10L), rep(22L, 10L))   # 10 tiny + 10 large clusters
  rep(seq_len(nc), times = sizes)           # n = 240, deliberately lumpy
}

# -----------------------------------------------------------------------------
# PRIMARY GUARD: on UNBALANCED clusters mboot() must equal the cluster-SUM
# definition AND differ materially from the old cluster-MEAN definition.
# (Kills: revert of cluster sums -> means, and/or SE scaling -> /sqrt(nc).)
# -----------------------------------------------------------------------------
test_that("UNBALANCED clusters: mboot uses cluster SUMS (Remark 10), not means", {
  skip_on_cran()
  set.seed(101)
  cl <- .unbalanced_clusters()
  n  <- length(cl); k <- 4L; nc <- length(unique(cl))
  IF <- matrix(stats::rnorm(n * k), n, k)
  dp <- .make_dp_unbal(cl, nboot = 6000)

  set.seed(2024)
  r <- mboot(IF, dp)

  # (1) implemented SE == cluster-SUM definition  bSigma_sum * sqrt(nc) / n
  csum <- rowsum(IF, cl, reorder = TRUE)
  set.seed(2024)
  bres_sum <- sqrt(nc) * BMisc::multiplier_bootstrap(csum, 6000)
  se_sum   <- as.numeric(.iqr_sd(bres_sum)) * sqrt(nc) / n
  expect_equal(as.numeric(r$se), se_sum, tolerance = 1e-10)

  # (2) implemented SE is MATERIALLY DIFFERENT from the old cluster-MEAN form
  #     bSigma_mean / sqrt(nc)  -- this is what makes the test informative.
  cmean <- rowsum(IF, cl, reorder = TRUE) / as.vector(table(cl))
  set.seed(2024)
  bres_mean <- sqrt(nc) * BMisc::multiplier_bootstrap(cmean, 6000)
  se_mean   <- as.numeric(.iqr_sd(bres_mean)) / sqrt(nc)
  rel <- abs(as.numeric(r$se) - se_mean) / se_mean
  expect_gt(max(rel), 0.05)   # sums vs means must diverge on unbalanced clusters
})

# -----------------------------------------------------------------------------
# INTEGRATION: end-to-end ddd() on unbalanced clusters -- analytic cluster-sum
# CRVE must approximate the clustered bootstrap (which derives clusters from the
# data independently). This is the unbalanced analogue of test-cluster-analytic
# CHECK 1, run on the regime the port actually changes.
# -----------------------------------------------------------------------------
test_that("UNBALANCED clusters: analytic clustered SEs approx bootstrap clustered SEs", {
  skip_on_cran()
  set.seed(123)
  d   <- gen_dgp_mult_periods(size = 800, dgp_type = 1)[["data"]]
  ids <- sort(unique(d$id))
  set.seed(124)
  cl  <- sample(seq_len(40), length(ids), replace = TRUE,
                prob = c(rep(8, 4), rep(1, 36)))   # unbalanced
  d   <- merge(d, data.frame(id = ids, cl = cl), by = "id")

  o_a <- suppressWarnings(ddd(
    yname = "y", tname = "time", idname = "id", gname = "state",
    pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
    data = d, control_group = "nevertreated", base_period = "varying",
    est_method = "dr", cluster = "cl", boot = FALSE))
  expect_false(isTRUE(o_a$argu$boot))   # analytic path reached

  o_b <- suppressWarnings(ddd(
    yname = "y", tname = "time", idname = "id", gname = "state",
    pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
    data = d, control_group = "nevertreated", base_period = "varying",
    est_method = "dr", cluster = "cl", boot = TRUE, nboot = 6000, cband = TRUE))

  ok  <- !is.na(o_a$se) & !is.na(o_b$se)
  rel <- abs(o_a$se[ok] - o_b$se[ok]) / o_a$se[ok]
  expect_lt(stats::median(rel), 0.10)

  # and the analytic clustered SEs must NOT collapse to the i.i.d. SEs:
  o_iid <- ddd(
    yname = "y", tname = "time", idname = "id", gname = "state",
    pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
    data = d, control_group = "nevertreated", base_period = "varying",
    est_method = "dr", cluster = NULL, boot = FALSE)
  ok2 <- !is.na(o_a$se) & !is.na(o_iid$se)
  expect_gt(max(abs(o_a$se[ok2] - o_iid$se[ok2]) / o_iid$se[ok2]), 0.02)
})
