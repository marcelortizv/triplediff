# =============================================================================
# Tests that the per-unit cluster vector built in att_gt() is CORRECTLY ALIGNED
# to the rows of inf_func_mat -- the porting guide calls this "the crux".
#
# Why this file exists (adversarial gap found in the original suite):
#   - CHECK 2 / CHECK 5b compute their "manual" CRVE from o$cluster_vector, the
#     object's OWN field, so they are self-referential: a deterministically
#     WRONG cluster vector passes them (implementation and check share the error).
#   - CHECK 7 only tests shuffle-INVARIANCE (stability), which a deterministic
#     misalignment satisfies, so it cannot detect a wrong-but-stable mapping.
#   - The only original test that catches a scrambled id->cluster map is the
#     analytic-vs-bootstrap CHECK 1, which is skip_on_cran() with a 10% tolerance.
#
# These tests pin the reported SE against an INDEPENDENT ground-truth cluster
# vector, reconstructed from the user's known id->cl map joined onto the unit
# ordering carried by first_period_dta (built in preprocessing, NOT by the
# att_gt() cluster_vec line under test). They fail if att_gt() misaligns,
# permutes, or otherwise corrupts the cluster vector.
# =============================================================================

.fit_clustered <- function(d, cg = "nevertreated") {
  suppressWarnings(ddd(
    yname = "y", tname = "time", idname = "id", gname = "state",
    pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
    data = d, control_group = cg, base_period = "varying",
    est_method = "dr", cluster = "cl", boot = FALSE))
}

# independent ground-truth cluster vector aligned to inf_func_mat rows:
# join the KNOWN id->cl map onto first_period_dta$id (unit row order).
.truth_cluster_vec <- function(o, clmap) {
  clmap$cl[match(o$first_period_dta$id, clmap$id)]
}
.crve <- function(IF, cv, n) sqrt(diag(crossprod(rowsum(IF, cv))) / n) / sqrt(n)

.make_data <- function(seed = 123, size = 600, n_clusters = 25, unbalanced = FALSE) {
  set.seed(seed)
  d <- gen_dgp_mult_periods(size = size, dgp_type = 1)[["data"]]
  ids <- sort(unique(d$id))
  if (unbalanced) {
    set.seed(seed + 1)
    cl <- sample(seq_len(n_clusters), length(ids), replace = TRUE,
                 prob = c(rep(8, 3), rep(1, n_clusters - 3)))
  } else {
    cl <- rep(seq_len(n_clusters), length.out = length(ids))
  }
  clmap <- data.frame(id = ids, cl = cl)
  list(data = merge(d, clmap, by = "id"), clmap = clmap)
}

# -----------------------------------------------------------------------------
# (1) cluster_vector and SE match an INDEPENDENT ground-truth partition.
#     Kills: misaligned / permuted / scrambled cluster vector in att_gt().
# -----------------------------------------------------------------------------
test_that("att_gt() cluster vector is aligned to the user's id->cl map (balanced & unbalanced)", {
  for (unb in c(FALSE, TRUE)) {
    dd  <- .make_data(seed = 123, size = 600, n_clusters = 25, unbalanced = unb)
    o   <- .fit_clustered(dd$data)
    IF  <- as.matrix(o$inf_func_mat); n <- o$n
    truth <- .truth_cluster_vec(o, dd$clmap)

    # the object's cluster vector must equal the independent ground truth
    expect_identical(as.numeric(o$cluster_vector), as.numeric(truth))

    # and the reported SE must be the CRVE of that independent partition
    manual <- .crve(IF, truth, n)
    ok <- is.finite(o$se) & is.finite(manual)
    expect_equal(as.numeric(o$se[ok]), as.numeric(manual[ok]), tolerance = 1e-8)
  }
})

# -----------------------------------------------------------------------------
# (2) row-permutation: reported SE matches the SAME independent ground truth
#     (not merely "is invariant"). Kills order-dependent misalignment that a
#     pure invariance check would miss.
# -----------------------------------------------------------------------------
test_that("clustered SE matches the ground-truth CRVE even after input rows are shuffled", {
  dd <- .make_data(seed = 55, size = 600, n_clusters = 25, unbalanced = TRUE)
  set.seed(777)
  d_shuf <- dd$data[sample(seq_len(nrow(dd$data))), ]

  o   <- .fit_clustered(d_shuf)
  IF  <- as.matrix(o$inf_func_mat); n <- o$n
  truth <- .truth_cluster_vec(o, dd$clmap)         # independent of row order

  expect_identical(as.numeric(o$cluster_vector), as.numeric(truth))
  manual <- .crve(IF, truth, n)
  ok <- is.finite(o$se) & is.finite(manual)
  expect_equal(as.numeric(o$se[ok]), as.numeric(manual[ok]), tolerance = 1e-8)
})

# -----------------------------------------------------------------------------
# (3) sensitivity: relabeling clusters (a bijection on labels) must NOT change
#     the SE; reassigning units to a DIFFERENT partition MUST change it. Guards
#     against the SE being accidentally independent of the clustering.
# -----------------------------------------------------------------------------
test_that("clustered SE is invariant to cluster relabeling but sensitive to repartitioning", {
  dd <- .make_data(seed = 7, size = 600, n_clusters = 20, unbalanced = TRUE)
  base <- .fit_clustered(dd$data)

  # (a) relabel: cl -> (n_clusters + 1 - cl). Same partition, different labels.
  d_relabel <- dd$data
  d_relabel$cl <- (max(dd$clmap$cl) + 1L) - d_relabel$cl
  o_relabel <- .fit_clustered(d_relabel)
  ok <- is.finite(base$se) & is.finite(o_relabel$se)
  expect_equal(as.numeric(base$se[ok]), as.numeric(o_relabel$se[ok]), tolerance = 1e-8)

  # (b) repartition: give every unit its own cluster -> clustered SE collapses
  #     to the i.i.d. analytic SE and DIFFERS from the grouped-cluster SE.
  d_singleton <- dd$data
  d_singleton$cl <- d_singleton$id
  o_singleton <- .fit_clustered(d_singleton)
  o_iid <- ddd(
    yname = "y", tname = "time", idname = "id", gname = "state",
    pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
    data = dd$data, control_group = "nevertreated", base_period = "varying",
    est_method = "dr", cluster = NULL, boot = FALSE)
  ok2 <- is.finite(o_singleton$se) & is.finite(o_iid$se)
  expect_equal(as.numeric(o_singleton$se[ok2]), as.numeric(o_iid$se[ok2]), tolerance = 1e-8)

  ok3 <- is.finite(base$se) & is.finite(o_singleton$se)
  expect_gt(max(abs(base$se[ok3] - o_singleton$se[ok3]) / o_singleton$se[ok3]), 0.02)
})
