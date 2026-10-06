# Item 1 (0.3.0): not-yet-treated GMM cells whose comparison set leaves out a treated cohort.

test_that("not-yet-treated GMM cells equal the GMM of two-period fits (item 1)", {
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

test_that("not-yet-treated results do not depend on how units are numbered (item 1)", {
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
