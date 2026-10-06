# Tests for plot.agg_ddd() and plot.ddd() (structural: class, layers, mapped
# data, labels, colours and errors)
skip_if_not_installed("ggplot2")

#-----------------------------------------------------------------------------
# Fixtures
#-----------------------------------------------------------------------------
set.seed(1234)
dat <- gen_dgp_mult_periods(size = 500, dgp_type = 1)[["data"]]
xf  <- ~cov1 + cov2 + cov3 + cov4
fit_u <- ddd(yname = "y", tname = "time", idname = "id", gname = "state", pname = "partition",
             xformla = xf, data = dat, control_group = "nevertreated",
             base_period = "universal", est_method = "dr")
set.seed(42)
fit_b <- ddd(yname = "y", tname = "time", idname = "id", gname = "state", pname = "partition",
             xformla = xf, data = dat, control_group = "nevertreated",
             base_period = "universal", est_method = "dr",
             boot = TRUE, nboot = 199, cband = TRUE)
set.seed(1234)
d2 <- gen_dgp_2periods(size = 1000, dgp_type = 1)[["data"]]
fit_2 <- ddd(yname = "y", tname = "time", idname = "id", gname = "state", pname = "partition",
             xformla = xf, data = d2, est_method = "dr")

es  <- agg_ddd(fit_u, type = "eventstudy")
grp <- agg_ddd(fit_u, type = "group")
cal <- agg_ddd(fit_u, type = "calendar")

COL_PRE  <- "#414F5D"
COL_POST <- "#E8743B"
TOL      <- 1e-10

#-----------------------------------------------------------------------------
# Helpers
#-----------------------------------------------------------------------------
geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], "")
vrange <- function(p, axis) {
  ggplot2::ggplot_build(p)$layout$panel_params[[1]][[paste0(axis, ".range")]]
}
layer_of <- function(p, geom) ggplot2::layer_data(p, which(geoms(p) == geom))
legend_pos <- function(p) p$theme$legend.position
# Layer data of a coefficient plot in stored order (first stored cell on top).
coef_layer <- function(p, geom) {
  ld <- layer_of(p, geom)
  ld[order(-as.numeric(ld$y)), , drop = FALSE]
}

#-----------------------------------------------------------------------------
# T1: class, visibility, no printing
#-----------------------------------------------------------------------------
test_that("plot() returns a ggplot visibly and prints nothing", {
  for (obj in list(es, grp, cal, fit_u, fit_2)) {
    devs_before <- grDevices::dev.list()
    p <- plot(obj)
    expect_s3_class(p, "ggplot")
    expect_true(withVisible(plot(obj))$visible)
    expect_silent(plot(obj))
    expect_identical(grDevices::dev.list(), devs_before)
  }
})

#-----------------------------------------------------------------------------
# T2: event study (pointwise)
#-----------------------------------------------------------------------------
test_that("event-study plot has the right layers, data, colours and labels", {
  a <- es$aggte_ddd
  p <- plot(es)

  expect_identical(unname(geoms(p)), c("GeomHline", "GeomVline", "GeomErrorbar", "GeomPoint"))
  expect_equal(p$layers[[1]]$aes_params$colour, COL_PRE)

  hl <- ggplot2::layer_data(p, 1)
  expect_equal(unique(hl$yintercept), 0)
  expect_identical(unique(hl$linetype), "dashed")
  vl <- ggplot2::layer_data(p, 2)
  expect_equal(unique(vl$xintercept), -0.5)
  expect_identical(unique(vl$linetype), "dotted")

  fin <- is.finite(a$att.egt)
  pts <- ggplot2::layer_data(p, 4)
  expect_equal(pts$x, as.numeric(a$egt[fin]), tolerance = TOL)
  expect_equal(pts$y, as.numeric(a$att.egt[fin]), tolerance = TOL)
  expect_identical(pts$colour, ifelse(a$egt[fin] < 0, COL_PRE, COL_POST))

  eb <- ggplot2::layer_data(p, 3)
  z <- stats::qnorm(0.975)
  has_se <- is.finite(a$se.egt)
  expect_true(any(has_se))
  m <- match(a$egt[has_se], eb$x)
  expect_equal(eb$ymin[m], a$att.egt[has_se] - z * a$se.egt[has_se], tolerance = TOL)
  expect_equal(eb$ymax[m], a$att.egt[has_se] + z * a$se.egt[has_se], tolerance = TOL)

  expect_identical(p$labels$x, "Event time (e)")
  expect_identical(p$labels$y, "ATT")
  expect_null(p$labels$title)
  expect_identical(p$labels$subtitle, "Pointwise 95% confidence intervals (analytical SEs)")

  expect_equal(ggplot2::layer_scales(p)$x$get_breaks(), sort(unique(a$egt)))
  expect_false(identical(legend_pos(p), "none"))
  expect_identical(legend_pos(p), "bottom")
})

#-----------------------------------------------------------------------------
# T3: group and calendar
#-----------------------------------------------------------------------------
test_that("group and calendar plots are one-colour paths that reach zero", {
  for (case in list(list(obj = grp, xlab = "Group enabling period (g)"),
                    list(obj = cal, xlab = "Calendar period (t)"))) {
    a <- case$obj$aggte_ddd
    p <- plot(case$obj)
    expect_identical(unname(geoms(p)), c("GeomHline", "GeomErrorbar", "GeomPoint"))
    pts <- ggplot2::layer_data(p, 3)
    expect_true(all(pts$colour == COL_POST))
    expect_equal(pts$x, as.numeric(a$egt), tolerance = TOL)
    expect_identical(legend_pos(p), "none")
    expect_identical(p$labels$x, case$xlab)
    expect_true(all(a$att.egt > 10))
    yr <- vrange(p, "y")
    expect_true(yr[1] <= 0 && yr[2] >= 0)
    expect_equal(ggplot2::layer_scales(p)$x$get_breaks(), sort(unique(a$egt)))
  }
})

#-----------------------------------------------------------------------------
# T4: uniform band follows the aggregation's cband
#-----------------------------------------------------------------------------
test_that("uniform band uses the stored critical value; cband = FALSE is pointwise", {
  set.seed(7)
  es_b <- agg_ddd(fit_b, type = "eventstudy")
  a <- es_b$aggte_ddd
  expect_true(isTRUE(a$argu$cband))
  expect_gt(a$crit.val.egt, stats::qnorm(0.975))

  p <- plot(es_b)
  eb <- ggplot2::layer_data(p, which(geoms(p) == "GeomErrorbar"))
  has_se <- is.finite(a$se.egt)
  m <- match(a$egt[has_se], eb$x)
  expect_equal(eb$ymin[m], a$att.egt[has_se] - a$crit.val.egt * a$se.egt[has_se], tolerance = TOL)
  expect_equal(eb$ymax[m], a$att.egt[has_se] + a$crit.val.egt * a$se.egt[has_se], tolerance = TOL)
  expect_identical(p$labels$subtitle,
                   "Uniform 95% confidence band (multiplier bootstrap, 199 draws)")

  set.seed(7)
  es_p <- agg_ddd(fit_b, type = "eventstudy", cband = FALSE)
  ap <- es_p$aggte_ddd
  pp <- plot(es_p)
  expect_identical(pp$labels$subtitle,
                   "Pointwise 95% confidence intervals (multiplier bootstrap SEs, 199 draws)")
  ebp <- ggplot2::layer_data(pp, which(geoms(pp) == "GeomErrorbar"))
  has_se <- is.finite(ap$se.egt)
  m <- match(ap$egt[has_se], ebp$x)
  z <- stats::qnorm(0.975)
  expect_equal(ebp$ymin[m], ap$att.egt[has_se] - z * ap$se.egt[has_se], tolerance = TOL)
  expect_equal(ebp$ymax[m], ap$att.egt[has_se] + z * ap$se.egt[has_se], tolerance = TOL)
})

#-----------------------------------------------------------------------------
# T5: cband fallback to the pointwise critical value
#-----------------------------------------------------------------------------
test_that("cband = TRUE with a pointwise critical value is labelled pointwise", {
  set.seed(7)
  es_b <- agg_ddd(fit_b, type = "eventstudy")
  es_f <- es_b
  es_f$aggte_ddd$crit.val.egt <- stats::qnorm(0.975)
  expect_true(isTRUE(es_f$aggte_ddd$argu$cband))
  p <- plot(es_f)
  expect_true(startsWith(p$labels$subtitle, "Pointwise 95%"))
})

#-----------------------------------------------------------------------------
# T6: alpha and cluster text
#-----------------------------------------------------------------------------
test_that("subtitle follows alpha and cluster", {
  es_90 <- agg_ddd(fit_u, type = "eventstudy", alpha = 0.10)
  a <- es_90$aggte_ddd
  p <- plot(es_90)
  expect_true(startsWith(p$labels$subtitle, "Pointwise 90%"))
  eb <- ggplot2::layer_data(p, which(geoms(p) == "GeomErrorbar"))
  has_se <- is.finite(a$se.egt)
  m <- match(a$egt[has_se], eb$x)
  z <- stats::qnorm(0.95)
  expect_equal(eb$ymin[m], a$att.egt[has_se] - z * a$se.egt[has_se], tolerance = TOL)
  expect_equal(eb$ymax[m], a$att.egt[has_se] + z * a$se.egt[has_se], tolerance = TOL)

  set.seed(7)
  es_c <- agg_ddd(fit_b, type = "eventstudy")
  es_c$aggte_ddd$argu$cluster <- "cluster"
  pc <- plot(es_c)
  expect_true(endsWith(pc$labels$subtitle, ", clustered by cluster)"))
  expect_true(startsWith(pc$labels$subtitle, "Uniform 95%"))
})

test_that("interval_subtitle() formats the confidence level and draws", {
  expect_identical(interval_subtitle(FALSE, 0.05, FALSE, NULL, NULL),
                   "Pointwise 95% confidence intervals (analytical SEs)")
  expect_identical(interval_subtitle(FALSE, 0.025, FALSE, NULL, NULL),
                   "Pointwise 97.5% confidence intervals (analytical SEs)")
  expect_identical(interval_subtitle(TRUE, 0.05, TRUE, 999, "cluster"),
                   "Uniform 95% confidence band (multiplier bootstrap, 999 draws, clustered by cluster)")
  expect_identical(interval_subtitle(FALSE, 0.1, TRUE, 199, NULL),
                   "Pointwise 90% confidence intervals (multiplier bootstrap SEs, 199 draws)")
})

#-----------------------------------------------------------------------------
# T7: ylim zoom
#-----------------------------------------------------------------------------
test_that("ylim zooms the value axis, keeps zero and keeps the data", {
  p <- plot(es, ylim = c(5, 10))
  yr <- vrange(p, "y")
  expect_true(yr[1] <= 0 && yr[2] >= 10)
  pts <- ggplot2::layer_data(p, which(geoms(p) == "GeomPoint"))
  expect_equal(nrow(pts), sum(is.finite(es$aggte_ddd$att.egt)))

  expect_error(plot(es, ylim = c(1, NA)), "`ylim` must be NULL", fixed = TRUE)
  expect_error(plot(es, ylim = 1), "`ylim` must be NULL", fixed = TRUE)
  expect_error(plot(es, ylim = "a"), "`ylim` must be NULL", fixed = TRUE)
})

#-----------------------------------------------------------------------------
# T8: NA estimates
#-----------------------------------------------------------------------------
test_that("NA estimates are dropped with a warning; all NA is an error", {
  a <- es$aggte_ddd
  es_na <- es
  k <- which(a$egt == 0)
  expect_length(k, 1)
  es_na$aggte_ddd$att.egt[k] <- NA
  expect_warning(p_na <- plot(es_na), "1 estimate(s) are NA", fixed = TRUE)
  pts <- ggplot2::layer_data(p_na, which(geoms(p_na) == "GeomPoint"))
  expect_equal(nrow(pts), sum(is.finite(a$att.egt)) - 1)
  expect_false(0 %in% pts$x)

  es_all <- es
  es_all$aggte_ddd$att.egt[] <- NA
  expect_error(plot(es_all), "Every estimate is NA", fixed = TRUE)

  fit_na <- fit_u
  fit_na$ATT[length(fit_na$ATT)] <- NA
  expect_warning(pd <- plot(fit_na), "1 estimate(s) are NA", fixed = TRUE)
  expect_equal(nrow(ggplot2::layer_data(pd, 3)), length(fit_u$ATT) - 1)
})

#-----------------------------------------------------------------------------
# T9: argument and type errors; label arguments
#-----------------------------------------------------------------------------
test_that("invalid arguments and simple aggregations error clearly", {
  expect_error(plot(agg_ddd(fit_u, type = "simple")), "simple", fixed = TRUE)
  expect_error(plot(es, subtitle = 1),
               "`subtitle` must be \"auto\", NULL, or a single character string.", fixed = TRUE)
  expect_error(plot(es, title = c("a", "b")),
               "`title` must be NULL or a single character string.", fixed = TRUE)
  expect_error(plot(es, xlab = NA),
               "`xlab` must be NULL or a single character string.", fixed = TRUE)
  expect_error(plot(fit_u, ylab = 1),
               "`ylab` must be NULL or a single character string.", fixed = TRUE)
  expect_error(plot(fit_u, xlim = c(0, Inf)),
               "`xlim` must be NULL or a numeric vector of length 2 with finite values.",
               fixed = TRUE)

  es_bad <- es
  es_bad$aggte_ddd$type <- "other"
  expect_error(plot(es_bad), "plot() supports aggregations of type", fixed = TRUE)

  p <- plot(es, subtitle = NULL)
  expect_null(p$labels$subtitle)
  p <- plot(es, title = "T", subtitle = "S", xlab = "", ylab = NULL)
  expect_identical(p$labels$title, "T")
  expect_identical(p$labels$subtitle, "S")
  expect_identical(p$labels$x, "")
  expect_null(p$labels$y)
})

#-----------------------------------------------------------------------------
# T10: event study without pre-treatment periods
#-----------------------------------------------------------------------------
test_that("event study with only post-treatment periods has no separator or legend", {
  es0 <- agg_ddd(fit_u, type = "eventstudy", min_e = 0)
  p <- plot(es0)
  expect_false("GeomVline" %in% geoms(p))
  expect_identical(legend_pos(p), "none")
  pts <- ggplot2::layer_data(p, which(geoms(p) == "GeomPoint"))
  expect_true(all(pts$colour == COL_POST))
})

#-----------------------------------------------------------------------------
# T11: plot.ddd, multiple periods
#-----------------------------------------------------------------------------
test_that("plot.ddd draws one row per ATT(g,t) in stored order", {
  p <- plot(fit_u)
  expect_identical(unname(geoms(p)), c("GeomVline", "GeomLinerange", "GeomPoint"))
  vl <- ggplot2::layer_data(p, 1)
  expect_equal(unique(vl$xintercept), 0)
  expect_identical(unique(vl$linetype), "dashed")

  pts <- coef_layer(p, "GeomPoint")
  expect_equal(pts$x, as.numeric(fit_u$ATT), tolerance = TOL)

  lr <- coef_layer(p, "GeomLinerange")
  expect_equal(lr$xmin, as.numeric(fit_u$lci), tolerance = TOL)
  expect_equal(lr$xmax, as.numeric(fit_u$uci), tolerance = TOL)

  lims <- ggplot2::ggplot_build(p)$layout$panel_scales_y[[1]]$get_limits()
  expect_identical(lims, rev(paste0("ATT(", fit_u$groups, ",", fit_u$periods, ")")))

  raw <- ggplot2::layer_data(p, 3)
  expect_equal(raw$x[which.max(as.numeric(raw$y))], as.numeric(fit_u$ATT[1]), tolerance = TOL)

  expect_identical(pts$colour, ifelse(fit_u$periods < fit_u$groups, COL_PRE, COL_POST))
  expect_identical(p$labels$x, "Estimate")
  expect_null(p$labels$y)
  xr <- vrange(p, "x")
  expect_true(xr[1] <= 0 && xr[2] >= 0)
  expect_identical(p$labels$subtitle, "Pointwise 95% confidence intervals (analytical SEs)")
  expect_identical(legend_pos(p), "bottom")

  ref <- which(is.na(fit_u$se))
  expect_true(length(ref) > 0)
  expect_equal(pts$x[ref], as.numeric(fit_u$ATT[ref]), tolerance = TOL)
  expect_true(all(is.na(lr$xmin[ref])))
  expect_true(all(is.na(lr$xmax[ref])))
})

#-----------------------------------------------------------------------------
# T12: plot.ddd with cband, and far from zero
#-----------------------------------------------------------------------------
test_that("plot.ddd uses the stored band and always shows zero", {
  p <- plot(fit_b)
  lr <- coef_layer(p, "GeomLinerange")
  expect_equal(lr$xmin, as.numeric(fit_b$lci), tolerance = TOL)
  expect_equal(lr$xmax, as.numeric(fit_b$uci), tolerance = TOL)
  expect_identical(p$labels$subtitle,
                   "Uniform 95% confidence band (multiplier bootstrap, 199 draws)")

  far <- fit_u
  keep <- is.finite(far$ATT) & far$ATT > 5
  expect_true(any(keep))
  for (f in c("ATT", "se", "lci", "uci", "groups", "periods")) far[[f]] <- far[[f]][keep]
  pf <- plot(far)
  xr <- vrange(pf, "x")
  expect_true(xr[1] <= 0 && xr[2] >= max(far$ATT))
  expect_identical(legend_pos(pf), "none")

  pz <- plot(fit_u, xlim = c(5, 30))
  xr <- vrange(pz, "x")
  expect_true(xr[1] <= 0 && xr[2] >= 30)
})

#-----------------------------------------------------------------------------
# T13: two-period fit
#-----------------------------------------------------------------------------
test_that("plot.ddd draws a single ATT row for a two-period fit", {
  p <- plot(fit_2)
  pts <- ggplot2::layer_data(p, which(geoms(p) == "GeomPoint"))
  expect_equal(nrow(pts), 1)
  expect_equal(pts$x, as.numeric(fit_2$ATT), tolerance = TOL)
  lr <- ggplot2::layer_data(p, which(geoms(p) == "GeomLinerange"))
  expect_equal(lr$xmin, as.numeric(fit_2$lci), tolerance = TOL)
  expect_equal(lr$xmax, as.numeric(fit_2$uci), tolerance = TOL)
  lims <- ggplot2::ggplot_build(p)$layout$panel_scales_y[[1]]$get_limits()
  expect_identical(lims, "ATT")
  expect_identical(pts$colour, COL_POST)
  expect_identical(legend_pos(p), "none")
  xr <- vrange(p, "x")
  expect_true(xr[1] <= 0 && xr[2] >= 0)
})

#-----------------------------------------------------------------------------
# T14: ggplot2 missing
#-----------------------------------------------------------------------------
test_that("a clear error is raised when ggplot2 is not available", {
  local_mocked_bindings(ddd_has_ggplot2 = function() FALSE)
  expect_error(plot(es), "needs the ggplot2 package", fixed = TRUE)
  expect_error(plot(fit_u), "needs the ggplot2 package", fixed = TRUE)
  expect_error(plot(es, ylim = "a"), "needs the ggplot2 package", fixed = TRUE)
})
