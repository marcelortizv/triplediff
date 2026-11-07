# Plot methods for ddd and agg_ddd objects
NULL

# Columns of the internal plotting data frames, mapped inside ggplot2::aes().
# `colour` is the scaled colour column that after_scale reads in plot_bar_aes().
utils::globalVariables(c("xval", "est", "lo", "hi", "phase", "label", "colour"))

#-----------------------------------------------------------------------------
# Constants: legend labels, colours and sizes
#-----------------------------------------------------------------------------
PLOT_LAB_PRE  <- "Pre-treatment"
PLOT_LAB_POST <- "Post-treatment"

PLOT_COL_POST  <- "#154780" # post-treatment points and bars
PLOT_COL_PRE   <- "#B28C1C" # pre-treatment points and bars
PLOT_COL_ZERO  <- "#415D56" # zero line
PLOT_COL_ONSET <- "#93A69F" # event-study line at e = -0.5
PLOT_COL_MUTED <- "#5C746C" # axis lines, axis text, axis titles and subtitle
PLOT_COL_GRID  <- "#E3E9E5" # vertical grid lines of plot.ddd()
PLOT_COL_INK   <- "#1E2B27" # legend text and row labels of plot.ddd()

PLOT_BAR_ALPHA        <- 0.6  # share of the colour in the bar tint, over white
PLOT_BAR_LINEWIDTH    <- 1.1  # interval bars, thick lines with no caps
PLOT_POINT_SIZE       <- 2.4
PLOT_POINT_SHAPE      <- 17   # triangle, in both phases
PLOT_REF_LINEWIDTH    <- 0.4  # zero line and event-study line
PLOT_THIN_LINEWIDTH   <- 0.3  # axis lines and grid lines
PLOT_BASE_SIZE        <- 11
PLOT_AXIS_TEXT_SIZE   <- 9.5
PLOT_AXIS_TITLE_SIZE  <- 10
PLOT_SUBTITLE_SIZE    <- 9.5
PLOT_LEGEND_TEXT_SIZE <- 10.5
PLOT_TEXT_GAP         <- 4    # pt, below the subtitle and between a legend key and its label
PLOT_LEGEND_KEY_WIDTH <- 10   # pt
PLOT_LEGEND_KEY_GAP   <- 36   # pt between two legend entries
PLOT_LEGEND_BOX_GAP   <- 9    # pt between the legend and the panel

PLOT_ES_SEPARATOR   <- -0.5 # event-study line between pre and post periods
PLOT_X_EXPAND_ADD   <- 0.5  # additive x expansion of agg plots
PLOT_MAX_ALL_BREAKS <- 20   # up to this many periods, every period gets a tick
PLOT_PRETTY_N       <- 10   # target number of ticks beyond that
PLOT_SIGNIF_LEVEL   <- 4    # significant digits of the confidence level text

#-----------------------------------------------------------------------------
# ggplot2 guard
#-----------------------------------------------------------------------------
ddd_has_ggplot2 <- function() requireNamespace("ggplot2", quietly = TRUE)

check_ggplot2 <- function() {
  if (!ddd_has_ggplot2()) {
    stop("plot() for triplediff objects needs the ggplot2 package. ",
         "Install it with install.packages(\"ggplot2\").", call. = FALSE)
  }
  invisible(TRUE)
}

#-----------------------------------------------------------------------------
# Argument checks
#-----------------------------------------------------------------------------
check_plot_label <- function(value, name) {
  if (!is.null(value) &&
      !(is.character(value) && length(value) == 1 && !is.na(value))) {
    stop("`", name, "` must be NULL or a single character string.", call. = FALSE)
  }
  invisible(TRUE)
}

check_plot_lim <- function(value, name) {
  if (!is.null(value) &&
      !(is.numeric(value) && length(value) == 2 && all(is.finite(value)))) {
    stop("`", name, "` must be NULL or a numeric vector of length 2 with finite values.",
         call. = FALSE)
  }
  invisible(TRUE)
}

check_plot_subtitle <- function(value) {
  if (!is.null(value) &&
      !(is.character(value) && length(value) == 1 && !is.na(value))) {
    stop("`subtitle` must be \"auto\", NULL, or a single character string.", call. = FALSE)
  }
  invisible(TRUE)
}

#-----------------------------------------------------------------------------
# Helpers
#-----------------------------------------------------------------------------

# Subtitle that names the interval drawn in the plot.
interval_subtitle <- function(uniform, alpha, boot, nboot, cluster) {
  pct   <- paste0(signif(100 * (1 - alpha), PLOT_SIGNIF_LEVEL), "%")
  draws <- if (isTRUE(boot) && !is.null(nboot)) paste0(", ", nboot, " draws") else ""
  cl    <- if (length(cluster) > 0) {
    paste0(", clustered by ", paste(cluster, collapse = ", "))
  } else {
    ""
  }
  if (uniform) {
    paste0("Uniform ", pct, " confidence band (multiplier bootstrap", draws, cl, ")")
  } else if (isTRUE(boot)) {
    paste0("Pointwise ", pct, " confidence intervals (multiplier bootstrap SEs",
           draws, cl, ")")
  } else {
    paste0("Pointwise ", pct, " confidence intervals (analytical SEs", cl, ")")
  }
}

# Ticks on the period axis: every drawn period, or integer pretty breaks when
# there are many periods.
period_breaks <- function(v) {
  u <- sort(unique(v))
  if (length(u) <= PLOT_MAX_ALL_BREAKS) return(u)
  b <- pretty(range(u), n = PLOT_PRETTY_N)
  b[b == round(b)]
}

# Drop rows with a missing estimate, with a warning; stop if none is left.
drop_na_estimates <- function(df) {
  keep <- is.finite(df$est)
  if (!any(keep)) {
    stop("Every estimate is NA; there is nothing to plot.", call. = FALSE)
  }
  if (any(!keep)) {
    warning(sprintf("%d estimate(s) are NA and are not drawn.", sum(!keep)),
            call. = FALSE)
  }
  df[keep, , drop = FALSE]
}

# Colour scale of the two phases (pre- and post-treatment).
plot_phase_scale <- function() {
  ggplot2::scale_colour_manual(
    values = stats::setNames(c(PLOT_COL_PRE, PLOT_COL_POST),
                             c(PLOT_LAB_PRE, PLOT_LAB_POST)),
    drop = TRUE
  )
}

# Opaque tint of a colour: the share `alpha` of the colour over white. The bars
# are opaque, so they hide the zero line and the grid lines behind them.
plot_tint <- function(col, alpha = PLOT_BAR_ALPHA) {
  m <- grDevices::col2rgb(col) / 255
  grDevices::rgb(t(alpha * m + (1 - alpha)))
}

# Mapping of the interval bars: the colour of the phase, tinted after the scale.
# The scale still sees `phase`, so a colour scale added by the user also recolours
# the bars. ggplot2 3.5 reads after_scale only in a call to the bare name stage(),
# so bind that name here. A call written ggplot2::stage() is not read in that version.
plot_bar_aes <- function(...) {
  stage <- ggplot2::stage
  ggplot2::aes(..., colour = stage(phase, after_scale = plot_tint(colour)))
}

# Shared theme. value_axis is the axis that carries the estimates. Event study,
# group and calendar plots (value_axis = "y") have axis lines and no grid. The
# coefficient plot (value_axis = "x") has vertical grid lines and no axis lines.
ddd_plot_theme <- function(value_axis) {
  th <- ggplot2::theme_minimal(base_size = PLOT_BASE_SIZE) +
    ggplot2::theme(
      plot.background      = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background     = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major     = ggplot2::element_blank(),
      panel.grid.minor     = ggplot2::element_blank(),
      panel.border         = ggplot2::element_blank(),
      axis.ticks           = ggplot2::element_blank(),
      axis.text            = ggplot2::element_text(colour = PLOT_COL_MUTED,
                                                   size = PLOT_AXIS_TEXT_SIZE),
      axis.title           = ggplot2::element_text(colour = PLOT_COL_MUTED,
                                                   size = PLOT_AXIS_TITLE_SIZE),
      plot.title.position  = "plot",
      plot.subtitle        = ggplot2::element_text(colour = PLOT_COL_MUTED,
                                                   size = PLOT_SUBTITLE_SIZE,
                                                   margin = ggplot2::margin(b = PLOT_TEXT_GAP)),
      legend.position      = "top",
      legend.direction     = "horizontal",
      legend.justification = "center",
      legend.title         = ggplot2::element_blank(),
      legend.text          = ggplot2::element_text(colour = PLOT_COL_INK,
                                                   size = PLOT_LEGEND_TEXT_SIZE,
                                                   margin = ggplot2::margin(l = PLOT_TEXT_GAP)),
      legend.key           = ggplot2::element_blank(),
      legend.key.width     = ggplot2::unit(PLOT_LEGEND_KEY_WIDTH, "pt"),
      legend.key.spacing.x = ggplot2::unit(PLOT_LEGEND_KEY_GAP, "pt"),
      legend.margin        = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing   = ggplot2::unit(PLOT_LEGEND_BOX_GAP, "pt")
    )
  if (value_axis == "y") {
    th + ggplot2::theme(axis.line = ggplot2::element_line(colour = PLOT_COL_MUTED,
                                                          linewidth = PLOT_THIN_LINEWIDTH,
                                                          lineend = "square"))
  } else {
    th + ggplot2::theme(panel.grid.major.x = ggplot2::element_line(colour = PLOT_COL_GRID,
                                                                   linewidth = PLOT_THIN_LINEWIDTH),
                        axis.text.y        = ggplot2::element_text(colour = PLOT_COL_INK,
                                                                   size = PLOT_AXIS_TEXT_SIZE))
  }
}

#-----------------------------------------------------------------------------
# plot.agg_ddd
#-----------------------------------------------------------------------------

#' Plot aggregated DDD effects
#'
#' @description
#' Draws the path of an event-study, group or calendar aggregation computed by
#' [agg_ddd()]. Each estimate is a triangle with its confidence interval as a
#' thick bar in a lighter tint of the same color. Pre-treatment estimates are
#' gold and post-treatment estimates are dark blue. A dashed horizontal line
#' marks zero.
#'
#' @param x An `agg_ddd` object returned by [agg_ddd()] with
#'   `type = "eventstudy"`, `"group"` or `"calendar"`.
#' @param ylim `NULL` (default) or a numeric vector of length 2 with finite
#'   values. It zooms the vertical (estimate) axis without removing data. The
#'   range shown is always widened to include zero.
#' @param xlab `NULL` (default) or a single character string. `NULL` uses the
#'   default label of the aggregation type: `"Event time"`,
#'   `"Group enabling period"` or `"Calendar period"`.
#' @param ylab `NULL` or a single character string for the vertical axis label.
#'   Default is `"Estimate"`. `NULL` removes the label.
#' @param title `NULL` (default, no title) or a single character string.
#' @param subtitle `"auto"` (default), `NULL`, or a single character string.
#'   `"auto"` describes the intervals drawn, for example
#'   `"Pointwise 95% confidence intervals (analytical SEs)"`. `NULL` removes the
#'   subtitle.
#' @param ... Not used.
#'
#' @details
#' The intervals are the ones stored in the object: the uniform confidence band
#' when `cband = TRUE` was used, pointwise intervals otherwise. The plot does not
#' recompute them with another critical value.
#'
#' In an event study, estimates at event times `e < 0` are pre-treatment and
#' estimates at `e >= 0` are post-treatment. A dashed vertical line at
#' `e = -0.5` separates the two. The line appears only
#' when pre- and post-treatment estimates are present. Group and calendar
#' aggregations average post-treatment cells only, so all their estimates are
#' post-treatment. The legend is at the top and appears only when both kinds of
#' estimates are present.
#'
#' The vertical axis always includes zero, with or without `ylim`. Estimates
#' that are `NA` are not drawn, with a warning.
#'
#' ggplot2 version 3.5.0 or later must be installed. It is a suggested package.
#' The result is a regular `ggplot` object, so it can be changed with `+`, for
#' example `plot(es) + ggplot2::labs(title = "Event study")`. A color scale
#' added with `+` changes the colors of the triangles and of the bars.
#'
#' @return A `ggplot` object. It is not printed by the function.
#'
#' @seealso [agg_ddd()], [ddd()]
#'
#' @examples
#' \dontshow{.old_dt_threads <- data.table::setDTthreads(2)}
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1234)
#'   data <- gen_dgp_mult_periods(size = 500, dgp_type = 1)[["data"]]
#'   out <- ddd(yname = "y", tname = "time", idname = "id", gname = "state",
#'              pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
#'              data = data, control_group = "nevertreated",
#'              base_period = "universal", est_method = "dr")
#'   print(plot(agg_ddd(out, type = "eventstudy")))
#'   print(plot(agg_ddd(out, type = "group")))
#' }
#' \dontshow{data.table::setDTthreads(.old_dt_threads)}
#'
#' @export
plot.agg_ddd <- function(x, ylim = NULL, xlab = NULL, ylab = "Estimate", title = NULL,
                         subtitle = "auto", ...) {
  check_ggplot2()
  check_plot_lim(ylim, "ylim")
  check_plot_label(xlab, "xlab")
  check_plot_label(ylab, "ylab")
  check_plot_label(title, "title")
  check_plot_subtitle(subtitle)

  a <- x$aggte_ddd
  type <- a$type
  if (identical(type, "simple")) {
    stop("A \"simple\" aggregation has one overall ATT and no path to plot. ",
         "Use summary(), or aggregate with type = \"eventstudy\", \"group\" or \"calendar\".",
         call. = FALSE)
  }
  if (!(is.character(type) && length(type) == 1 &&
        type %in% c("eventstudy", "group", "calendar"))) {
    stop("plot() supports aggregations of type \"eventstudy\", \"group\" and \"calendar\".",
         call. = FALSE)
  }

  # Inference: draw exactly the interval stored by the aggregation
  alpha <- a$argu$alpha
  z <- stats::qnorm(1 - alpha / 2)
  crit <- a$crit.val.egt
  uniform <- isTRUE(a$argu$cband) && isTRUE(is.finite(crit)) &&
    isTRUE(crit > z + sqrt(.Machine$double.eps))

  df <- data.frame(xval = as.numeric(a$egt),
                   est  = as.numeric(a$att.egt),
                   se   = as.numeric(a$se.egt))
  df$lo <- df$est - crit * df$se
  df$hi <- df$est + crit * df$se
  df <- drop_na_estimates(df)

  is_es <- type == "eventstudy"
  pre <- if (is_es) df$xval < 0 else rep(FALSE, nrow(df))
  df$phase <- factor(ifelse(pre, PLOT_LAB_PRE, PLOT_LAB_POST),
                     levels = c(PLOT_LAB_PRE, PLOT_LAB_POST))
  two_phases <- length(unique(df$phase)) == 2

  if (is.null(xlab)) {
    xlab <- switch(type,
                   eventstudy = "Event time",
                   group      = "Group enabling period",
                   calendar   = "Calendar period")
  }
  if (identical(subtitle, "auto")) {
    subtitle <- interval_subtitle(uniform, alpha, a$argu$boot, a$argu$nboot, a$argu$cluster)
  }

  b <- period_breaks(df$xval)

  # The event-study line comes first, so the zero line is drawn over it where they cross.
  p <- ggplot2::ggplot(df, ggplot2::aes(x = xval, y = est, colour = phase))
  if (is_es && two_phases) {
    p <- p + ggplot2::geom_vline(xintercept = PLOT_ES_SEPARATOR, colour = PLOT_COL_ONSET,
                                 linewidth = PLOT_REF_LINEWIDTH, linetype = "dashed")
  }
  p <- p +
    ggplot2::geom_hline(yintercept = 0, colour = PLOT_COL_ZERO,
                        linewidth = PLOT_REF_LINEWIDTH, linetype = "dashed") +
    ggplot2::geom_linerange(plot_bar_aes(ymin = lo, ymax = hi),
                            linewidth = PLOT_BAR_LINEWIDTH, lineend = "butt",
                            na.rm = TRUE, show.legend = FALSE) +
    ggplot2::geom_point(shape = PLOT_POINT_SHAPE, size = PLOT_POINT_SIZE) +
    plot_phase_scale() +
    ggplot2::scale_x_continuous(breaks = b,
                                labels = format(b, trim = TRUE, scientific = FALSE),
                                expand = ggplot2::expansion(add = PLOT_X_EXPAND_ADD))
  if (!is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(ylim = range(c(ylim, 0)))
  }
  p <- p +
    ggplot2::labs(x = xlab, y = ylab, title = title, subtitle = subtitle, colour = NULL) +
    ddd_plot_theme("y")
  if (!two_phases) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  p
}

#-----------------------------------------------------------------------------
# plot.ddd
#-----------------------------------------------------------------------------

#' Plot group-time DDD effects
#'
#' @description
#' Draws a horizontal coefficient plot of the group-time average treatment
#' effects ATT(g,t) estimated by [ddd()], with one row per ATT(g,t) labelled
#' `"ATT(g,t)"`. Each estimate is a triangle with its confidence interval as a
#' thick horizontal bar in a lighter tint of the same color. Pre-treatment
#' cells (`t < g`) are gold and post-treatment cells (`t >= g`) are dark blue.
#' A dashed vertical line marks zero. For a two-period fit, the plot has a
#' single row, `"ATT"`.
#'
#' @param x A `ddd` object returned by [ddd()].
#' @param xlim `NULL` (default) or a numeric vector of length 2 with finite
#'   values. It zooms the horizontal (estimate) axis without removing data. The
#'   range shown is always widened to include zero.
#' @param xlab `NULL` or a single character string for the horizontal axis
#'   label. Default is `"Estimate"`. `NULL` removes the label.
#' @param ylab `NULL` (default, no label) or a single character string for the
#'   vertical axis label.
#' @param title `NULL` (default, no title) or a single character string.
#' @param subtitle `"auto"` (default), `NULL`, or a single character string.
#'   `"auto"` describes the intervals drawn, for example
#'   `"Pointwise 95% confidence intervals (analytical SEs)"`. `NULL` removes the
#'   subtitle.
#' @param ... Not used.
#'
#' @details
#' The intervals are the ones stored in the object: the uniform confidence band
#' when `cband = TRUE` was used, pointwise intervals otherwise.
#'
#' The rows keep the order stored in the object (groups in ascending order,
#' then periods in ascending order), with the first cell at the top. A cell is
#' pre-treatment when its period is before its group's enabling period
#' (`t < g`) and post-treatment otherwise; the legend appears only when both
#' kinds are present. With a universal base period, the reference cell is
#' drawn as a point at zero with no interval.
#'
#' The horizontal axis always includes zero, with or without `xlim`, so the
#' zero line is always in view. Estimates that are `NA` are not drawn, with a
#' warning.
#'
#' ggplot2 version 3.5.0 or later must be installed. It is a suggested package.
#' The result is a regular `ggplot` object, so it can be changed with `+`, for
#' example `plot(out) + ggplot2::labs(title = "ATT(g,t)")`. A color scale added
#' with `+` changes the colors of the triangles and of the bars. To plot an
#' event study, group or calendar path, aggregate first with [agg_ddd()] and
#' plot the result.
#'
#' @return A `ggplot` object. It is not printed by the function.
#'
#' @seealso [agg_ddd()], [ddd()]
#'
#' @examples
#' \dontshow{.old_dt_threads <- data.table::setDTthreads(2)}
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1234)
#'   data <- gen_dgp_mult_periods(size = 500, dgp_type = 1)[["data"]]
#'   out <- ddd(yname = "y", tname = "time", idname = "id", gname = "state",
#'              pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
#'              data = data, control_group = "nevertreated",
#'              base_period = "universal", est_method = "dr")
#'   plot(out)
#' }
#' \dontshow{data.table::setDTthreads(.old_dt_threads)}
#'
#' @export
plot.ddd <- function(x, xlim = NULL, xlab = "Estimate", ylab = NULL, title = NULL,
                     subtitle = "auto", ...) {
  check_ggplot2()
  check_plot_lim(xlim, "xlim")
  check_plot_label(xlab, "xlab")
  check_plot_label(ylab, "ylab")
  check_plot_label(title, "title")
  check_plot_subtitle(subtitle)

  # One row per stored ATT(g,t), or a single row for a two-period fit
  if (isTRUE(x$argu$multiple_periods)) {
    g <- as.numeric(x$groups)
    t <- as.numeric(x$periods)
    df <- data.frame(est = as.numeric(x$ATT),
                     lo  = as.numeric(x$lci),
                     hi  = as.numeric(x$uci))
    df$label <- paste0("ATT(", format(g, trim = TRUE, scientific = FALSE), ",",
                       format(t, trim = TRUE, scientific = FALSE), ")")
    pre <- t < g
  } else {
    df <- data.frame(est = as.numeric(x$ATT),
                     lo  = as.numeric(x$lci),
                     hi  = as.numeric(x$uci),
                     label = "ATT")
    pre <- FALSE
  }
  df$phase <- factor(ifelse(pre, PLOT_LAB_PRE, PLOT_LAB_POST),
                     levels = c(PLOT_LAB_PRE, PLOT_LAB_POST))
  df <- drop_na_estimates(df)

  # First stored cell at the top: on a discrete y axis the first level is at
  # the bottom, so reverse the levels.
  df$label <- factor(df$label, levels = rev(unique(df$label)))

  uniform <- isTRUE(x$argu$cband) && isTRUE(x$argu$boot)
  if (identical(subtitle, "auto")) {
    subtitle <- interval_subtitle(uniform, x$argu$alpha, x$argu$boot, x$argu$nboot,
                                  x$argu$cluster)
  }
  two_phases <- length(unique(df$phase)) == 2

  p <- ggplot2::ggplot(df, ggplot2::aes(x = est, y = label, colour = phase)) +
    ggplot2::geom_vline(xintercept = 0, colour = PLOT_COL_ZERO,
                        linewidth = PLOT_REF_LINEWIDTH, linetype = "dashed") +
    ggplot2::geom_linerange(plot_bar_aes(xmin = lo, xmax = hi), orientation = "y",
                            linewidth = PLOT_BAR_LINEWIDTH, lineend = "butt",
                            na.rm = TRUE, show.legend = FALSE) +
    ggplot2::geom_point(shape = PLOT_POINT_SHAPE, size = PLOT_POINT_SIZE) +
    plot_phase_scale()
  if (!is.null(xlim)) {
    p <- p + ggplot2::coord_cartesian(xlim = range(c(xlim, 0)))
  }
  p <- p +
    ggplot2::labs(x = xlab, y = ylab, title = title, subtitle = subtitle, colour = NULL) +
    ddd_plot_theme("x")
  if (!two_phases) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  p
}
