# Plot methods for ddd and agg_ddd objects
NULL

# Columns of the internal plotting data frames, mapped inside ggplot2::aes().
utils::globalVariables(c("xval", "est", "lo", "hi", "phase", "label"))

#-----------------------------------------------------------------------------
# Constants: legend labels and the panelty palette
#-----------------------------------------------------------------------------
PLOT_LAB_PRE  <- "Pre-treatment"
PLOT_LAB_POST <- "Post-treatment"

PLOT_COL_POST    <- "#E8743B" # post-treatment points and intervals
PLOT_COL_PRE     <- "#414F5D" # pre-treatment points and intervals
PLOT_COL_REF     <- "#414F5D" # zero line, e = -0.5 line, axis lines and ticks
PLOT_COL_GRID    <- "#E4E8ED" # value-axis grid lines
PLOT_COL_SUBTEXT <- "#555555" # subtitle text

PLOT_REF_LINEWIDTH  <- 0.4  # reference lines, axis lines, ticks and grid
PLOT_BAR_LINEWIDTH  <- 0.6  # error bars and interval lines
PLOT_BAR_WIDTH      <- 0.15 # error-bar cap width (agg plots), in x units
PLOT_POINT_SIZE     <- 2.4
PLOT_BASE_SIZE      <- 11
PLOT_SUBTITLE_REL   <- 0.85
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

# Shared theme; value_axis is the axis that carries the estimates.
ddd_plot_theme <- function(value_axis) {
  th <- ggplot2::theme_minimal(base_size = PLOT_BASE_SIZE) +
    ggplot2::theme(
      panel.grid.minor    = ggplot2::element_blank(),
      axis.line           = ggplot2::element_line(colour = PLOT_COL_REF,
                                                  linewidth = PLOT_REF_LINEWIDTH),
      axis.ticks          = ggplot2::element_line(colour = PLOT_COL_REF,
                                                  linewidth = PLOT_REF_LINEWIDTH),
      plot.title.position = "plot",
      plot.subtitle       = ggplot2::element_text(colour = PLOT_COL_SUBTEXT,
                                                  size = ggplot2::rel(PLOT_SUBTITLE_REL)),
      legend.position     = "bottom",
      legend.title        = ggplot2::element_blank(),
      plot.background     = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background    = ggplot2::element_rect(fill = "white", colour = NA)
    )
  grid_line <- ggplot2::element_line(colour = PLOT_COL_GRID, linewidth = PLOT_REF_LINEWIDTH)
  if (value_axis == "y") {
    th + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                        panel.grid.major.y = grid_line)
  } else {
    th + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                        panel.grid.major.x = grid_line)
  }
}

#-----------------------------------------------------------------------------
# plot.agg_ddd
#-----------------------------------------------------------------------------

#' Plot aggregated DDD effects
#'
#' @description
#' Draws the path of an event-study, group or calendar aggregation computed by
#' [agg_ddd()]. Each estimate is a point with its confidence interval as an
#' error bar. Pre-treatment estimates are dark slate and post-treatment
#' estimates are orange. A dashed horizontal line marks zero.
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
#' estimates at `e >= 0` are post-treatment. A dotted vertical line at
#' `e = -0.5` separates the two, and a legend names them; both appear only
#' when pre- and post-treatment estimates are present. Group and calendar
#' aggregations average post-treatment cells only, so all their estimates are
#' post-treatment and no legend is drawn.
#'
#' The vertical axis always includes zero, with or without `ylim`. Estimates
#' that are `NA` are not drawn, with a warning.
#'
#' ggplot2 must be installed; it is a suggested package. The result is a
#' regular `ggplot` object, so it can be changed with `+`, for example
#' `plot(es) + ggplot2::labs(title = "Event study")`.
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

  p <- ggplot2::ggplot(df, ggplot2::aes(x = xval, y = est, colour = phase)) +
    ggplot2::geom_hline(yintercept = 0, colour = PLOT_COL_REF,
                        linewidth = PLOT_REF_LINEWIDTH, linetype = "dashed")
  if (is_es && two_phases) {
    p <- p + ggplot2::geom_vline(xintercept = PLOT_ES_SEPARATOR, colour = PLOT_COL_REF,
                                 linewidth = PLOT_REF_LINEWIDTH, linetype = "dotted")
  }
  p <- p +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = PLOT_BAR_WIDTH,
                           linewidth = PLOT_BAR_LINEWIDTH, na.rm = TRUE) +
    ggplot2::geom_point(size = PLOT_POINT_SIZE) +
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
#' `"ATT(g,t)"`. Each estimate is a point with its confidence interval as a
#' horizontal line. Pre-treatment cells (`t < g`) are dark slate and
#' post-treatment cells (`t >= g`) are orange. A dashed vertical line marks
#' zero. For a two-period fit, the plot has a single row, `"ATT"`.
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
#' ggplot2 must be installed; it is a suggested package. The result is a
#' regular `ggplot` object, so it can be changed with `+`, for example
#' `plot(out) + ggplot2::labs(title = "ATT(g,t)")`. To plot an event study,
#' group or calendar path, aggregate first with [agg_ddd()] and plot the result.
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
    ggplot2::geom_vline(xintercept = 0, colour = PLOT_COL_REF,
                        linewidth = PLOT_REF_LINEWIDTH, linetype = "dashed") +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lo, xmax = hi), orientation = "y",
                            linewidth = PLOT_BAR_LINEWIDTH, na.rm = TRUE) +
    ggplot2::geom_point(size = PLOT_POINT_SIZE) +
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
