# Main function for aggregations steps in triplediff
NULL
#' Aggregate Group-Time Average Treatment Effects in Staggered Triple-Differences Designs
#' Aggregation procedures
#' @description
#' \code{agg_ddd} is a function that take group-time average treatment effects
#'  and aggregate them into a smaller number of summary parameters in staggered triple differences designs.
#'  There are several possible aggregations including "simple", "eventstudy", "group",
#'  and "calendar."  Default is \code{"eventstudy"}.
#'
#' @param ddd_obj a ddd object (i.e., the results of the [ddd()] function)
#' @param type Which type of aggregated treatment effect parameter to compute.
#'   \code{"simple"} just computes a weighted average of all
#'   group-time average treatment effects with weights proportional to group
#'   size.
#'   \code{"eventstudy"} computes average effects across
#'   different lengths of exposure to the treatment (event times). Here the overall effect averages the effect of the
#'   treatment across the positive lengths of exposure. This is the default option;
#'   \code{"group"} computes average treatment effects across different groups/cohorts; here
#'   the overall effect averages the effect across different groups using group size as weights;
#'   \code{"calendar"} computes average treatment effects across different
#'   time periods, with weights proportional to the group size; here the overall effect averages the effect across each
#'   time period.
#' @param balance_e If set (and if one computes event study), it balances
#'  the sample with respect to event time.  For example, if `balance_e=2`,
#'  `agg_ddd` will drop groups that are not exposed to treatment for
#'  at least three periods, the initial period `e=0` as well as the
#'  next two periods, `e=1` and `e=2`.  This ensures that
#'  the composition of groups does not change when event time changes.
#' @param min_e For event studies, this is the smallest event time to compute
#'  dynamic effects for.  By default, `min_e = -Inf` so that effects at
#'  all lengths of exposure are computed.
#' @param max_e For event studies, this is the largest event time to compute
#'  dynamic effects for.  By default, `max_e = Inf` so that effects at
#'  all lengths of exposure are computed.
#' @param na.rm Logical value if we are to remove missing Values from analyses. Defaults is FALSE.
#' @param boot Boolean for whether or not to compute standard errors using
#'  the multiplier bootstrap.  If standard errors are clustered, then one
#'  must set `boot=TRUE`. Default is value set in the ddd object.  If `boot = FALSE`, then analytical
#'  standard errors are reported.
#' @param nboot The number of bootstrap iterations to use.  The default is the value set in the ddd object,
#'  and this is only applicable if `boot=TRUE`.
#' @param alpha The level of confidence for the confidence intervals.  The default is 0.05. Otherwise, it will
#' use the value set in the ddd object.
#' @param cband Boolean for whether or not to compute a uniform confidence
#'  band that covers all of the group-time average treatment effects
#'  with fixed probability `0.95`.  In order to compute uniform confidence
#'  bands, `boot` must also be set to `TRUE`.  The default is
#'  the value set in the ddd object
#'
#' @return A object (list) of class [`agg_ddd`] that holds the results from the
#'  aggregation step.  The object contains the following elements:
#'
#' @examples
#' #----------------------------------------------------------
#' # Triple Diff with multiple time periods
#' #----------------------------------------------------------
#'
#' data <- gen_dgp_mult_periods(size = 1000, dgp_type = 1)[["data"]]
#'
#' out <- ddd(yname = "y", tname = "time", idname = "id", dname = NULL,
#'             gname = "state", pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
#'             data = data, control_group = "nevertreated", base_period = "varying",
#'             est_method = "trad")
#' # Simple aggregation
#' agg_ddd(out, type = "simple", alpha = 0.10)
#'
#' # Event study aggregation
#' agg_ddd(out, type = "eventstudy", alpha = 0.10)
#'
#' # Group aggregation
#' agg_ddd(out, type = "group", alpha = 0.10)
#'
#' # Calendar aggregation
#' agg_ddd(out, type = "calendar", alpha = 0.10)
#'
#'
#' @export


agg_ddd <- function(ddd_obj,
                   type = "eventstudy",
                   balance_e = NULL,
                   min_e = -Inf,
                   max_e = Inf,
                   na.rm = FALSE,
                   boot = NULL,
                   nboot = NULL,
                   cband = NULL,
                   alpha = 0.05) {

 call <- match.call()
 # Record all arguments used in the function
 # args <- mget(names(formals()), sys.frame(sys.nframe()))
 # argu <- list(
 #   boot = args$boot,
 #   nboot = args$nboot,
 #   cband = args$cband,
 #   alpha = args$alpha
 # )

 aggte_ddd <- compute_aggregation(ddd_obj = ddd_obj,
                                  type = type,
                                  balance_e = balance_e,
                                  min_e = min_e,
                                  max_e = max_e,
                                  na.rm = na.rm,
                                  boot = boot,
                                  nboot = nboot,
                                  cband = cband,
                                  alpha = alpha
                                )


 ret <- list(aggte_ddd = aggte_ddd,
             call.params = call)#,
             #argu = argu)

 # define a class
 class(ret) <- "agg_ddd"

 return(ret)

}
