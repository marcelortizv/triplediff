# Main function for aggregations steps in triplediff
NULL
#' Aggregate Group-Time Average Treatment Effects.
#' Aggregation procedures are based on Callaway, Brantly and Pedro H.C. Sant'Anna.
#' "Difference-in-Differences with Multiple Time Periods." Journal of Econometrics,
#' Vol. 225, No. 2, pp. 200-230, 2021.
#' @description
#' \code{aggddd} is a function that take group-time average treatment effects
#'  and aggregate them into a smaller number of parameters.  This method is only valid
#'  when there are multiple time periods and staggered treatment adoption. There are
#'  several possible aggregations including "simple", "eventstudy", "group",
#'  and "calendar."  Default is \code{"simple"}.
#'
#' @param ddd_obj a ddd object (i.e., the results of the [ddd()] function)
#' @param type Which type of aggregated treatment effect parameter to compute.
#'   \code{"simple"} just computes a weighted average of all
#'   group-time average treatment effects with weights proportional to group
#'   size.  This is the default option;
#'   \code{"eventstudy"} computes average effects across
#'   different lengths of exposure to the treatment and is similar to an
#'   "event study". Here the overall effect averages the effect of the
#'   treatment across all positive lengths of exposure;
#'   \code{"group"} computes average treatment effects across different groups; here
#'   the overall effect averages the effect across different groups; and
#'   \code{"calendar"} computes average treatment effects across different
#'   time periods; here the overall effect averages the effect across each
#'   time period.
#' @param balance_e If set (and if one computes event study), it balances
#'  the sample with respect to event time.  For example, if `balance.e=2`,
#'  `aggddd` will drop groups that are not exposed to treatment for
#'  at least three periods. (the initial period when `e=0` as well as the
#'  next two periods when `e=1` and the `e=2`).  This ensures that
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
#'
#' @param cband Boolean for whether or not to compute a uniform confidence
#'  band that covers all of the group-time average treatment effects
#'  with fixed probability `0.95`.  In order to compute uniform confidence
#'  bands, `boot` must also be set to `TRUE`.  The default is
#'  the value set in the ddd object
#'
#' @return A object (list) of class [`aggddd`] that holds the results from the
#'  aggregation step.  The object contains the following elements:
#'
#' @examples
#' #----------------------------------------------------------
#' # Triple Diff with multiple time periods
#' #----------------------------------------------------------
#' \dontrun{
#' data <- gen_dgp_mult_periods(size = 1000, tperiods = 4, dgp_type = 1)
#'
#' out <- ddd(yname = "Y", tname = "period", idname = "id", dname = NULL,
#'             gname = "G", partition_name = "L", xformla = ~X,
#'             data = data, control_group = "nevertreated", base_period = "varying",
#'             est_method = "trad")
#' ** Simple aggregation **
#' aggddd(out, type = "simple")
#'
#' ** Event study aggregation **
#' aggddd(out, type = "eventstudy")
#'
#' ** Group aggregation **
#' aggddd(out, type = "group")
#'
#' ** Calendar aggregation **
#' aggddd(out, type = "calendar")
#' }
#'
#' @export


aggddd <- function(ddd_obj,
                   type = "simple",
                   balance_e = NULL,
                   min_e = -Inf,
                   max_e = Inf,
                   na.rm = FALSE,
                   boot = FALSE,
                   nboot = NULL,
                   cband = NULL) {

 call <- match.call()
 # Record all arguments used in the function
 args <- mget(names(formals()), sys.frame(sys.nframe()))
 argu <- list(
   boot = args$boot,
   nboot = args$nboot,
   cband = args$cband
 )

 aggte <- compute_aggregation(ddd_obj = ddd_obj,
                            type = type,
                            balance_e = balance_e,
                            min_e = min_e,
                            max_e = max_e,
                            na.rm = na.rm,
                            boot = boot,
                            nboot = nboot,
                            cband = cband
                            )

 ret <- list(aggte = aggte,
             call.params = call,
             argu = argu)

 # define a class
 class(ret) <- "aggddd"

 return(ret)

}
