# Main function for triplediff
NULL
#' Doubly Robust DDD estimators for the group-time average treatment effects.
#'
#' @description
#' \code{ddd} is the main function for computing the Doubly Robust DDD estimators for the ATT, with balanced panel data.
#' It can be used with covariates and/or under multiple time periods. At its core, \code{triplediff} employs
#' the doubly robust estimator for the ATT, which is a combination of the propensity score weighting and the outcome regression.
#' Furthermore, this package supports the application of machine learning methods for the estimation of the nuisance parameters.
#' @param yname The name of the outcome variable.
#' @param tname The name of the column containing the time periods.
#' @param idname The name of the column containing the unit id.
#' @param gname The name of the column containing the first period when a particular observation is treated. It is a positive number
#' for treated units and defines which group the unit belongs to. It takes value 0 or Inf for untreated units.
#' @param pname The name of the column containing the partition variable (e.g., the subgroup identifier). This is an indicator variable that is 1 for
#' the units eligible for treatment and 0 otherwise.
#' @param xformla The formula for the covariates to be included in the model. It should be of the form \code{~ x1 + x2}.
#' Default is \code{xformla = ~1} (no covariates).
#' @param data A data frame or data table containing the data.
#' @param control_group Valid for multiple periods only. The control group to be used in the estimation. Default is \code{control_group = "notyettreated"} which sets as control group the units that have not yet participated in the treatment.
#' The alternative is \code{control_group = "nevertreated"} which sets as control group the units that never participate in the treatment and does not change across groups or time periods.
#' @param base_period Valid for multiple periods. Choose between a "varying" or "universal" base period. Both yield the same post-treatment ATT(g,t) estimates.
#' Varying base period: Computes pseudo-ATT in pre-treatment periods by comparing outcome changes for a group to its comparison group from t-1 to t, repeatedly changing t.
#' Universal base period: Fixes the base period to (g-1), reporting average changes from t to (g-1) for a group relative to its comparison group, similar to event study regressions.
#' Varying base period reports ATT(g,t) right before treatment. Universal base period normalizes the estimate before treatment to be 0, adding one extra estimate in an earlier period.
#' @param est_method The estimation method to be used. Default is \code{"dr"} (doubly robust). It computes propensity score using logistic regression
#' and outcome regression using OLS. The alternative are \code{c("reg", "ipw")}.
#' @param panel Logical. If \code{TRUE} (default), the data is treated as panel data where each unit is observed in all time periods.
#' If \code{FALSE}, the data is treated as repeated cross-sections (RCS) where each observation may represent a different unit.
#' For RCS data, \code{idname} can be omitted or set to \code{NULL}, and the function will automatically create unique IDs for each observation.
#' @param allow_unbalanced_panel Logical. If \code{TRUE}, allows for unbalanced panel data where units may not be observed in all time periods.
#' Default is \code{FALSE}. Note: This parameter requires \code{panel = TRUE} and a valid \code{idname}.
#' @param weightsname The name of the column containing the weights. Default is \code{NULL}. As part of data processing, weights are enforced to be normalized
#' and have mean 1 across all observations.
#' @param boot Logical. If \code{TRUE}, the function computes standard errors using the multiplier bootstrap. Default is \code{FALSE}.
#' @param nboot The number of bootstrap samples to be used. Default is \code{NULL}. If \code{boot = TRUE}, the default is \code{nboot = 999}.
#' @param cluster The name of the variable to be used for clustering. The maximum number of cluster variables is 1. Default is \code{NULL}.
#' If \code{boot = TRUE}, the function computes the bootstrap standard errors clustering at the unit level setting as cluster variable the one in \code{idname}.
#' @param cband Logical. If \code{TRUE}, the function computes a uniform confidence band that covers all of the average treatment effects
#' with fixed probability \code{1-alpha}.  In order to compute uniform confidence bands, \code{boot} must also be set to `TRUE`.  The default is `FALSE`.
#' @param alpha The level of significance for the confidence intervals.  Default is \code{0.05}.
#' @param use_parallel Logical. If \code{TRUE}, the function runs in parallel processing. Valid only when \code{boot = TRUE}. Default is \code{FALSE}.
#' @param cores The number of cores to be used in the parallel processing. Default is \code{cores = 1}.
#' @param inffunc Logical. If \code{TRUE}, the function returns the influence function. Default is \code{FALSE}.
#' @param skip_data_checks Logical. If \code{TRUE}, the function skips data validation checks and proceeds directly to estimation.
#' This can improve performance when you are confident the data is correctly formatted. Default is \code{FALSE}.
#' Use with caution as skipping checks may lead to unexpected errors if data is malformed.
#'
#' @return A `ddd` object with the following basic elements:
#' \item{ATT}{The average treatment effect on the treated.}
#' \item{se}{The standard error of the ATT.}
#' \item{uci}{The upper confidence interval of the ATT.}
#' \item{lci}{The lower confidence interval of the ATT.}
#' \item{inf_func}{The estimate of the influence function.}
#'
#'
#' @examples
#' #----------------------------------------------------------
#' # Triple Diff with covariates and 2 time periods
#' #----------------------------------------------------------
#' set.seed(1234) # Set seed for reproducibility
#' # Simulate data for a two-periods DDD setup
#' df <- gen_dgp_2periods(size = 5000, dgp_type = 1)$data
#'
#' head(df)
#'
#' att_22 <- ddd(yname = "y", tname = "time", idname = "id", gname = "state",
#'               pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
#'              data = df, control_group = "nevertreated", est_method = "dr")
#'
#' summary(att_22)
#'
#' # Performing clustered standard errors with mutiplier bootstrap
#'
#' att_cluster <-  ddd(yname = "y", tname = "time", idname = "id", gname = "state",
#' pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
#' data = df, control_group = "nevertreated",
#' base_period = "universal", est_method = "dr", cluster = "cluster")
#'
#' summary(att_cluster)
#'
#' #----------------------------------------------------------
#' # Triple Diff with multiple time periods
#' #----------------------------------------------------------
#' data <- gen_dgp_mult_periods(size = 1000, dgp_type = 1)[["data"]]
#'
#' ddd(yname = "y", tname = "time", idname = "id",
#'      gname = "state", pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
#'      data = data, control_group = "nevertreated", base_period = "varying",
#'      est_method = "dr")
#'
#' @export
ddd <- function(yname,
                tname,
                idname = NULL,
                gname,
                pname,
                xformla,
                data,
                control_group = NULL,
                base_period = NULL,
                est_method = "dr",
                panel = TRUE,
                allow_unbalanced_panel = FALSE,
                weightsname = NULL,
                boot = FALSE,
                nboot = NULL,
                cluster = NULL,
                cband = FALSE,
                alpha = 0.05,
                use_parallel = FALSE,
                cores = 1,
                inffunc = FALSE,
                skip_data_checks = FALSE) {


  #------------------------------------------
  # Running initial arguments validation
  #------------------------------------------

  # Check if 'dta' is a data.table
  if (!"data.table" %in% class(data)) {
    # converting data to data.table
    dta <- data.table::as.data.table(data)
  } else {
    dta <- data
  }

  # Check what's the setting: 2 time periods or multiple time periods
  if (is.null(gname)) {
    stop("gname should be provided")
  } else {
    unique_gname_values <- dta[, uniqueN(get(gname))]
    unique_tname_values <- dta[, uniqueN(get(tname))]
    if (max(unique_gname_values,unique_tname_values) == 2) {
      multiple_periods <- FALSE
    } else if (max(unique_gname_values,unique_tname_values)  > 2) {
      multiple_periods <- TRUE
    } else {
      stop("Invalid gname. Please check your arguments.")
    }
  }

  #------------------------------------------
  # Validate idname parameter
  #------------------------------------------

  # If idname is not provided, assume RCS data
  if (is.null(idname)) {
    # Check for incompatible parameter combinations
    if (panel) {
      stop("idname is required when panel = TRUE. For repeated cross-section data, set panel = FALSE.")
    }
    if (allow_unbalanced_panel) {
      stop("idname is required when allow_unbalanced_panel = TRUE. For repeated cross-section data, set panel = FALSE and leave idname as NULL.")
    }
    # Set a placeholder - preprocessing will create the actual ID
    idname <- ".rcs_id"
  } else {
    # User provided idname - check if column exists
    if (!idname %in% names(dta)) {
      stop(paste0("Column '", idname, "' not found in data."))
    }
  }

  # Flag for est_method
  if (!(est_method %in% c("reg", "ipw", "dr"))) {
    warning("est_method = ", est_method, " is invalid or not supported. Using 'dr'.")
    est_method <- "dr"
  }

  # Flag for control group in multiple periods
  if (multiple_periods && is.null(control_group)) {
      warning("control_group should be provided for multiple time periods. Using 'nevertreated'")
      control_group <- "nevertreated"
  }

  # Flag for base period in multiple periods
  if (multiple_periods && is.null(base_period)) {
      warning("base_period should be provided for multiple time periods. Using 'varying'")
      base_period <- "varying"
  }

  #------------------------------------------
  # Run preprocess and validation checks
  #------------------------------------------

  if (skip_data_checks){
    if (!multiple_periods){
      dp <- run_nopreprocess_2periods(yname = yname,
                                      tname = tname,
                                      idname = idname,
                                      gname = gname,
                                      pname = pname,
                                      xformla = xformla,
                                      dta = dta,
                                      control_group = NULL,
                                      est_method = est_method,
                                      panel = panel,
                                      allow_unbalanced_panel = allow_unbalanced_panel,
                                      learners = NULL,
                                      n_folds = NULL,
                                      weightsname = weightsname,
                                      boot = boot,
                                      nboot = nboot,
                                      cluster = cluster,
                                      cband = cband,
                                      alpha = alpha,
                                      use_parallel = use_parallel,
                                      cores = cores,
                                      inffunc = inffunc)
    } else {
      stop("Triple Diff with multiple time periods is not yet supported.")
    }

  } else {
    if ((multiple_periods) && (est_method %in% c("dr", "reg", "ipw"))) {
      dp <- run_preprocess_multPeriods(yname = yname,
                                       tname = tname,
                                       idname = idname,
                                       gname = gname,
                                       pname = pname,
                                       xformla = xformla,
                                       dta = dta,
                                       control_group = control_group,
                                       base_period = base_period,
                                       est_method = est_method,
                                       panel = panel,
                                       allow_unbalanced_panel = allow_unbalanced_panel,
                                       learners = NULL,
                                       n_folds = NULL,
                                       weightsname = weightsname,
                                       boot = boot,
                                       nboot = nboot,
                                       cluster = cluster,
                                       cband = cband,
                                       alpha = alpha,
                                       use_parallel = use_parallel,
                                       cores = cores,
                                       inffunc = inffunc)

    } #else if ((multiple_periods) && (est_method == "dml")) {
      # dp <- run_preprocess_multPeriods(yname = yname,
      #                                  tname = tname,
      #                                  idname = idname,
      #                                  gname = gname,
      #                                  pname = pname,
      #                                  xformla = xformla,
      #                                  dta = dta,
      #                                  control_group = control_group,
      #                                  base_period = base_period,
      #                                  est_method = "dml",
      #                                  learners = learners,
      #                                  n_folds = n_folds,
      #                                  weightsname = weightsname,
      #                                  boot = boot,
      #                                  nboot = nboot,
      #                                  inffunc = inffunc)
      #stop("Triple Diff with multiple time periods and DML is not yet supported")}
    else if ((!multiple_periods) && (est_method %in% c("dr", "reg", "ipw"))) {
      dp <- run_preprocess_2Periods(yname = yname,
                                    tname = tname,
                                    idname = idname,
                                    gname = gname,
                                    pname = pname,
                                    xformla = xformla,
                                    dta = dta,
                                    control_group = NULL,
                                    est_method = est_method,
                                    panel = panel,
                                    allow_unbalanced_panel = allow_unbalanced_panel,
                                    learners = NULL,
                                    n_folds = NULL,
                                    weightsname = weightsname,
                                    boot = boot,
                                    nboot = nboot,
                                    cluster = cluster,
                                    cband = cband,
                                    alpha = alpha,
                                    use_parallel = use_parallel,
                                    cores = cores,
                                    inffunc = inffunc)
    } # else if ((!multiple_periods) && (est_method == "dml")) {
      # dp <- run_preprocess_2Periods(yname = yname,
      #                               tname = tname,
      #                               idname = idname,
      #                               gname = gname,
      #                               pname = pname,
      #                               xformla = xformla,
      #                               dta = dta,
      #                               control_group = NULL,
      #                               est_method = "dml",
      #                               learners = learners,
      #                               n_folds = n_folds,
      #                               weightsname = weightsname,
      #                               boot = boot,
      #                               nboot = nboot,
      #                               cluster = cluster,
      #                               cband = cband,
      #                               alpha = alpha,
      #                               use_parallel = use_parallel,
      #                               cores = cores,
      #                               inffunc = inffunc)}
  }

  #------------------------------------------
  # Run the estimation
  #------------------------------------------

  # multiple time periods case: dr or dml
  if (multiple_periods){
    # RUN DR for multiple periods
    if (est_method %in% c("dr", "reg", "ipw")){
      att_gt_dr <- att_gt(dp)
    } #else {
      # # RUN DML for multiple time periods
      # stop("DML for multiple time periods is not yet supported")
      # # TODO: IMPLEMENT att_gt_dml PROCEDURE AND ADJUST PARAMETERS
      # # att_gt_dml <- att_gt_dml(dp)
    #}
  } else {
    if (est_method %in% c("dml")){
        # RUN DML for 2 time periods
        # att_dml <- att_dml(dp)
      stop("DML estimation method is not yet supported.")
    } else {
      # Use true_repeated_cross_sections flag to determine estimation method
      # Extract with fallback for backward compatibility
      true_rcs <- if (!is.null(dp$true_repeated_cross_sections)) {
        dp$true_repeated_cross_sections
      } else {
        !dp$panel  # Fallback: if flag missing, assume RCS only when panel=FALSE
      }

      if (dp$panel && !true_rcs){
        # RUN DR for 2 time periods (balanced panel only)
        att_dr <- att_dr(dp)
      } else {
        # RUN DR for RCS 2 time periods (true RCS or unbalanced panel)
        att_dr <- att_dr_rc(dp)
      }
    }
  }


  #------------------------------------------
  # Return the results
  #------------------------------------------

  # record the call
  call.params <- match.call()
  # Record all arguments used in the function
  # Capture all arguments except 'data'
  arg_names <- setdiff(names(formals()), "data")
  args <- mget(arg_names, sys.frame(sys.nframe()))
  argu <- list(
    yname = args$yname,
    pname = args$pname,
    control_group = args$control_group,
    est_method = est_method,
    multiple_periods = multiple_periods,
    # learners = args$learners,
    # n_folds = args$n_folds,
    panel = args$panel, # getting from args because it could change in the pre process
    allow_unbalanced_panel = args$allow_unbalanced_panel, # getting from args because it could change in the pre process
    cband = dp$cband, # getting from dp because it could change in the pre process
    cluster = args$cluster,
    boot = dp$boot, # getting from dp because it could change in the pre process
    alpha = dp$alpha, # getting from dp because it could change in the pre process
    nboot = dp$nboot, # getting from dp because it could change in the pre process
    cores = dp$cores # getting from dp because it could change in the pre process
  )

  if (!multiple_periods){
    if (est_method %in% c("dr", "reg", "ipw")){
        ret <- list(
          ATT = att_dr$ATT,
          se = att_dr$se,
          lci = att_dr$lci,
          uci = att_dr$uci,
          nboot = att_dr$nboot, # this is not included under argu because it could change in the estimation process
          bT = att_dr$bT,
          att_inf_func = att_dr$inf_func,
          subgroup_counts = att_dr$subgroup_counts,
          call.params = call.params,
          argu = argu
        )
  } #else {
   #      ret <- list(
   #        ATT = att_dml$ATT,
   #        se = att_dml$se,
   #        lci = att_dml$lci,
   #        uci = att_dml$uci,
   #        subgroup_counts = att_dml$subgroup_counts,
   #        call.params = call.params,
   #        argu = argu
   #      )
   # }
  }# RETURNING LIST FOR 2 PERIODS CASE with DML

  # multiple time periods case: dr or dml
  if (multiple_periods){
    if (est_method %in% c("dr", "reg", "ipw")){
      ret <- list(
        ATT = att_gt_dr$ATT, # this is a vector of attgt
        se = att_gt_dr$se, # this is a vector of std. error for each attgt
        lci = att_gt_dr$lci, # this a vector of lower bound of CI for each attgt
        uci = att_gt_dr$uci, # this is a vector of upper bound of CI for each attgt
        groups = att_gt_dr$groups,
        periods = att_gt_dr$periods,
        tlist = att_gt_dr$tlist,
        glist = att_gt_dr$glist,
        cohort_size = att_gt_dr$cohort_size,
        n = att_gt_dr$n,
        bT = att_gt_dr$bT,
        inf_func_mat = att_gt_dr$inf_func_mat,
        first_period_dta = att_gt_dr$first_period_dta,
        call.params = call.params,
        argu = argu
      )
    }
    # } else {
    #   ret <- list(
    #     ATT = att_gt_dml$ATT,
    #     se = att_gt_dml$se,
    #     lci = att_gt_dml$lci,
    #     uci = att_gt_dml$uci,
    #     att.inf.func = att_gt_dml$inf.func,
    #     call.params = call.params,
    #     argu = argu
    #   )
    # }
  }#RETURNING LIST FOR MULTIPLE PERIODS CASE

  # define a class
  class(ret) <- "ddd"

  return(ret)
  }


