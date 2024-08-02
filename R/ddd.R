# Main function for triplediff
NULL
#' Doubly Robust DDD estimators for the ATT and beyond
#'
#' @description
#' \code{ddd} is the main function for computing the Doubly Robust DDD estimators for the ATT, with panel data.
#' It can be used with covariates and/or under multiple time periods. At its core, \code{triplediff} employs
#' the doubly robust estimator for the ATT, which is a combination of the propensity score weighting and the outcome regression.
#' Furthermore, this package supports the application of machine learning methods for the estimation of the nuisance functions.
#'
#' @param yname The name of the outcome variable.
#' @param tname The name of the column containing the time periods.
#' @param idname The name of the column containing the unit id.
#' @param dname Valid for 2 time periods only. The name of the column containing the treatment indicator dummy (1 if treated in the post-treatment period, 0 otherwise).
#' It is mutually exclusive with \code{gname}. If \code{dname} is provided, \code{gname} is ignored.
#' @param gname Valid for multiple periods only. The name of the column containing the first period when a particular observation is treated. It is a positive number
#' for treated units and defines which group the unit belongs to. It takes value 0 or Inf for untreated units. If \code{gname} is specified,
#' we assume that the treatment is staggered. It is mutually exclusive with \code{dname}. If \code{gname} is provided, \code{dname} is ignored.
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
#' and outcome regression using OLS. The alternative is \code{"dml"} (double machine learning). It allows the user to compute propensity score using a
#' machine learning algorithm and outcome regression using a different machine learning algorithm. We provide some wrappers for popular learners but
#' the user can also provide their own learner. See our vignette on How to construct a use-provided learner for more details (# TODO).
#' @param learners A list of learners to be used in the estimation. It should be a list of two elements,
#' the first element being the learner for the propensity score and the second element being the learner for the outcome regression.
#' Default is \code{NULL}, then OLS and MLE Logit is used to estimate nuisances parameters. If \code{est_method = "dml"}, user have to specify \code{learners}.
#' @param n_folds The number of folds to be used in the cross-fitting. Default is \code{NULL}. If \code{est_method = "dml"}, user have to specify \code{n_folds}.
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
#' @param skip_data_checks Logical. If \code{TRUE}, the function skips the data checks and go straight to estimation. Default is \code{FALSE}.
#'
#' @return A `ddd` object with the following basic elements:
#' \item{ATT}{The average treatment effect on the treated.}
#' \item{se}{The standard error of the ATT.}
#' \item{uci}{The upper confidence interval of the ATT.}
#' \item{lci}{The lower confidence interval of the ATT.}
#' \item{inf_func}{The estimate of the influence function.}
#' \item{...}{Other elements that are specific to the estimation method.}
#'
#'
#' @examples
#' #----------------------------------------------------------
#' # Triple Diff with covariates and 2 time periods
#' #----------------------------------------------------------
#' seed = 123
#' num_ids = 500
#' time = 2
#' initial.year = 2019
#' treatment.year = 2020
#'
#' sim_data = generate_test_panel(seed = seed,
#'                               num_ids = num_ids,
#'                               time = time,
#'                               initial.year = initial.year,
#'                               treatment.year = treatment.year)
#'
#' ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
#'     gname = NULL, pname = "partition", xformla = ~x1 + x2,
#'     data = sim_data, control_group = NULL,
#'     est_method = "dr")
#'
#' #----------------------------------------------------------
#' # DML Triple Diff with covariates and 2 time periods
#' #----------------------------------------------------------
#' library(mlr3)
#' library(mlr3learners)
#'
#' learner_rf <- lrn("classif.ranger", predict_type = "prob", num.trees = 100,
#'                  min.node.size = 1, importance = 'impurity')
#'
#' learner_regr <- lrn("regr.xgboost")
#'
#' learners <- list(ml_pa = learner_rf, ml_md = learner_regr)
#'
#' ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
#'     gname = NULL, pname = "partition", xformla = ~x1 + x2,
#'     data = sim_data, control_group = NULL,
#'     est_method = "dml", learners = learners, n_folds = 3)
#'
#' #----------------------------------------------------------
#' # Triple Diff with multiple time periods
#' #----------------------------------------------------------
#' data <- gen_dgp_mult_periods(size = 1000, dgp_type = 1)[["data"]]
#'
#' ddd(yname = "y", tname = "time", idname = "id", dname = NULL,
#'      gname = "state", pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
#'      data = data, control_group = "nevertreated", base_period = "varying",
#'      est_method = "dr")
#'
#' #----------------------------------------------------------
#' # DML Triple Diff with multiple time periods
#' #----------------------------------------------------------
#' # TBA
#' @export

ddd <- function(yname,
                tname,
                idname,
                dname,
                gname,
                pname,
                xformla,
                data,
                control_group = NULL,
                base_period = NULL,
                est_method = "dr",
                learners = NULL,
                n_folds = NULL,
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

  # Check what's the setting: 2 time periods or multiple time periods
  if (is.null(dname) && is.null(gname)) {
    stop("Either dname or gname should be provided")
  } else if (!is.null(dname)){
    multiple_periods <- FALSE
    gname <- NULL
  } else if (!is.null(gname)){
    multiple_periods <- TRUE
    dname <- NULL
  }

  # Flag for est_method
  if ((est_method != "dml") && (est_method != "dr")) {
    warning("est_method = ", est_method, " is not supported. Using 'dr'.")
    est_method <- "dr"
  }

  # Flag for control group in multiple periods
  if (multiple_periods && is.null(control_group)) {
      warning("control_group should be provided for multiple time periods. Using 'notyettreated'")
      control_group <- "notyettreated"
  }

  #------------------------------------------
  # Run preprocess and validation checks
  #------------------------------------------

  if (skip_data_checks){
    if (!multiple_periods){
      dp <- run_nopreprocess_2periods(yname = yname,
                                      tname = tname,
                                      idname = idname,
                                      dname = dname,
                                      gname = NULL,
                                      pname = pname,
                                      xformla = xformla,
                                      data = data,
                                      control_group = NULL,
                                      est_method = "dr",
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
      stop("Triple Diff with multiple time periods is not yet supported")
    }

  } else {
    if ((multiple_periods) && (est_method =="dr")) {
      dp <- run_preprocess_multPeriods(yname = yname,
                                       tname = tname,
                                       idname = idname,
                                       dname = NULL,
                                       gname = gname,
                                       pname = pname,
                                       xformla = xformla,
                                       data = data,
                                       control_group = control_group,
                                       base_period = base_period,
                                       est_method = "dr",
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

    } else if ((multiple_periods) && (est_method == "dml")) {
      # dp <- run_preprocess_multPeriods(yname = yname,
      #                                  tname = tname,
      #                                  idname = idname,
      #                                  dname = NULL,
      #                                  gname = gname,
      #                                  pname = pname,
      #                                  xformla = xformla,
      #                                  data = data,
      #                                  control_group = control_group,
      #                                  base_period = base_period,
      #                                  est_method = "dml",
      #                                  learners = learners,
      #                                  n_folds = n_folds,
      #                                  weightsname = weightsname,
      #                                  boot = boot,
      #                                  nboot = nboot,
      #                                  inffunc = inffunc)
      stop("Triple Diff with multiple time periods and DML is not yet supported")
    } else if ((!multiple_periods) && (est_method == "dr")) {
      dp <- run_preprocess_2Periods(yname = yname,
                                    tname = tname,
                                    idname = idname,
                                    dname = dname,
                                    gname = NULL,
                                    pname = pname,
                                    xformla = xformla,
                                    data = data,
                                    control_group = NULL,
                                    est_method = "dr",
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
    } else if ((!multiple_periods) && (est_method == "dml")) {
      dp <- run_preprocess_2Periods(yname = yname,
                                    tname = tname,
                                    idname = idname,
                                    dname = dname,
                                    gname = NULL,
                                    pname = pname,
                                    xformla = xformla,
                                    data = data,
                                    control_group = NULL,
                                    est_method = "dml",
                                    learners = learners,
                                    n_folds = n_folds,
                                    weightsname = weightsname,
                                    boot = boot,
                                    nboot = nboot,
                                    cluster = cluster,
                                    cband = cband,
                                    alpha = alpha,
                                    use_parallel = use_parallel,
                                    cores = cores,
                                    inffunc = inffunc)
    }
  }

  #------------------------------------------
  # Run the estimation
  #------------------------------------------

  # multiple time periods case: dr or dml
  if (multiple_periods){
    if (est_method == "dr"){
      att_gt_dr <- att_gt(dp)
    } else {
      stop("DML for multiple time periods is not yet supported")
      # TODO: IMPLEMENT att_gt_dml PROCEDURE AND ADJUST PARAMETERS
      # att_gt_dml <- att_gt_dml(dp)
    }#RUN DML for multiple time periods
  } else {
    if (est_method == "dr"){
      att_dr <- att_dr(dp)
    } else {
      att_dml <- att_dml(dp)
    }#RUN DML for 2 time periods
  }#RUN DR for 2 time periods

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
    learners = args$learners,
    n_folds = args$n_folds,
    cband = args$cband,
    cluster = args$cluster,
    boot = dp$boot, # getting from dp because it could change in the pre process
    alpha = dp$alpha, # getting from dp because it could change in the pre process
    nboot = dp$nboot, # getting from dp because it could change in the pre process
    cores = dp$cores # getting from dp because it could change in the pre process
  )

  if (!multiple_periods){
    if (est_method == "dr"){
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
  } else {
        ret <- list(
          ATT = att_dml$ATT,
          se = att_dml$se,
          lci = att_dml$lci,
          uci = att_dml$uci,
          subgroup_counts = att_dml$subgroup_counts,
          call.params = call.params,
          argu = argu
        )
   }
  }# RETURNING LIST FOR 2 PERIODS CASE

  # multiple time periods case: dr or dml
  if (multiple_periods){
    if (est_method == "dr"){
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


