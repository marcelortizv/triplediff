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
#' @param dname The name of the column containing the treatment indicator dummy (1 if treated in the post-treatment period, 0 otherwise).
#' It is mutually exclusive with \code{gname}.
#' @param gname The name of the column containing the first period when a particular observation is treated. It is a positive number
#' for treated units and defines which group the unit belongs to. It takes value 0 for untreated units. If \code{gname} is specified,
#' we assume that the treatment is staggered. It is mutually exclusive with \code{dname}.
#' @param partition.name The name of the column containing the partition variable (e.g., the subgroup identifier). This is an indicator variable that is 1 for
#' the units targeted for treatment and 0 otherwise.
#' @param xformla The formula for the covariates to be included in the model. It should be of the form \code{~ x1 + x2}.
#' Default is \code{xformla = ~1} (no covariates).
#' @param data A data frame or data table containing the data.
#' @param control.group The control group to be used in the estimation. Default is \code{control.group = "nevertreated"} which sets as control group
#' the units that never participate in the treatment and does not change across groups or time periods. The alternative is
#' \code{control.group = "notyettreated"} which sets as control group the units that have not yet participated in the treatment. This includes never
#' treated units and units that will be treated in the future.
#' @param estMethod The estimation method to be used. Default is \code{"trad"} (traditional). It computes propensity score using logistic regression
#' and outcome regression using OLS. The alternative is \code{"dml"} (double machine learning). It allows the user to compute propensity score using a
#' machine learning algorithm and outcome regression using a different machine learning algorithm. We provide some wrappers for popular learners but
#' the user can also provide their own learner. See our vignette on How to construct a use-provided learner for more details (TODO).
#' @param learners A list of learners to be used in the estimation. It should be a list of two elements,
#' the first element being the learner for the propensity score and the second element being the learner for the outcome regression.
#' Default is \code{NULL}.
#' @param weightsname The name of the column containing the weights. Default is \code{NULL}.
#' @param boot Logical. If \code{TRUE}, the function computes the bootstrap standard errors. Default is \code{FALSE}.
#' @param boot.type The type of bootstrap to be used. Default is \code{"multiplier"}. The alternative is \code{"weighted"}.
#' @param nboot The number of bootstrap samples to be used. Default is \code{NULL}. If \code{boot = TRUE}, the default is \code{nboot = 999}.
#' @param inffunc Logical. If \code{TRUE}, the function returns the influence function. Default is \code{FALSE}.
#'
#' @return A list with the following elements:
#' \item{ATT}{The average treatment effect on the treated.}
#' \item{se}{The standard error of the ATT.}
#' \item{uci}{The upper confidence interval of the ATT.}
#' \item{lci}{The lower confidence interval of the ATT.}
#' \item{inf.func}{The estimate of the influence function.}
#'
#' @details
#'
#' TBA
#'
#' @examples
#' #----------------------------------------------------------
#' # Triple Diff with covariates and 2 time periods
#' #----------------------------------------------------------
#' seed = 123
#' num_ids = 100
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
#'     gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
#'     data = sim_data, control.group = NULL,
#'     estMethod = "trad", learners = NULL, weightsname = NULL, boot = FALSE,
#'     boot.type = "multiplier", nboot = NULL, inffunc = FALSE)
#'
#' #----------------------------------------------------------
#' # Triple Diff with multiple time periods
#' #----------------------------------------------------------
#'
#' \dontrun{
#' ddd(yname = "outcome", tname = "year", idname = "id", dname = NULL,
#'     gname = "group", partition.name = "partition", xformla = ~x1 + x2,
#'     data = sim_data, control.group = "notyettreated",
#'     estMethod = "trad", learners = NULL, weightsname = NULL, boot = FALSE,
#'     boot.type = "multiplier", nboot = NULL, inffunc = FALSE)
#' }
#'
#'
#' @export

ddd <- function(yname, tname, idname, dname, gname, partition.name, xformla,
                data, control.group = NULL,
                estMethod = "trad", learners = NULL, weightsname = NULL,
                boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE) {


  #------------------------------------------
  # Running initial arguments validation
  #------------------------------------------

  # Check what's the setting: 2 time periods or multiple time periods
  if (is.null(dname) && is.null(gname)) {
    stop("Either dname or gname should be provided")
  } else if (!is.null(dname) && !is.null(gname)) {
    stop("Only one of dname or gname should be provided")
  } else if (!is.null(dname)){
    multiple.periods <- FALSE
  } else if(!is.null(gname)){
    multiple.periods <- TRUE
  }

  # Flag for estMethod
  if ((estMethod!="dml") && (estMethod!="trad")) {
    warning("estMethod = ", estMethod, " is not supported. Using 'trad'.")
    estMethod <- "trad"
  }

  #------------------------------------------
  # Run preprocess and validation checks
  #------------------------------------------

  if ((multiple.periods) && (estMethod=="trad")) {
    # dp <- run_preprocess_multPeriods(yname = yname,
    #                                  tname = tname,
    #                                  idname = idname,
    #                                  dname = NULL,
    #                                  gname = gname,
    #                                  partition.name = partition.name,
    #                                  xformla = xformla,
    #                                  data = data,
    #                                  control.group = control.group,
    #                                  estMethod = "trad",
    #                                  learners = NULL,
    #                                  weightsname = weightsname,
    #                                  boot = boot,
    #                                  boot.type = boot.type,
    #                                  nboot = nboot,
    #                                  inffunc = inffunc)
    stop("Triple Diff with multiple time periods is not yet supported")
  } else if ((multiple.periods) && (estMethod=="dml")) {
    # dp <- run_preprocess_multPeriods(yname = yname,
    #                                  tname = tname,
    #                                  idname = idname,
    #                                  dname = NULL,
    #                                  gname = gname,
    #                                  partition.name = partition.name,
    #                                  xformla = xformla,
    #                                  data = data,
    #                                  control.group = control.group,
    #                                  estMethod = "dml",
    #                                  learners = learners,
    #                                  weightsname = weightsname,
    #                                  boot = boot,
    #                                  boot.type = boot.type,
    #                                  nboot = nboot,
    #                                  inffunc = inffunc)
    stop("Triple Diff with multiple time periods and DML is not yet supported")
  } else if ((!multiple.periods) && (estMethod=="trad")) {
    dp <- run_preprocess_2Periods(yname = yname,
                                 tname = tname,
                                 idname = idname,
                                 dname = dname,
                                 gname = NULL,
                                 partition.name = partition.name,
                                 xformla = xformla,
                                 data = data,
                                 control.group = NULL,
                                 estMethod = "trad",
                                 learners = NULL,
                                 weightsname = weightsname,
                                 boot = boot,
                                 boot.type = boot.type,
                                 nboot = nboot,
                                 inffunc = inffunc)
  } else if ((!multiple.periods) && (estMethod=="dml")) {
    # dp <- run_preprocess_2Periods(yname = yname,
    #                              tname = tname,
    #                              idname = idname,
    #                              dname = dname,
    #                              gname = NULL,
    #                              partition.name = partition.name,
    #                              xformla = xformla,
    #                              data = data,
    #                              control.group = NULL,
    #                              estMethod = "dml",
    #                              learners = learners,
    #                              weightsname = weightsname,
    #                              boot = boot,
    #                              boot.type = boot.type,
    #                              nboot = nboot,
    #                              inffunc = inffunc)
    stop("Triple Diff with 2 periods and DML is not yet supported")
  }

  #------------------------------------------
  # Run the estimation
  #------------------------------------------

  # 2 time periods case: trad or dml
  if (multiple.periods == FALSE){
    if (estMethod == "trad"){
      att_dr <- att_dr(dp)
    }
  } else {
      #TODO: IMPLEMENT att_dml PROCEDURE AND ADJUST PARAMETERS
      # att_dml <- att_dml(dp)
  }#RUN DML ESTIMATION

  # multiple time periods case: trad or dml
  if (multiple.periods == TRUE){
    if (estMethod == "trad"){
      #TODO: IMPLEMENT att_gt_dr PROCEDURE AND ADJUST PARAMETERS
      # att_gt_dr <- att_gt_dr(dp)
    }
  } else {
      #TODO: IMPLEMENT att_gt_dml PROCEDURE AND ADJUST PARAMETERS
      # att_gt_dml <- att_gt_dml(dp)
  }#RUN DML ESTIMATION

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
    partition.name = args$partition.name,
    control.group = args$control.group,
    estMethod = estMethod,
    multiple.periods = multiple.periods,
    learners = args$learners,
    boot = args$boot,
    boot.type = args$boot.type
  )

  if (multiple.periods == FALSE){
    if (estMethod == "trad"){
        ret <- list(
          ATT = att_dr$ATT,
          se = att_dr$se,
          lci = att_dr$lci,
          uci = att_dr$uci,
          nboot = att_dr$nboot,
          att.inf.func = att_dr$inf.func,
          call.params = call.params,
          argu = argu
        )
  } #else {
   #      ret <- list(
   #        ATT = att_dml$ATT,
   #        se = att_dml$se,
   #        lci = att_dml$lci,
   #        uci = att_dml$uci,
   #        att.inf.func = att_dml$inf.func,
   #        call.params = call.params,
   #        argu = argu
   #      )
   # }
  }# RETURNING LIST FOR 2 PERIODS CASE

  # multiple time periods case: trad or dml
  # if (multiple.periods == TRUE){
  #   if (estMethod == "trad"){
  #     ret <- list(
  #       ATT = att_gt_dr$ATT,
  #       se = att_gt_dr$se,
  #       lci = att_gt_dr$lci,
  #       uci = att_gt_dr$uci,
  #       att.inf.func = att_gt_dr$inf.func,
  #       call.params = call.params,
  #       argu = argu
  #     )
  #   } else {
  #     ret <- list(
  #       ATT = att_gt_dml$ATT,
  #       se = att_gt_dml$se,
  #       lci = att_gt_dml$lci,
  #       uci = att_gt_dml$uci,
  #       att.inf.func = att_gt_dml$inf.func,
  #       call.params = call.params,
  #       argu = argu
  #     )
  #   }
  # }#RETURNING LIST FOR MULTIPLE PERIODS CASE

  # define a class
  class(ret) <- "ddd"

  return(ret)
  }


