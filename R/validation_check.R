# This process run all the error check before enter to the preprocessing

#' @importFrom Rcpp sourceCpp

validate_args_2Periods <- function(args, dta){
  # get args
  yname <- args$yname
  tname <- args$tname
  idname <- args$idname
  dname <- args$dname
  partition.name <- args$partition.name
  xformla <- args$xformla
  estMethod <- args$estMethod
  learners <- args$learners
  n_folds <- args$n_folds
  weightsname <- args$weightsname
  boot <- args$boot
  boot.type <- args$boot.type
  nboot <- args$nboot
  inffunc <- args$inffunc


  # Flag for yname
  if (!is.element(yname, base::colnames(dta))) {
    stop("yname = ",yname,  " could not be found in the data provided.")
  }

  # Flag for tname
  if (!is.element(tname, base::colnames(dta))) {
    stop("tname = ",tname,  " could not be found in the data provided.")
  }

  # check if times periods are numeric
  if (!all(sapply(dta[, ..tname], is.numeric))) {
    stop("tname = ",tname,  " is not numeric. Please convert it")
  }

  # Check if there is only 2 time periods
  tlist <- unique(dta[[tname]])[base::order(unique(dta[[tname]]))]
  if (length(tlist) != 2) {
    stop("The type of ddd specified only allow for two time periods (pre and post). Change type of ddd for multiple time periods")
  }

  # Flag for dname
  if ( !is.element(dname, base::colnames(dta))) {
    stop("dname = ",dname,  " could not be found in the data provided.")
  }

  # Check if there is only 2 groups
  dlist <- unique(dta[[dname]])[base::order(unique(dta[[dname]]))]
  if (length(dlist) != 2) {
    stop("The type of ddd specified only allow for two groups (treated and untreated). Change type of ddd for multiple groups")
  }

  # Flag for partition.name
  if ( !is.element(partition.name, base::colnames(dta))) {
    stop("partition.name = ",partition.name,  " could not be found in the data provided.")
  }


  # check if partition.name is numeric
  if (!all(sapply(dta[, ..partition.name], is.numeric))) {
    stop("partition.name = ",partition.name,  " is not numeric. Please convert it")
  }

  # Check if partition values are binary
  plist <- unique(dta[[partition.name]])[base::order(unique(dta[[partition.name]]))]
  if (length(plist) != 2) {
    stop("partition.name =", partition.name, " must have only two values (0 and 1). Please check partition.name")
  }

  # Check if idname is in the data
  if ( !is.element(idname, base::colnames(dta))) {
    stop("idname = ",idname,  " could not be found in the data provided.")
  }

  #  check if idname is numeric
  if (!all(sapply(dta[, ..idname], is.numeric))) {
    stop("data[, idname] must be numeric. Please convert it.")
  }

  # Check if any combination of idname and tname is duplicated
  n_id_year = anyDuplicated(dta[, .(idname, tname)])
  # If any combination is duplicated, stop execution and throw an error
  if (n_id_year > 0) {
    stop("The value of idname must be unique (by tname)")
  }

  # Check if partition is unique by idname
  checkPartitionUniqueness(dta, idname, partition.name)

  # Check if dname is unique by idname
  checkTreatmentUniqueness(dta, idname, dname)

  # Flag for weightsname
  if(!is.null(weightsname)){
    if ( !is.element(weightsname, base::colnames(dta))) {
      stop("weightsname = ",weightsname,  " could not be found in the data provided.")
    }
  }

  # FLAGS FOR DML ESTIMATION
  if (estMethod=="dml") {
    if (is.null(learners)) {
      stop("learners should be provided when estMethod = 'dml'")
    }

    if (is.null(n_folds)) {
      stop("n_folds should be provided when estMethod = 'dml'")
    }

    # check if there 2 learners provided in learners list
    if (length(learners) != 2) {
      stop("learners must be a list of two learners: ml_pa and ml_md")
    }

    # check if learner is a classifier type
    if (!inherits(learners$ml_pa, "LearnerClassif")) {
      stop("The learner must be a classification learner.")
    }

    # check if learner  is a regression type
    if (!inherits(learners$ml_md, "LearnerRegr")) {
      stop("The learner must be a regression learner.")
    }

  }

}

# TODO: ADD ARGUMENTS VALIDATION FOR MULTIPLE PERIODS SETTING
# validate_args_multPeriods <- function(args, dta){
#   # TODO: Identify what code is repetitive and generate validate_general_args() to use in both functions
#   # add code
# }
