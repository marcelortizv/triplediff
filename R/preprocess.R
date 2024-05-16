#' @import data.table
#' @import stats
#' @import BMisc
NULL
#--------------------------------------------------
# Function the pre-process the data to use on ddd estimator

# Preprocess function for 2 periods case
run_preprocess_2Periods <- function(yname,
                                   tname,
                                   idname,
                                   dname,
                                   gname = NULL,
                                   partition_name,
                                   xformla = ~1,
                                   data,
                                   control_group = NULL,
                                   est_method = "trad",
                                   learners = NULL,
                                   n_folds = NULL,
                                   weightsname = NULL,
                                   boot = FALSE,
                                   boot_type = "multiplier",
                                   nboot = NULL,
                                   inffunc = FALSE){


  #-------------------------------------
  # Error checking
  #-------------------------------------

  # Flag for boot_type
  if (boot){
    if (boot_type!="multiplier") {
      warning("boot_type = ",boot_type,  " is not supported. Using 'multiplier'.")
      boot_type <- "multiplier"
    }
  }

  # Flag for est_method
  if (est_method!="trad" && est_method!="dml") {
    warning("est_method = ",est_method,  " is not supported. Using 'trad'.")
    est_method <- "trad"
  }

  # Check if 'dta' is a data.table
  if (!"data.table" %in% class(data)) {
    # converting data to data.table
    dta <- data.table::as.data.table(data)
  }
  dta<-data


  # Capture all arguments except 'data'
  arg_names <- setdiff(names(formals()), "data")
  args <- mget(arg_names, sys.frame(sys.nframe()))
  # Run argument checks
  validate_args_2Periods(args, dta)

  # set weights
  weights <- base::ifelse(is.null(weightsname), rep(1,nrow(dta)), dta[,weightsname])
  dta$weights <- weights

  # Check if weights are unique by idname
  checkWeightsUniqueness(dta, idname)


  # Flag for xformla
  if (is.null(xformla)) {
    xformla <- ~1
  }

  tryCatch(
    expr = {
      # attempt to convert xfomrla argument to a formula
      xformla <- as.formula(xformla)
    },
    error = function(e) {
      # if an error is thrown, stop execution and throw an error
      stop("xformla is not a valid formula. Please check xformla", e$message)
    }
  )

  # Creating a post dummy variable based on tlist[2] (second period = post treatment)
  tlist <- unique(dta[[args$tname]])[base::order(unique(dta[[args$tname]]))]
  dta$post <- as.numeric(dta[[tname]] == tlist[2])
  # Checking if covariates are time invariant in panel data case
  # Create the model matrix
  cov_pre <- model.matrix(as.formula(xformla), data = dta[dta[["post"]] == 0])
  cov_post <- model.matrix(as.formula(xformla), data = dta[dta[["post"]] == 1])
  # Check if covariates are time invariant
  if (!all(cov_pre == cov_post)) {
    stop("Covariates are not time invariant. Please check xformla")
  }

  # sort data based on idnam and tname and make it balanced
  data.table::setorderv(dta, c(idname, tname), c(1,1))
  dta <- BMisc::makeBalancedPanel(dta, idname, tname)


  #-------------------------------------
  # Generate the output for estimation
  # -------------------------------------

  cleaned_data <- data.table::data.table(id = dta[[idname]],
                                         y = dta[[yname]],
                                         post = dta$post,
                                         treat = dta[[dname]],
                                         partition = dta[[partition_name]],
                                         weights = dta$weights)
  # creating subgroup variable
  # 4 if (partition ==1 & treat == 1); 3 if (partition ==0 & treat == 1); 2 if (partition ==1 & treat == 0); 1 if (partition ==0 & treat == 0)
  cleaned_data$subgroup <- ifelse((cleaned_data$partition == 1) & (cleaned_data$treat == 1), 4,
                          ifelse((cleaned_data$partition == 0) & (cleaned_data$treat == 1), 3,
                          ifelse((cleaned_data$partition == 1) & (cleaned_data$treat == 0), 2, 1)))

  # Flag for not enough observations for each subgroup
  # Calculate the size of each subgroup in the 'subgroup' column
  subgroup_counts <- cleaned_data[, .N/2, by = subgroup]
  # Check if each subgroup has at least 5 observations. Check this threshold if needed.
  sufficient_obs <- all(subgroup_counts$N >= 5)
  # Stop the code if not all subgroups have at least 5 observations
  if (!sufficient_obs) {
    stop("Not enough observations in each subgroup. Please check the data.")
  }

  # Flag for small groups for inference
  # Calculate the size of each group in the 'treat' column
  # Subset 'gsize' to get the small groups based on "required size"
  gsize <- cleaned_data[, if(.N < length(BMisc::rhs.vars(xformla)) + 5) .N, by = "treat"]
  # If there are any small groups, stop execution and print an error message
  if (nrow(gsize) > 0) {
    stop("Either treatment or the comparison group in your dataset is very small. Inference is not feasible.")
  }

  # adding covariates
  if (!is.null(xformla)) {
    cleaned_data <- cbind(cleaned_data, stats::model.matrix(xformla,
                                                          data = dta,
                                                          na.action = na.pass))
  }

  # Check for missing values in several columns
  #missing_flags <- lapply(cleaned_data[, .(treat, post, y, weights)], anyNA)
  missing_flags <- lapply(cleaned_data[, .SD, .SDcols = c("treat", "post", "y", "weights")], anyNA)
  missing_X_flag <- base::anyNA(cleaned_data[,-c(1:7)]) # all except id, y, post, treat, partition, subgroup, weights

  # Print warning messages if any missing values are found
  if (any(unlist(missing_flags))) {
    # Get the names of the TRUE values
    var.with.na <- names(missing_flags)[which(unlist(missing_flags))]
    # Print warning messages
    warning(paste("Some", var.with.na, "contain missing values(NA). We are dropping those observations"))
  }
  if(missing_X_flag) warning("Some covariates are missing (NA). We are dropping those observations")

  # Remove rows with missing values
  cleaned_data <- na.omit(cleaned_data)

  # Remove collinear variables
  # Convert the remaining columns to a matrix
  cov_m <- as.matrix(cleaned_data[, -c(1:7)])
  # Use the qr() function to detect collinear columns
  qr_m <- qr(cov_m, tol = 1e-6)
  # Get the rank of the matrix
  rank_m <- qr_m$rank
  # Get the indices of the non-collinear columns
  non_collinear_indices <- qr_m$pivot[seq_len(rank_m)]
  # Drop the collinear columns from the data.table
  cleaned_data <- cleaned_data[, c(seq(1,7,1), non_collinear_indices + 7), with = FALSE]

  # drop the intercept
  cleaned_data[, 8 := NULL]

  # TODO: Add subgroup_counts to the output to present information about the subgroups

  out <- list(preprocessed_data = cleaned_data,
              xformula = xformla,
              est_method = est_method,
              learners = learners,
              n_folds = n_folds,
              boot = boot,
              boot_type = boot_type,
              nboot = nboot,
              inffunc = inffunc,
              subgroup_counts = subgroup_counts)

  return(out)
}

# Preprocess function for multiple periods case.
# TODO: Implement the function
# run_preprocess_multPeriods <-function(yname,
#                                       tname,
#                                       idname,
#                                       dname = NULL,
#                                       gname,
#                                       partition_name,
#                                       xformla = ~1,
#                                       data,
#                                       control_group,
#                                       est_method = "trad",
#                                       learners = NULL,
#                                       weightsname = NULL,
#                                       boot = FALSE,
#                                       boot_type = "multiplier",
#                                       nboot = NULL,
#                                       inffunc = FALSE) {
#   # add code here
#   out <- list(preprocessed_data = cleaned_data,
#               xformula = xformla,
#               est_method = est_method,
#               learners = learners,
#               weightsname = weightsname,
#               boot = boot,
#               boot_type = boot_type,
#               nboot = nboot,
#               inffunc = inffunc)
#
#   return(out)
# }


