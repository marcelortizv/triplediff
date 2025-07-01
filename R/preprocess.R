#' @import data.table
#' @importFrom stats model.matrix
#' @import BMisc
#' @export
NULL
#--------------------------------------------------
# Function to pre-process the data to use on ddd estimator

run_nopreprocess_2periods <- function(yname,
                                      tname,
                                      idname,
                                      gname,
                                      pname,
                                      xformla = ~1,
                                      dta,
                                      control_group = NULL,
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
                                      cores = NULL,
                                      inffunc = FALSE){

  arg_names <- setdiff(names(formals()), "dta")
  args <- mget(arg_names, sys.frame(sys.nframe()))

  # Flag for alpha > 0.10
  if (alpha > 0.10) {
    warning("alpha = ", alpha, " is too high. Using alpha = 0.05 as default.")
    alpha <- 0.05
    args$alpha <- alpha
  }

  # For dml, we only allow analytical standard errors. User needs to specify the number fo folds.
  if (est_method == "dml"){
    # allowing bootstrap only for dr.
    if (boot == TRUE){
      warning("Bootstrapping is not allowed for DML. Setting boot = FALSE.")
      boot <- FALSE
      args$boot <- boot
    }
    # number of folds
    if (n_folds < 2 || is.null(n_folds)) {
      warning("For DML estimation, n_folds must be at least 2. Setting n_folds = 2.")
      n_folds <- 2
      args$n_folds <- n_folds
    }
  }

  # setting default bootstrap reps
  if (boot == TRUE){
    if (is.null(nboot)){
      warning("Number of bootstrap samples not specified. Defaulting to 999 reps.")
      nboot <- 999
      args$nboot <- nboot
    }
  }

  # Flags for cluster variable
  if (!is.null(cluster)) {
    # dropping idname from cluster
    if (idname %in% cluster) {
      cluster <- setdiff(cluster, idname)
    }

    # flag if cluster variable is in dataset
    if (!is.element(cluster, base::colnames(dta))) {
      stop("cluster = ",cluster,  " could not be found in the data provided. Please check your arguments")
    }

    # check if user is providing more than 2 cluster variables (different than idname)
    if (length(cluster) > 1) {
      stop("You can only provide 1 cluster variable additionally to the one provided in idname. Please check your arguments")
    }
  }


  # set weights
  base::ifelse(is.null(weightsname), weights <- rep(1, nrow(dta)), weights <- dta[[weightsname]])
  # Check for missing values in the weights vector
  if (anyNA(weights)) {
    stop("There are missing values in the weights column. Please check them.")
  }
  # enforcing normalization of weights
  weights <- weights/mean(weights)
  dta$weights <- weights


  # Flag for xformla
  if (is.null(xformla)) {
    xformla <- ~1
  }

  # Creating a post dummy variable based on tlist[2] (second period = post treatment)
  tlist <- dta[, sort(unique(get(tname)))]
  dta[, post := as.numeric(get(tname) == tlist[2])]

  # sort data based on idnam and tname and make it balanced
  data.table::setorderv(dta, c(idname, tname), c(1,1))
  dta <- BMisc::makeBalancedPanel(dta, idname, tname)

  #-------------------------------------
  # Generate the output for estimation
  # -------------------------------------

  cleaned_data <- data.table::data.table(id = dta[[idname]],
                                         y = dta[[yname]],
                                         post = dta$post,
                                         treat = dta[[gname]],
                                         period = dta[[tname]],
                                         partition = dta[[pname]],
                                         weights = dta$weights)
  idx_static_vars <- 8 # useful to perform the elimination of collinear variables
  # Add cluster column if cluster argument is provided
  if (!is.null(cluster)) {
    cleaned_data[, cluster := dta[[cluster]]]
    idx_static_vars <- 9
  }

  # creating subgroup variable
  # 4 if (partition ==1 & treat == 1); 3 if (partition ==0 & treat == 1); 2 if (partition ==1 & treat == 0); 1 if (partition ==0 & treat == 0)
  # cleaned_data$subgroup <- ifelse((cleaned_data$partition == 1) & (cleaned_data$treat == 1), 4,
  #                                 ifelse((cleaned_data$partition == 0) & (cleaned_data$treat == 1), 3,
  #                                        ifelse((cleaned_data$partition == 1) & (cleaned_data$treat == 0), 2, 1)))
  cleaned_data[, subgroup := fifelse(partition == 1 & treat == 1, 4,
                                     fifelse(partition == 0 & treat == 1, 3,
                                             fifelse(partition == 1 & treat == 0, 2, 1)))]

  # Flag for not enough observations for each subgroup
  # Calculate the size of each subgroup in the 'subgroup' column
  subgroup_counts <- cleaned_data[, .N/2, by = subgroup][order(-subgroup)]

  # adding covariates
  if (!is.null(xformla)) {
    cleaned_data <- cbind(cleaned_data, stats::model.matrix(xformla,
                                                            data = dta,
                                                            na.action = na.pass))
  }

  # Remove rows with missing values
  cleaned_data <- na.omit(cleaned_data)

  # Remove collinear variables
  # Convert the remaining columns to a matrix
  cov_m <- as.matrix(cleaned_data[, -c(1:idx_static_vars), with = FALSE])
  # Use the qr() function to detect collinear columns
  qr_m <- qr(cov_m, tol = 1e-6)
  # Get the rank of the matrix
  rank_m <- qr_m$rank
  # Get the indices of the non-collinear columns
  non_collinear_indices <- qr_m$pivot[seq_len(rank_m)]
  # Drop the collinear columns from the data.table
  cleaned_data <- cleaned_data[, c(seq(1,idx_static_vars,1), non_collinear_indices + idx_static_vars), with = FALSE]

  # drop the intercept
  #cleaned_data[, (idx_static_vars+1) := NULL]
  cleaned_data[, "(Intercept)" := NULL]

  out <- list(preprocessed_data = cleaned_data,
              xformula = xformla,
              est_method = est_method,
              learners = learners,
              n_folds = n_folds,
              boot = boot,
              nboot = nboot,
              cluster = cluster,
              cband = cband,
              alpha = alpha,
              use_parallel = use_parallel,
              cores = cores,
              inffunc = inffunc,
              subgroup_counts = subgroup_counts)

  return(out)

}


# Preprocess function for 2 periods case
run_preprocess_2Periods <- function(yname,
                                   tname,
                                   idname,
                                   gname,
                                   pname,
                                   xformla = ~1,
                                   dta,
                                   control_group = NULL,
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
                                   cores = NULL,
                                   inffunc = FALSE){

  # Capture all arguments except 'data'
  arg_names <- setdiff(names(formals()), "dta")
  args <- mget(arg_names, sys.frame(sys.nframe()))

  #-------------------------------------
  # Error checking
  #-------------------------------------

  # Flag for parallel and cores
  if (boot){
    if ((use_parallel) && (is.null(cores))) {
      warning("Parallel processing is enabled but the number of cores is not specified. Using 1 core as default.")
      cores <- 1
      args$cores <- cores
    }
  }

  # Flag for alpha > 0.10
  if (alpha > 0.10) {
    warning("alpha = ", alpha, " is too high. Using alpha = 0.05 as default.")
    alpha <- 0.05
    args$alpha <- alpha
  }

  # For dml, we only allow analytical standard errors.
  if (est_method == "dml" & boot == TRUE){
    warning("Bootstrapping is not allowed for DML. Setting boot = FALSE.")
    boot <- FALSE
    args$boot <- boot
  }

  # setting default bootstrap reps
  if (boot == TRUE){
    if (is.null(nboot)){
      warning("Number of bootstrap samples not specified. Defaulting to 999 reps.")
      nboot <- 999
      args$nboot <- nboot
    }
  }

  # Run argument checks
  validate_args_2Periods(args, dta)

  # Flags for cluster variable
  if (!is.null(cluster)) {
    # dropping idname from cluster
    if (idname %in% cluster) {
      cluster <- setdiff(cluster, idname)
    }

    # flag if cluster variable is in dataset
    if (!is.element(cluster, base::colnames(dta))) {
      stop("cluster = ",cluster,  " could not be found in the data provided. Please check your arguments.")
    }

    # check if user is providing more than 2 cluster variables (different than idname)
    if (length(cluster) > 1) {
      stop("You can only provide 1 cluster variable additionally to the one provided in idname. Please check your arguments.")
    }

    # Check that cluster variables do not vary over time within each unit
    if (length(cluster) > 0) {
      # Efficiently check for time-varying cluster variables
      clust_tv <- dta[, lapply(.SD, function(col) length(unique(col)) == 1), by = id, .SDcols = cluster]
      # If any cluster variable varies over time within any unit, stop execution
      if (!all(unlist(clust_tv[, -1, with = FALSE]))) {
        stop("triplediff cannot handle time-varying cluster variables at the moment. Please check your cluster variable.")
      }
    }
  }

  # set weights
  base::ifelse(is.null(weightsname), weights <- rep(1, nrow(dta)), weights <- dta[[weightsname]])
  # Check for missing values in the weights vector
  if (anyNA(weights)) {
    stop("There are missing values in the weights column. Please check them.")
  }
  # enforcing normalization of weights
  weights <- weights/mean(weights)
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
  tlist <- dta[, sort(unique(get(tname)))]
  dta[, post := as.numeric(get(tname) == tlist[2])]

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
                                         treat = dta[[gname]],
                                         period = dta[[tname]],
                                         partition = dta[[pname]],
                                         weights = dta$weights)

  idx_static_vars <- 8 # useful to perform the elimination of collinear variables
  # Add cluster column if cluster argument is provided
  if (!is.null(cluster)) {
    cleaned_data[, cluster := dta[[cluster]]]
    idx_static_vars <- 9
  }

  # creating subgroup variable
  # 4 if (partition ==1 & treat == 1); 3 if (partition ==0 & treat == 1); 2 if (partition ==1 & treat == 0); 1 if (partition ==0 & treat == 0)
  cleaned_data[, subgroup := fifelse(partition == 1 & treat == 1, 4,
                                     fifelse(partition == 0 & treat == 1, 3,
                                             fifelse(partition == 1 & treat == 0, 2, 1)))]

  # Flag for not enough observations for each subgroup
  # Calculate the size of each subgroup in the 'subgroup' column
  subgroup_counts <- cleaned_data[, .N/2, by = subgroup][order(-subgroup)]
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
  # all except id, y, post, treat, partition, period, subgroup, weights, cluster (if any)
  missing_X_flag <- base::anyNA(cleaned_data[,-c(1:idx_static_vars), with = FALSE])
  # Print warning messages if any missing values are found
  if (any(unlist(missing_flags))) {
    # Get the names of the TRUE values
    var.with.na <- names(missing_flags)[which(unlist(missing_flags))]
    # Print warning messages
    warning(paste("Some", var.with.na, "contain missing values (NA). We are dropping those observations"))
  }
  if(missing_X_flag) warning("Some covariates are missing (NA). We are dropping those observations")

  # Remove rows with missing values
  cleaned_data <- na.omit(cleaned_data)

  # Remove collinear variables
  # Convert the remaining columns to a matrix
  cov_m <- as.matrix(cleaned_data[, -c(1:idx_static_vars), with = FALSE])
  # Use the qr() function to detect collinear columns
  qr_m <- qr(cov_m, tol = 1e-6)
  # Get the rank of the matrix
  rank_m <- qr_m$rank
  # Get the indices of the non-collinear columns
  non_collinear_indices <- qr_m$pivot[seq_len(rank_m)]
  # Drop the collinear columns from the data.table
  cleaned_data <- cleaned_data[, c(seq(1,idx_static_vars,1), non_collinear_indices + idx_static_vars), with = FALSE]

  # drop the intercept
  cleaned_data[, "(Intercept)" := NULL]

  out <- list(preprocessed_data = cleaned_data,
              xformula = xformla,
              tname = tname,
              est_method = est_method,
              learners = learners,
              n_folds = n_folds,
              boot = boot,
              nboot = nboot,
              cluster = cluster,
              cband = cband,
              alpha = alpha,
              use_parallel = use_parallel,
              cores = cores,
              inffunc = inffunc,
              subgroup_counts = subgroup_counts)

  return(out)
}

# Preprocess function for multiple periods case.
# @keywords internal
# @export
run_preprocess_multPeriods <- function(yname,
                                       tname,
                                       idname,
                                       gname,
                                       pname,
                                       xformla = ~1,
                                       dta,
                                       control_group,
                                       base_period,
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
                                       cores = NULL,
                                       inffunc = FALSE){

  # Capture all arguments except 'data'
  arg_names <- setdiff(names(formals()), "dta")
  args <- mget(arg_names, sys.frame(sys.nframe()))

  #-------------------------------------
  # Error checking
  #-------------------------------------

  # Flag for parallel and cores
  if (boot){
    if ((use_parallel) && (is.null(cores))) {
      warning("Parallel processing is enabled but the number of cores is not specified. Using 1 core.")
      cores <- 1
      args$cores <- cores
    }
  }

  # Flag for alpha > 0.10
  if (alpha > 0.10) {
    warning("alpha = ", alpha, " is too high. Using alpha = 0.05 as default.")
    alpha <- 0.05
    args$alpha <- alpha
  }

  # For dml, we only allow analytical standard errors.
  if (est_method == "dml" & boot == TRUE){
    warning("Bootstrapping is not currently allowed for DML. Setting boot = FALSE.")
    boot <- FALSE
    args$boot <- boot
  }

  # setting default bootstrap reps
  if (boot == TRUE){
    if (is.null(nboot)){
      warning("Number of bootstrap samples not specified. Defaulting to 999 reps.")
      nboot <- 999
      args$nboot <- nboot
    }
  }

  # Run argument checks
  validate_args_multPeriods(args, dta)

  # Flags for cluster variable
  if (!is.null(cluster)) {
    # dropping idname from cluster
    if (idname %in% cluster) {
      cluster <- setdiff(cluster, idname)
    }

    # flag if cluster variable is in dataset
    if (!is.element(cluster, base::colnames(dta))) {
      stop("cluster = ",cluster,  " could not be found in the data provided. Please check your arguments.")
    }

    # check if user is providing more than 2 cluster variables (different than idname)
    if (length(cluster) > 1) {
      stop("You can only provide 1 cluster variable additionally to the one provided in idname. Please check your arguments.")
    }

    # Check that cluster variables do not vary over time within each unit
    if (length(cluster) > 0) {
      # Efficiently check for time-varying cluster variables
      clust_tv <- dta[, lapply(.SD, function(col) length(unique(col)) == 1), by = id, .SDcols = cluster]
      # If any cluster variable varies over time within any unit, stop execution
      if (!all(unlist(clust_tv[, -1, with = FALSE]))) {
        stop("triplediff cannot handle time-varying cluster variables at the moment. Please check your cluster variable.")
      }
    }
  }

  # set in-blank xformla if no covariates are provided
  if (is.null(xformla)) {
    xformla <- ~1
  }

  # keep relevant columns in data
  cols_to_keep <- c(idname, tname, yname, gname, pname, weightsname, cluster)

  model_frame <- model.frame(xformla, data = dta, na.action = na.pass)
  # Subset the data.table to keep only relevant columns
  dta <- dta[, ..cols_to_keep]

  # Bind the model frame columns
  dta <- cbind(dta, as.data.table(model_frame))

  # Check if any covariates were missing
  n_orig <- dta[, .N]
  dta <- dta[complete.cases(dta)]
  n_new <- dta[, .N]
  n_diff <- n_orig - n_new
  if (n_diff != 0) {
    warning(paste0("Dropped ", n_diff, " rows from original data due to missing data."))
  }

  # set weights
  base::ifelse(is.null(weightsname), weights <- rep(1, n_new), weights <- dta[[weightsname]])
  # enforcing normalization of weights. At this point we already drop any missing values in weights.
  weights <- weights/mean(weights)
  dta$weights <- weights

  # get a list of dates from min to max
  tlist <- dta[, sort(unique(get(tname)))]

  # Coerce control group identified with Inf as zero
  dta[get(gname) == Inf, (gname) := 0]

  # Identify groups with treatment time bigger than the maximum treatment time
  # calculate the maximum treatment time
  max_treatment_time <- max(tlist, na.rm = TRUE)
  dta[, asif_never_treated := (get(gname) > max_treatment_time)]
  # replace NA values with FALSE in the logical vector
  dta[is.na(asif_never_treated), asif_never_treated := FALSE]
  # set gname to 0 for those groups considered as never treated
  dta[asif_never_treated == TRUE, (gname) := 0]
  # remove the temporary column
  dta[, asif_never_treated := NULL]

  # get list of treated groups (by time) from min to max
  glist <- sort(unique(dta[[gname]]))

  # Check if there is a never treated group
  if (0 %in% glist == FALSE) {
    if (control_group == "nevertreated") {
      stop("There is no available never-treated group")
    } else {
      # Drop all time periods with time periods >= latest treated
      latest_treated_time <- max(glist)
      # with anticipation?
      #latest_treated_time <- max(glist) - anticipation
      dta <- dta[get(tname) < latest_treated_time]

      tlist <- sort(unique(dta[[tname]]))
      glist <- sort(unique(dta[[gname]]))

      # don't compute ATT(g,t) for groups that are only treated at end
      # and only play a role as a comparison group
      glist <- glist[glist < max(glist)]
    }
  }

  # Get the first period
  first_period <- tlist[1]

  # Filter the treated groups
  # glist <- glist[glist > 0 & glist > (first_period + anticipation)]
  glist <- glist[glist > 0 & glist > first_period]

  # Check for groups treated in the first periods and drop them
  # identify groups treated in the first period
  dta[, treated_first_period := (get(gname) <= first_period) & (get(gname) != 0)]
  dta[is.na(treated_first_period), treated_first_period := FALSE]

  # count the number of units treated in the first period
  nfirstperiod <- uniqueN(dta[treated_first_period == TRUE, get(idname)])

  # handle units treated in the first period
  if (nfirstperiod > 0) {
    warning(paste0("Dropped ", nfirstperiod, " units that were already treated in the first period."))
    dta <- dta[get(gname) %in% c(0, glist)]

    # update tlist and glist
    tlist <- dta[, sort(unique(get(tname)))]
    glist <- dta[, sort(unique(get(gname)))]
    glist <- glist[glist > 0]

    # Drop groups treated in the first period or before
    first_period <- tlist[1]
    #glist <- glist[glist > first_period + anticipation]
    glist <- glist[glist > first_period]
  }
  # remove the temporary column
  dta[, treated_first_period := NULL]

  # ----------------------------------------------------------------
  # panel data only
  # ----------------------------------------------------------------

  # Check for complete cases
  keepers <- complete.cases(dta)
  n <- uniqueN(dta[[idname]])
  n_keep <- uniqueN(dta[keepers, get(idname)])

  if (nrow(dta[keepers, ]) < nrow(dta)) {
    warning(paste0("Dropped ", (n - n_keep), " observations that had missing data."))
    dta <- dta[keepers, ]
  }

  # Make it a balanced data set
  n_old <- uniqueN(dta[[idname]])
  row_orig <- dta[, .N]
  dta <- BMisc::makeBalancedPanel(dta, idname, tname)
  n <- uniqueN(dta[[idname]])
  row_new <- dta[, .N]
  row_diff <- row_orig - row_new


  if (n < n_old) {
    warning(paste0("Dropped observations from ", n_old - n, " units ", "(", row_diff ," rows)" ," while converting to balanced panel."))
  }

  #-------------------------------------
  # Generate the output for estimation
  # ------------------------------------
  cleaned_data <- data.table::data.table(id = dta[[idname]],
                                         y = dta[[yname]],
                                         first_treat = dta[[gname]],
                                         period = dta[[tname]],
                                         partition = dta[[pname]],
                                         weights = dta$weights)

  # Add cluster column if cluster argument is provided
  if (!is.null(cluster)) {
    cleaned_data[, cluster := dta[[cluster]]]
  }

  # If drop all data, you do not have a panel.
  if (nrow(cleaned_data) == 0) {
    stop("All observations dropped to convert data to balanced panel. Please check your argument in 'idname'.")
  }

  n <- nrow(cleaned_data[period == tlist[1], ])

  # Check if groups is empty
  if(length(glist)==0){
    stop("No valid groups. The variable in 'gname' should be expressed as the time a unit is first treated (0 if never-treated).")
  }

  # Check for small comparison groups
  # Calculate the size of each group in the 'treat' column
  gsize <- cleaned_data[, .N / length(tlist), by = first_treat][order(-first_treat)]
  # Calculate the required size
  reqsize <- length(BMisc::rhs.vars(xformla)) + 5
  # Identify groups to warn about
  small_groups <- gsize[V1 < reqsize]
  # Warn if some groups are small
  if (nrow(small_groups) > 0) {
    gpaste <- paste(small_groups[, treat], collapse = ",")
    warning(paste0("Be aware that there are some small groups in your dataset.\n  Check groups: ", gpaste, "."))

    if (0 %in% small_groups[, treat] & control_group == "nevertreated") {
      stop("Never treated group is too small, try setting control_group = \"notyettreated\"")
    }
  }

  # How many time periods
  nT <- length(tlist)
  # How many treated groups
  nG <- length(glist)

  # adding covariates
  if (!is.null(xformla)) {
    cleaned_data <- cbind(cleaned_data, stats::model.matrix(xformla,
                                                            data = dta,
                                                            na.action = na.pass))
  }
  # drop the intercept
  cleaned_data[, "(Intercept)" := NULL]

  # order dataset wrt idname and tname
  setorder(cleaned_data, "id", "period")

  out <- list(preprocessed_data = cleaned_data,
              xformula = xformla,
              tname = tname,
              est_method = est_method,
              boot = boot,
              nboot = nboot,
              cluster = cluster,
              cband = cband,
              alpha = alpha,
              use_parallel = use_parallel,
              cores = cores,
              control_group = control_group,
              base_period = base_period,
              n = n,
              cohorts = nG,
              time_periods = nT,
              cohort_size = gsize,
              tlist = tlist,
              glist = glist)

  return(out)
}

# @keywords internal
# @export
# Process results inside att_gt_dr function
process_attgt <- function(attgt_list){
  groups <- length(unique(unlist(BMisc::getListElement(attgt_list, "group"))))
  time_periods <- length(unique(unlist(BMisc::getListElement(attgt_list, "year"))))

  # empty vectors to hold results
  group <- c()
  att <- c()
  periods <- c()

  i <- 1
  # extract results and populate vectors
  for (g in 1:groups) {
    for (t in 1:time_periods) {
      group[i] <- attgt_list[[i]]$group
      periods[i] <- attgt_list[[i]]$year
      att[i] <- attgt_list[[i]]$att
      i <- i + 1
    }
  }

  return(list(group = group, att = att, periods = periods))
}


