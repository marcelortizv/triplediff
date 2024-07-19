#' Function to compute nuisances parameters and did
#' @importFrom stats update
#' @import parglm
#' @import speedglm
#' @import data.table
#' @import mlr3
#' @import mlr3learners
#' @import mlr3tuning
#' @noRd
#--------------------------------------------------

# utility function to generate equation to estimate pscore
get_formula_pscore <- function(xformula, weights = TRUE){
  # Get a formula object for pscore estimation
  if (weights) {
    formula_obj <- stats::update(xformula, PA4 ~ . + weights)
  } else {
    formula_obj <- stats::update(xformula, PA4 ~ .)
  }
  return(formula_obj)
}

# utility function to generate equation to estimate outcome regression
get_formula_reg <- function(xformula, weights = TRUE){
  # Get a formula object for outcome regression
  if (weights) {
    formula_obj <- stats::update(xformula, deltaY ~ . + weights)
  } else {
    formula_obj <- stats::update(xformula, deltaY ~ .)
  }
  return(formula_obj)
}

# Utility function to reshape the dataset from long to wide format
get_wide_data <- function(dt) {
  # Ensure the input is a data.table
  if (!inherits(dt, "data.table")) {
    dt <- as.data.table(dt)
  }

  # Extract the names of the invariant columns
  invariant_cols <- setdiff(names(dt), c("id", "y", "post"))

  # Create a formula for dcast that includes all invariant columns
  formula <- paste("id +", paste(invariant_cols, collapse = " +"), "~ post")

  # Reshape from long to wide format using dcast
  reshaped_dt <- dcast(dt, formula, value.var = "y")

  # rename the columns for clarity
  setnames(reshaped_dt, old = c("0", "1"), new = c("y0", "y1"))

  return(reshaped_dt)
}

# Utility function to get long vector of scores Psi hat
get_long_scores <- function(dmlddd_scores_hat_k){
  scores_list <- lapply(dmlddd_scores_hat_k, function(x) x$score_hat)
  return(unlist(scores_list))
}

# Define a function to extract the 'ddd_k' value and compute its mean
mean_ddd_k <- function(dmlddd_scores_hat_k) {
  # Use lapply to extract the 'ddd_k' values
  ddd_k_values <- lapply(dmlddd_scores_hat_k, function(x) x$ddd_k)
  return(mean(unlist(ddd_k_values)))
}

# Function to compute propensity scores using parglm for multiple subgroups
compute_pscore <- function(data, condition_subgroup, xformula) {
  # get formula for pscore estimation using covariates
  formula_pscore <- get_formula_pscore(xformula, weights = FALSE)
  # Subset data for condition_subgroup and subgroup == 4 or the given condition_subgroup
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]
  # Adding treatment variable P(1{PA4 = 4}|X)
  condition_data[, "PA4" := ifelse(condition_data$subgroup == 4, 1, 0)]
  uid_condition_data <- unique(condition_data, by = "id")

  # Fit logistic regression model using parglm
  model <- parglm::parglm(formula_pscore, data = uid_condition_data,
                          family = stats::binomial(),
                          weights = weights,
                          control = parglm.control(nthreads = data.table::getDTthreads()),
                          intercept = FALSE)

  # Flag for convergence of glm
  if (!model$converged) {
    warning(paste("Logistic regression model for subgroup", condition_subgroup,
                  "did not converge."))
  }

  # Flag for weird coefficients
  if(anyNA(model$coefficients)) {
    stop(paste("Pscore model coefficients for subgroup", condition_subgroup,
               "have NA components. Multicollinearity of covariates is a likely reason"))
  }

  # Compute propensity scores
  propensity_scores <- predict(model, newdata = uid_condition_data, type = "response")

  # Warning for overlap condition
  if (any(propensity_scores < 0.005)) {
    warning(paste("Propensity scores for comparison subgroup", condition_subgroup,
                  "have poor overlap."))
  }

  # Avoid divide by zero
  propensity_scores <- pmin(propensity_scores, 1 - 1e-16)

  # Save Hessian matrix for influence function
  hessian_matrix <- stats::vcov(model) * nrow(uid_condition_data)

  return(list(propensity_scores = propensity_scores, hessian_matrix = hessian_matrix))
}

# Function to compute outcome regression for multiple subgroups
compute_outcome_regression <- function(data, condition_subgroup, xformula){
  # Subset data for condition_subgroup and subgroup == 4 or the given condition_subgroup
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]
  # Subset data for the control group (subgroup != 4)
  control_data <- condition_data[condition_data$subgroup == condition_subgroup]
  y1_control = control_data[control_data$post == 1, y]
  y0_control = control_data[control_data$post == 0, y]
  # generate deltaY for control group
  deltaY_control = y1_control - y0_control
  # get covariates including the intercept (conditioning in pre-treatment period since covariates are invariant)
  cov_control <- stats::model.matrix(as.formula(xformula), data = control_data[control_data$post == 0])
  # get weights (conditioning in pre-treatment period since weights are invariant)
  i_weights = control_data[control_data$post == 0, weights]

  # get coefficients of outcome regression model for the control group (subgroup != 4)
  # Attempt to compute reg.coeff using speedglm
  try_speedglm <- tryCatch({
    reg.coeff <- stats::coef(speedglm::speedglm.wfit(y = deltaY_control,
                                                     X = as.matrix(cov_control),
                                                     intercept = FALSE,
                                                     weights = i_weights))
  }, error = function(e) {
    # If error, attempt to compute reg.coeff using lm.wfit
    reg.coeff <- tryCatch({
      stats::coef(stats::lm.wfit(x = as.matrix(cov_control),
                                 y = deltaY_control,
                                 w = i_weights))
    }, error = function(e2) {
      # If error with lm.wfit, stop the program and send an error message
      stop("Error in computing regression coefficients: Subgroup ", condition_subgroup, " may have insufficient data.")
    })
  })

  # Flag for NA coefficients
  if(anyNA(reg.coeff)) {
    stop(paste("Outcome regression model coefficients for subgroup", condition_subgroup,
               "have NA components. Multicollinearity of covariates is a likely reason"))
  }

  # compute regression adjustment
  y1 = condition_data[condition_data$post == 1, y]
  y0 = condition_data[condition_data$post == 0, y]
  deltaY = y1 - y0
  # get covariates including the intercept (conditioning in pre-treatment period since covariates are invariant)
  covX <- stats::model.matrix(xformula, data = condition_data[condition_data$post == 0])
  # compute OR delta
  or_delta = as.vector(tcrossprod(reg.coeff, as.matrix(covX)))
  # compute regression adjustment
  # reg_adjust = deltaY - or_delta
  return(list(deltaY = deltaY, or_delta = or_delta))
}

# Function to compute the average treatment effect for multiple subgroups
compute_did <- function(data, condition_subgroup, pscores, reg_adjustment, xformula){

  data <- unique(data, by = "id")
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]
  PA4 = ifelse(condition_data$subgroup == 4, 1, 0)
  PAa = ifelse(condition_data$subgroup == condition_subgroup, 1, 0)

  # Compute propensity scores
  if (condition_subgroup == 3) {
    pscore <- pscores[[1]]$propensity_scores
    hessian <- pscores[[1]]$hessian_matrix
    deltaY <- reg_adjustment[[1]]$deltaY
    or_delta <- reg_adjustment[[1]]$or_delta

  } else if (condition_subgroup == 2) {
    pscore <- pscores[[2]]$propensity_scores
    hessian <- pscores[[2]]$hessian_matrix
    deltaY <- reg_adjustment[[2]]$deltaY
    or_delta <- reg_adjustment[[2]]$or_delta

  } else if (condition_subgroup == 1) {
    pscore <- pscores[[3]]$propensity_scores
    hessian <- pscores[[3]]$hessian_matrix
    deltaY <- reg_adjustment[[3]]$deltaY
    or_delta <- reg_adjustment[[3]]$or_delta

  } else {
    stop("Invalid condition_subgroup")
  }

  # get weights
  i_weights = condition_data[condition_data$post == 0, weights]

  ################################
  # Get doubly-robust estimation #
  ################################

  w_treat = i_weights * PA4
  w_control = (i_weights * pscore * PAa) / (1 - pscore)
  riesz_treat = w_treat * (deltaY - or_delta)
  riesz_control = w_control * (deltaY - or_delta)
  att_treat = mean(riesz_treat, na.rm = TRUE) / mean(w_treat, na.rm = TRUE)
  att_control = mean(riesz_control, na.rm = TRUE) / mean(w_control, na.rm = TRUE)

  dr_att = att_treat - att_control

  ##########################
  # Get influence function #
  ##########################

  # Influence function related to the estimation of pscores
  covX = stats::model.matrix(xformula, data = condition_data)
  M2 <- base::colMeans(w_control * (deltaY - or_delta - att_control) * covX, na.rm = TRUE) # reg_adjust = deltaY - m_delta(x)

  score_ps <- i_weights * (PA4 - pscore) * covX
  # asymptotic linear representation of logit's beta
  score_ps_no_na <- na.omit(score_ps) # Exclude rows with NA values
  asy_lin_rep_ps <- score_ps_no_na %*% hessian
  inf_control_pscore <- asy_lin_rep_ps %*% as.matrix(M2)

  # Influence function related to the estimation of pscores

  M1 <- base::colMeans(w_treat * covX, na.rm = TRUE)
  M3 <- base::colMeans(w_control * covX, na.rm = TRUE)

  # Influence function related to the estimation of regression model
  or_x <- i_weights * PAa * covX
  or_ex <- i_weights * PAa * (deltaY - or_delta) * covX
  XpX <- crossprod(or_x, covX)/nrow(condition_data)

  #asymptotic linear representation of the beta
  asy_linear_or <- t(solve(XpX, t(or_ex)))

  #or for treat
  inf_treat_or <- -asy_linear_or %*% M1
  #or for control
  inf_cont_or <- -asy_linear_or %*% M3

  # Influence function from did
  inf_control_did <- riesz_control - w_control*att_control
  inf_treat_did <- riesz_treat - w_treat*att_treat

  # Influence function for the ATT
  inf_control <- (inf_control_did + inf_control_pscore + inf_cont_or)/mean(w_control)
  inf_treat <- (inf_treat_did + inf_treat_or)/mean(w_treat)
  # putting all together
  inf_func <- inf_treat - inf_control

  # fill zeros in influence function for observation outside the subgroup analyzed
  data[, "inf_func_result" := numeric(.N)]
  data[data$subgroup %in% c(condition_subgroup, 4), "inf_func_result" := inf_func]

  return(list(dr_att = dr_att, inf_func = data[["inf_func_result"]]))
}

# ---------------------------- #
# FUNCTIONS FOR DML
# ---------------------------- #

# Function to get components of a linear score
get_score_elements <- function(y, d, p_hat, m_hat){
  #' Compute the score elements for the ATT
  #' @param y Outcome variable
  #' @param d Treatment variable
  #' @param p_hat Propensity score
  #' @param m_hat Outcome regression

  # check if treatment variable is a factor
  if (is.factor(d)) {
    d <- as.numeric(as.character(d))
  }

  reg_adjustment <- y - m_hat
  weight_psi_a <- d / mean(d)
  weight_pscore <- (1 - d) * (p_hat / (1 - p_hat))
  weight_residual <- (d / mean(d)) - (weight_pscore / mean(weight_pscore))

  # score elements
  psi_a <- - weight_psi_a
  psi_b <- weight_residual * reg_adjustment

  return(list(psi_a = psi_a, psi_b = psi_b))
}

# Function to get the predicted score for each k fold
compute_scores <- function(ml_pa, ml_md, condition_subgroup){
  #' Function to compute scores based on first step nuisances predictions
  #' @param ml_pa MLR3 object for propensity score estimation
  #' @param ml_md MLR3 object for regression adjustment estimation
  #' @param condition_subgroup Subgroup to analyze
  #' @return: A list of scores and id for each k fold

  # Warning for overlap condition
  if (any(ml_pa$prediction()$prob[,2] < 0.01) | any(ml_pa$prediction()$prob[,2] > 0.99)) {
    stop(paste("Propensity scores for subgroup", condition_subgroup,
                  "have poor overlap."))
  }

  # Extract the predictions
  pscores_pred <- ml_pa$predictions()
  reg_pred <- ml_md$predictions()
  # Check that the list lengths are equal
  if (length(pscores_pred) != length(reg_pred)) {
    stop("Lists of predictions for pscore and regression must be of equal length.")
  }

  dml_scores <- vector("list", length(pscores_pred))

  for (k in seq_along(pscores_pred)) {
    # Get pscores
    pscores <- pscores_pred[[k]]$prob[,2]

    # Avoid divide by zero
    pscores <- pmin(pscores, 1 - 1e-16)
    pscores <- pmax(pscores, 1e-16)

    # Get the predictions for the k-th fold
    pscores_k_dt <- as.data.table(list(row_ids = pscores_pred[[k]]$row_ids,
                                       d = pscores_pred[[k]]$truth,
                                       p_hat = pscores))
    reg_k_dt <- as.data.table(list(row_ids = reg_pred[[k]]$row_ids,
                                   y = reg_pred[[k]]$truth,
                                   m_hat = reg_pred[[k]]$response))

    # Merge the predictions for each k-th fold
    predk_merged <- merge(pscores_k_dt, reg_k_dt, by = "row_ids", all = TRUE)

    # Compute the score elements
    scores_k <- get_score_elements(predk_merged$y, predk_merged$d, predk_merged$p_hat, predk_merged$m_hat)

    # Store the scores
    dml_scores[[k]] <- c(idk= list(predk_merged$row_ids), scores_k)
  }

  return(dml_scores)
}

# Function to compute dmlddd + Psi hat
compute_dml = function(scores) {
  psi_a <- scores$psi_a
  psi_b <- scores$psi_b
  N = length(psi_a)
  dml_ddd = -sum(psi_b) / sum(psi_a)
  psi = dml_ddd * psi_a + psi_b
  Psi = - psi / mean(psi_a)
  return(list(idk = scores$idk, ddd_k = dml_ddd, score_hat = Psi))
}


# Function to compute the DML ATT + score function
compute_dml_nuisances <- function(data, condition_subgroup, xformula, ml_pa, ml_md, n_folds) {
  #' Compute the propensity score for the given condition subgroup
  #' @param data A data.table containing the data processed
  #' @param condition_subgroup The condition subgroup for which to compute the propensity score
  #' @param xformula The formula for the propensity score model
  #' @param ml_pa The machine learning algorithm to use for the propensity score model
  #' @param ml_md The machine learning algorithm to use for the regression adjustment model
  #' @param n_folds The number of folds for cross-fitting
  #' @return A list containing ids, estimator for each k fold, and scores

  set.seed(123)

  formula_pscore <- get_formula_pscore(xformula, weights = TRUE)
  formula_reg <- get_formula_reg(xformula, weights = TRUE)

  # Subset data for subgroup == 4 or the given condition_subgroup
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]

  # get wide panel
  condition_data <- get_wide_data(condition_data)

  # Adding treatment variable P(1{PA4 = 4}|X)
  condition_data[, "PA4" := as.factor(ifelse(condition_data$subgroup == 4, 1, 0))]
  # Adding outcome variable E[deltaY | X]
  condition_data[, deltaY:= y1 - y0]

  pscore_condition_data <- as.data.table(stats::model.frame(formula_pscore, data = condition_data))
  reg_condition_data <- as.data.table(stats::model.frame(formula_reg, data = condition_data))

  # Initialize an empty vector to store propensity scores
  propensity_scores <- numeric(nrow(pscore_condition_data))
  # Initialize an empty vector to store the regression adjustment estimates
  regression_adjustment <- numeric(nrow(reg_condition_data))

  if ("weights" %in% ml_pa$properties){
    # Set weights for the classification task
    task_classif <- TaskClassif$new("pscore_task", backend = pscore_condition_data, target = "PA4")
    task_classif$set_col_roles("weights", "weight")
  } else {
    task_classif <- TaskClassif$new("pscore_task", backend = pscore_condition_data, target = "PA4")
    warning("The learner provided for propensity score estimation does not support weights.")
  }

  if ("weights" %in% ml_md$properties){
    # Set weights for the regression task
    task_regr <- TaskRegr$new("reg_task", backend = reg_condition_data, target = "deltaY")
    task_regr$set_col_roles("weights", "weight")
  } else {
    task_regr <- TaskRegr$new("reg_task", backend = reg_condition_data, target = "deltaY")
    warning("The learner provided for regression adjustment does not support weights.")
  }

  # Shared resampling method: k-fold cross-validation
  resampling <- rsmp("cv", folds = n_folds)
  resampling$instantiate(task_regr)  # Instantiate with one task, used for both

  # Perform cross-fitting for classification
  ml_classif <- resample(task_classif, ml_pa, resampling)
  # Perform cross-fitting for regression
  ml_regr <- resample(task_regr, ml_md, resampling)

  # Compute the scores for each k fold
  dml_scores <- compute_scores(ml_pa = ml_classif, ml_md = ml_regr, condition_subgroup = condition_subgroup)

  # Get a list of dmlddd estimator + scores, for each k-fold
  ddd_dml <- lapply(dml_scores, compute_dml)

  return(ddd_dml)
}

# Function to compute se for DML ATT
compute_se_dml <- function(dml_scores1, dml_scores2, dml_scores3) {
  scores <- c(dml_scores1, dml_scores2, dml_scores3)
  N = length(scores)
  sigma2 = stats::var(scores)
  se = sqrt(sigma2 / N)
  return(se)
}

# ---------------------------- #
# FUNCTIONS FOR AGGREGATION
# ---------------------------- #

#' @title Compute extra term in influence function due to estimating weights
#'
#' @description A function to compute the extra term that shows up in the
#'  influence function for aggregated treatment effect parameters
#'  due to estimating the weights
#'
#' @param keepers a vector of indices for which group-time average
#'  treatment effects are used to compute a particular aggregated parameter
#' @param pg a vector with same length as total number of group-time average
#'  treatment effects that contains the probability of being in particular group
#' @param weights additional sampling weights (nx1)
#' @param G vector containing which group a unit belongs to (nx1)
#' @param group vector of groups
#'
#' @return nxk influence function matrix
#'
#' @keywords internal
get_weight_influence <- function(keepers, pg, weights, G, group) {
  # note: weights are all of the form P(G=g|cond)/sum_cond(P(G=g|cond))
  # this is equal to P(G=g)/sum_cond(P(G=g)) which simplifies things here

  # effect of estimating weights in the numerator
  if1 <- sapply(keepers, function(k) {
    (weights * BMisc::TorF(G==group[k]) - pg[k]) /
      sum(pg[keepers])
  })
  # effect of estimating weights in the denominator
  if2 <- base::rowSums( sapply( keepers, function(k) {
    weights * BMisc::TorF(G==group[k]) - pg[k]
  })) %*%
    t(pg[keepers]/(sum(pg[keepers])^2))

  # return the influence function for the weights
  return(if1 - if2)
}

#' @title Get an influence function for particular aggregate parameters
#'
#' @title This is a generic internal function for combining influence
#'  functions across ATT(g,t)'s to return an influence function for
#'  various aggregated treatment effect parameters.
#'
#' @param att vector of group-time average treatment effects
#' @param inf_func influence function for all group-time average treatment effects
#'  (matrix)
#' @param whichones which elements of att will be used to compute the aggregated
#'  treatment effect parameter
#' @param weights_agg the weights to apply to each element of att(whichones);
#'  should have the same dimension as att(whichones)
#' @param wif extra influence function term coming from estimating the weights;
#'  should be n x k matrix where k is dimension of whichones
#'
#' @return nx1 influence function
#'
#' @keywords internal
get_agg_inf_func <- function(att, inf_func, whichones, weights_agg, wif=NULL) {
  # enforce weights are in matrix form
  weights_agg <- as.matrix(weights_agg)

  # multiplies influence function times weights and sums to get vector of weighted IF (of length n)
  agg_inf_func <- inf_func[,whichones] %*% weights_agg

  # Incorporate influence function of the weights
  if (!is.null(wif)) {
    agg_inf_func <- agg_inf_func + wif%*%as.matrix(att[whichones])
  }

  # return influence function
  return(agg_inf_func)
}

#' @title Take influence function and compute standard errors
#'
#' @description Function to take an nx1 influence function and return
#'  a standard error
#'
#' @param influence_function An influence function
#'
#' @return scalar standard error
#'
#' @keywords internal
compute_se_agg <- function(influence_function, DIDparams=NULL) {
  alpha <- .05
  boot <- FALSE
  n <- length(influence_function)

  # if (!is.null(DIDparams)) {
  #   boot <- DIDparams$boot
  #   alpha <- DIDparams$alpha
  #   cband <- DIDparams$cband
  # }

  if (boot) {
    # TODO; MAKE AVAILABLE MULTIPLIER BOOTSTRAP
    # bout <- mboot(influence_function, DIDparams)
    # return(bout$se)
    stop("Bootstrapping not yet implemented")
  } else {
    return(sqrt( mean((influence_function)^2)/n ))
  }
}

#' @title Compute Aggregated Treatment Effect Parameters
#'
#' @description Does the heavy lifting on computing aggregated group-time
#'  average treatment effects

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
#'  with fixed probability `1 - alpha`.  In order to compute uniform confidence
#'  bands, `boot` must also be set to `TRUE`.  The default is
#'  the value set in the ddd object
#'
#' @return Aggregation object (list) of class [`agg_ddd`]
#'
#' @keywords internal
#'
#' @export
compute_aggregation <- function(ddd_obj,
                                type = "simple",
                                balance_e = NULL,
                                min_e = -Inf,
                                max_e = Inf,
                                na.rm = FALSE,
                                boot = FALSE,
                                nboot = NULL,
                                cband = NULL,
                                alpha = 0.05){
  # check if the object is of class `ddd`
  if (!inherits(ddd_obj, "ddd")) {
    stop("Object must be of class `ddd`")
  }

  # check if the object is a multiple period object
  if (!(ddd_obj$argu$multiple_periods)) {
    stop("Object must be a multiple period object")
  }

  # validate types of aggregation
  if (!(type %in% c("simple", "eventstudy", "group", "calendar"))){
    stop("type must be one of 'simple', 'eventstudy', 'group' or 'calendar'")
  }

  # get parameters
  groups <- ddd_obj$groups
  periods <- ddd_obj$periods
  ATT <- ddd_obj$ATT
  inf_func_mat <- ddd_obj$inf_func_mat
  n <- ddd_obj$n
  data <- ddd_obj$data
  tlist <- ddd_obj$tlist
  glist <- ddd_obj$glist
  dta <- ddd_obj$first_period_dta
  yname <- ddd_obj$argu$yname
  partition_name <- ddd_obj$argu$partition_name
  control_group <- ddd_obj$argu$control_group

  if((na.rm == FALSE) && base::anyNA(ATT)) stop("Missing values at att_gt found. If you want to remove these, set `na.rm = TRUE'.")

  # removing NA values from ATT(g,t), groups and periods objects
  if(na.rm){
    notna <- !is.na(ATT)
    groups <- groups[notna]
    periods <- periods[notna]
    ATT <- ATT[notna]
    inf_func_mat <- inf_func_mat[, notna]
    glist <- sort(unique(groups))

    # Ensure we have non-missing post-treatment ATTs for each group if the type is "group"
    if(type == "group"){
      # Get the groups that have some non-missing ATT(g,t) in post-treatmemt periods
      gnotna <- sapply(glist, function(g) {
        # look at post-treatment periods for group g
        whichg <- which((groups == g) & (g <= periods))
        attg <- ATT[whichg]
        group_select <- !is.na(mean(attg))
        return(group_select)
      })
      gnotna <- glist[gnotna]
      # indicator for not all post-treatment ATT(g,t) missing
      not_all_na <- groups %in% gnotna
      # Re-do the na.rm thing to update the groups
      groups <- groups[not_all_na]
      periods <- periods[not_all_na]
      ATT <- ATT[not_all_na]
      inf_func_mat <- inf_func_mat[, not_all_na]
      glist <- sort(unique(groups))
    }
  }

  #-----------------------------------------------------------------------------
  # data manipulation
  #-----------------------------------------------------------------------------

  orig_periods <- periods
  orig_group <- groups
  orig_glist <- glist
  orig_tlist <- tlist
  # In case g's are not part of tlist
  orig_gtlist <- sort(unique(c(orig_tlist, orig_glist)))
  uniquet <- seq(1,length(unique(orig_gtlist)))

  # function to switch from "new" t values to  original t values
  t2orig <- function(t) {
    origt <- unique(c(orig_gtlist,0))[which(c(uniquet,0)==t)]
    return(origt)
  }

  # function to switch between "original" t values and new t values
  orig2t <- function(orig) {
    new_t <- c(uniquet,0)[which(unique(c(orig_gtlist,0))==orig)]
    out <- ifelse(length(new_t) == 0, NA, new_t)
    return(out)
  }

  t <- sapply(orig_periods, orig2t)
  group <- sapply(orig_group, orig2t)
  glist <- sapply(orig_glist, orig2t)
  tlist <- unique(t)
  maxT <- max(t)

  weights = dta[["weights"]]

  # compute probability of belonging of each group in glist
  pg <- sapply(orig_glist, function(g) mean(weights*(dta[["first_treat"]]==g)))
  pgg <- pg

  # making pg have the same length of ATT(g,t)
  pg <- pg[match(group, glist)]

  # which group time average treatment effects are post-treatment
  keepers <- which(group <= t & t <= (group + max_e))

  # n x 1 vector of group variable
  G <-  unlist(lapply(dta[["first_treat"]], orig2t))

  #-----------------------------------------------------------------------------
  # Compute the simple ATT summary
  #-----------------------------------------------------------------------------

  if (type == "simple") {

    # simple att
    # averages all post-treatment ATT(g,t) with weights
    # given by group size
    simple.att <- sum(ATT[keepers]*pg[keepers])/(sum(pg[keepers]))
    if(is.nan(simple.att)) simple.att <- NA

    # get the part of the influence function coming from estimated weights
    simple.wif <- get_weight_influence(keepers, pg, weights, G, group)

    # get the overall influence function
    simple.if <- get_agg_inf_func(att=ATT,
                                  inf_func=inf_func_mat,
                                  whichones=keepers,
                                  weights_agg=pg[keepers]/sum(pg[keepers]),
                                  wif=simple.wif)
    # Make it as vector
    simple.if <- as.numeric(simple.if)

    # get standard errors from overall influence function
    simple.se <- compute_se_agg(simple.if, dp)
    if(!is.na(simple.se)){
      if(simple.se <= sqrt(.Machine$double.eps)*10) simple.se <- NA
    }


    return(list(overall.att = simple.att,
                overall.se = simple.se,
                type = type,
                yname = yname,
                partition_name = partition_name,
                control_group = control_group,
                inf.function = list(simple.att = simple.if)))
  }

  #-----------------------------------------------------------------------------
  # Compute the group (i.e., selective) treatment timing estimators
  #-----------------------------------------------------------------------------

  if (type == "group") {

    # get group specific ATTs
    # note: there are no estimated weights here
    selective.att.g <- sapply(glist, function(g) {
      # look at post-treatment periods for group g
      whichg <- which( (group == g) & (g <= t) & (t <= (group + max_e))) # added last condition to allow for limit on longest period included in att
      attg <- ATT[whichg]
      mean(attg)
    })
    selective.att.g[is.nan(selective.att.g)] <- NA


    # get standard errors for each group specific ATT
    selective.se.inner <- lapply(glist, function(g) {
      whichg <- which( (group == g) & (g <= t) & (t <= (group + max_e)))  # added last condition to allow for limit on longest period included in att
      inf.func.g <- as.numeric(get_agg_inf_func(att=ATT,
                                                inf_func=inf_func_mat,
                                                whichones=whichg,
                                                weights_agg=pg[whichg]/sum(pg[whichg]),
                                                wif=NULL))
      se.g <- compute_se_agg(inf.func.g, dp)
      list(inf.func=inf.func.g, se=se.g)
    })

    # recover standard errors separately by group
    selective.se.g <- unlist(BMisc::getListElement(selective.se.inner, "se"))
    selective.se.g[selective.se.g <= sqrt(.Machine$double.eps)*10] <- NA

    # recover influence function separately by group
    selective.inf.func.g <- simplify2array(BMisc::getListElement(selective.se.inner, "inf.func"))

    # use multiplier bootstrap (across groups) to get critical value
    # for constructing uniform confidence bands
    selective.crit.val <- stats::qnorm(1 - alpha/2)
    # TODO; ALLOW FOR CONFIDENCE BAND AND MULTIPLIER BOOTSTRAP
    if ((!is.null(cband)) && (cband==TRUE)){
      # if(bstrap == FALSE){
      #   warning('Used bootstrap procedure to compute simultaneous confidence band.')
      # }
      # selective.crit.val <- mboot(selective.inf.func.g, dp)$crit.val
      #
      # if(is.na(selective.crit.val) | is.infinite(selective.crit.val)){
      #   warning('Simultaneous critival value is NA. This probably happened because we cannot compute t-statistic (std errors are NA). We then report pointwise conf. intervals.')
      #   selective.crit.val <- stats::qnorm(1 - alpha/2)
      #   dp$cband <- FALSE
      # }
      #
      # if(selective.crit.val < stats::qnorm(1 - alpha/2)){
      #   warning('Simultaneous conf. band is somehow smaller than pointwise one using normal approximation. Since this is unusual, we are reporting pointwise confidence intervals')
      #   selective.crit.val <- stats::qnorm(1 - alpha/2)
      #   dp$cband <- FALSE
      # }
      #
      # if(selective.crit.val >= 7){
      #   warning("Simultaneous critical value is arguably `too large' to be realible. This usually happens when number of observations per group is small and/or there is no much variation in outcomes.")
      # }
      warning("Confidence bands and multiplier bootstrap are not implemented yet. We are reporting pointwise confidence intervals.")
    }

    # get overall att under selective treatment timing
    # (here use pgg instead of pg because we can just look at each group)
    selective.att <- sum(selective.att.g * pgg)/sum(pgg)

    # account for having to estimate pgg in the influence function
    selective.wif <- get_weight_influence(keepers=1:length(glist),
                                          pg=pgg,
                                          weights=weights,
                                          G=G,
                                          group=group)

    # get overall influence function
    selective.inf.func <- get_agg_inf_func(att=selective.att.g,
                                           inf_func=selective.inf.func.g,
                                           whichones=(1:length(glist)),
                                           weights_agg=pgg/sum(pgg),
                                           wif=selective.wif)


    selective.inf.func <- as.numeric(selective.inf.func)
    # get overall standard error
    selective.se <- compute_se_agg(selective.inf.func, dp)
    if(!is.na(selective.se)){
      if((selective.se <= sqrt(.Machine$double.eps)*10)) selective.se <- NA
    }

    return(list(overall.att=selective.att,
                overall.se=selective.se,
                type=type,
                yname = yname,
                partition_name = partition_name,
                control_group = control_group,
                egt=orig_glist,
                att.egt=selective.att.g,
                se.egt=selective.se.g,
                crit.val.egt=selective.crit.val,
                inf.function = list(selective.inf.func.g = selective.inf.func.g,
                                    selective.inf.func = selective.inf.func)
    ))

  }

  #-----------------------------------------------------------------------------
  # calendar time effects
  #-----------------------------------------------------------------------------

  if (type == "calendar") {

    # drop time periods where no one is treated yet
    # (can't get treatment effects in those periods)
    minG <- min(group)
    calendar.tlist <- tlist[tlist>=minG]

    # calendar time specific atts
    calendar.att.t <- sapply(calendar.tlist, function(t1) {
      # look at post-treatment periods for group g
      whicht <- which( (t == t1) & (group <= t))
      attt <- ATT[whicht]
      pgt <- pg[whicht]/(sum(pg[whicht]))
      sum(pgt * attt)
    })

    # get standard errors and influence functions
    # for each time specific att
    calendar.se.inner <- lapply(calendar.tlist, function(t1) {
      whicht <- which( (t == t1) & (group <= t))
      pgt <- pg[whicht]/(sum(pg[whicht]))
      wif.t <- get_weight_influence(keepers=whicht,
                                    pg=pg,
                                    weights=weights,
                                    G=G,
                                    group=group)
      inf.func.t <- as.numeric(get_agg_inf_func(att=ATT,
                                                inf_func=inf_func_mat,
                                                whichones=whicht,
                                                weights_agg=pgt,
                                                wif=wif.t))
      se.t <- compute_se_agg(inf.func.t, dp)
      list(inf.func=inf.func.t, se=se.t)
    })

    # recover standard errors separately by time
    calendar.se.t <- unlist(BMisc::getListElement(calendar.se.inner, "se"))
    calendar.se.t[calendar.se.t <= sqrt(.Machine$double.eps)*10] <- NA
    # recover influence function separately by time
    calendar.inf.func.t <- simplify2array(BMisc::getListElement(calendar.se.inner, "inf.func"))

    # use multiplier boostrap (across groups) to get critical value
    # for constructing uniform confidence bands
    calendar.crit.val <-  stats::qnorm(1-alpha/2)
    # TODO; ALLOW FOR CONFIDENCE BAND AND MULTIPLIER BOOTSTRAP
    if ((!is.null(cband)) && (cband==TRUE)){
      # if(bstrap == FALSE){
      #   warning('Used bootstrap procedure to compute simultaneous confidence band')
      # }
      # calendar.crit.val <- mboot(calendar.inf.func.t, dp)$crit.val
      #
      # if(is.na(calendar.crit.val) | is.infinite(calendar.crit.val)){
      #   warning('Simultaneous critival value is NA. This probably happened because we cannot compute t-statistic (std errors are NA). We then report pointwise conf. intervals.')
      #   calendar.crit.val <- stats::qnorm(1 - alpha/2)
      #   dp$cband <- FALSE
      # }
      #
      # if(calendar.crit.val < stats::qnorm(1 - alpha/2)){
      #   warning('Simultaneous conf. band is somehow smaller than pointwise one using normal approximation. Since this is unusual, we are reporting pointwise confidence intervals')
      #   calendar.crit.val <- stats::qnorm(1 - alpha/2)
      #   dp$cband <- FALSE
      # }
      #
      # if(calendar.crit.val >= 7){
      #   warning("Simultaneous critical value is arguably `too large' to be realible. This usually happens when number of observations per group is small and/or there is no much variation in outcomes.")
      # }
      warning("Confidence bands and multiplier bootstrap are not implemented yet. We are reporting pointwise confidence intervals.")
    }

    # get overall att under calendar time effects
    # this is just average over all time periods
    calendar.att <- mean(calendar.att.t)

    # get overall influence function
    calendar.inf.func <- get_agg_inf_func(att=calendar.att.t,
                                          inf_func=calendar.inf.func.t,
                                          whichones=(1:length(calendar.tlist)),
                                          weights_agg=rep(1/length(calendar.tlist), length(calendar.tlist)),
                                          wif=NULL)
    calendar.inf.func <- as.numeric(calendar.inf.func)
    # get overall standard error
    calendar.se <- compute_se_agg(calendar.inf.func, dp)
    if(!is.na(calendar.se)){
      if (calendar.se <= sqrt(.Machine$double.eps)*10) calendar.se <- NA
    }
    return(list(overall.att=calendar.att,
                overall.se=calendar.se,
                type=type,
                yname = yname,
                partition_name = partition_name,
                control_group = control_group,
                egt=sapply(calendar.tlist,t2orig),
                att.egt=calendar.att.t,
                se.egt=calendar.se.t,
                crit.val.egt=calendar.crit.val,
                inf.function = list(calendar.inf.func.t = calendar.inf.func.t,
                                    calendar.inf.func = calendar.inf.func)
    ))
  }

  #-----------------------------------------------------------------------------
  # Compute the event-study estimators
  #-----------------------------------------------------------------------------

  if (type == "eventstudy") {

    # note: event times can be negative here.
    # note: event time = 0 corresponds to "on impact"
    # event times
    eseq <- unique(orig_periods - orig_group)
    eseq <- eseq[order(eseq)]

    # if the user specifies balance_e, then we are going to
    # drop some event times and some groups; if not, we just
    # keep everything (that is what this variable is for)
    include.balanced.gt <- rep(TRUE, length(orig_group))

    # if we balance the sample with respect to event time
    if (!is.null(balance_e)) {
      include.balanced.gt <- (t2orig(maxT) - orig_group >= balance_e)

      eseq <- unique(orig_periods[include.balanced.gt] - orig_group[include.balanced.gt])
      eseq <- eseq[order(eseq)]

      eseq <- eseq[ (eseq <= balance_e) & (eseq >= balance_e - t2orig(maxT) + t2orig(1))]

    }

    # only looks at some event times
    eseq <- eseq[ (eseq >= min_e) & (eseq <= max_e) ]

    # compute atts that are specific to each event time
    dynamic.att.e <- sapply(eseq, function(e) {
      # keep att(g,t) for the right g&t as well as ones that
      # are not trimmed out from balancing the sample
      whiche <- which( (orig_periods - orig_group == e) & (include.balanced.gt) )
      atte <- ATT[whiche]
      pge <- pg[whiche]/(sum(pg[whiche]))
      sum(atte*pge)
    })

    # compute standard errors for dynamic effects
    dynamic.se.inner <- lapply(eseq, function(e) {
      whiche <- which( (orig_periods - orig_group == e) & (include.balanced.gt) )
      pge <- pg[whiche]/(sum(pg[whiche]))
      wif.e <- get_weight_influence(whiche, pg, weights, G, group)
      inf.func.e <- as.numeric(get_agg_inf_func(att=ATT,
                                                inf_func=inf_func_mat,
                                                whichones=whiche,
                                                weights_agg=pge,
                                                wif=wif.e))
      se.e <- compute_se_agg(inf.func.e, dp)
      list(inf.func=inf.func.e, se=se.e)
    })

    dynamic.se.e <- unlist(BMisc::getListElement(dynamic.se.inner, "se"))
    dynamic.se.e[dynamic.se.e <= sqrt(.Machine$double.eps)*10] <- NA

    dynamic.inf.func.e <- simplify2array(BMisc::getListElement(dynamic.se.inner, "inf.func"))

    dynamic.crit.val <- stats::qnorm(1 - alpha/2)
    # TODO; ALLOW FOR CONFIDENCE BAND AND MULTIPLIER BOOTSTRAP
    if ((!is.null(cband)) && (cband==TRUE)){
      # if(bstrap == FALSE){
      #   warning('Used bootstrap procedure to compute simultaneous confidence band')
      # }
      # calendar.crit.val <- mboot(calendar.inf.func.t, dp)$crit.val
      #
      # if(is.na(calendar.crit.val) | is.infinite(calendar.crit.val)){
      #   warning('Simultaneous critival value is NA. This probably happened because we cannot compute t-statistic (std errors are NA). We then report pointwise conf. intervals.')
      #   calendar.crit.val <- stats::qnorm(1 - alpha/2)
      #   dp$cband <- FALSE
      # }
      #
      # if(calendar.crit.val < stats::qnorm(1 - alpha/2)){
      #   warning('Simultaneous conf. band is somehow smaller than pointwise one using normal approximation. Since this is unusual, we are reporting pointwise confidence intervals')
      #   calendar.crit.val <- stats::qnorm(1 - alpha/2)
      #   dp$cband <- FALSE
      # }
      #
      # if(calendar.crit.val >= 7){
      #   warning("Simultaneous critical value is arguably `too large' to be realible. This usually happens when number of observations per group is small and/or there is no much variation in outcomes.")
      # }
      warning("Confidence bands and multiplier bootstrap are not implemented yet. We are reporting pointwise confidence intervals.")
    }

    # get overall average treatment effect
    # by averaging over positive dynamics
    epos <- eseq >= 0
    dynamic.att <- mean(dynamic.att.e[epos])
    dynamic.inf.func <- get_agg_inf_func(att=dynamic.att.e[epos],
                                         inf_func=as.matrix(dynamic.inf.func.e[,epos]),
                                         whichones=(1:sum(epos)),
                                         weights_agg=(rep(1/sum(epos), sum(epos))),
                                         wif=NULL)

    dynamic.inf.func <- as.numeric(dynamic.inf.func)
    dynamic.se <- compute_se_agg(dynamic.inf.func, dp)
    if(!is.na(dynamic.se)){
      if (dynamic.se <= sqrt(.Machine$double.eps)*10) dynamic.se <- NA
    }

    return(list(overall.att=dynamic.att,
                overall.se=dynamic.se,
                type=type,
                yname = yname,
                partition_name = partition_name,
                control_group = control_group,
                egt=eseq,
                att.egt=dynamic.att.e,
                se.egt=dynamic.se.e,
                crit.val.egt=dynamic.crit.val,
                inf.function = list(dynamic.inf.func.e = dynamic.inf.func.e,
                                    dynamic.inf.func = dynamic.inf.func),
                min_e=min_e,
                max_e=max_e,
                balance_e=balance_e
    ))
  }

}
