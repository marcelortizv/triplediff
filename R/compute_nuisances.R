#' Function to compute nuisances parameters and did
#' @import stats
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
  if (any(propensity_scores < 0.01) | any(propensity_scores > 0.99)) {
    warning(paste("Propensity scores for subgroup", condition_subgroup,
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
  i.weights = control_data[control_data$post == 0, weights]

  # get coefficients of outcome regression model for the control group (subgroup != 4)
  # Attempt to compute reg.coeff using speedglm
  try_speedglm <- tryCatch({
    reg.coeff <- stats::coef(speedglm::speedglm.wfit(y = deltaY_control,
                                                     X = as.matrix(cov_control),
                                                     intercept = FALSE,
                                                     weights = i.weights))
  }, error = function(e) {
    # If error, attempt to compute reg.coeff using lm.wfit
    reg.coeff <- tryCatch({
      stats::coef(stats::lm.wfit(x = as.matrix(cov_control),
                                 y = deltaY_control,
                                 w = i.weights))
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
  cov <- stats::model.matrix(xformula, data = condition_data[condition_data$post == 0])
  # compute OR delta
  or_delta = as.vector(tcrossprod(reg.coeff, as.matrix(cov)))
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
  i.weights = condition_data[condition_data$post == 0, weights]

  ################################
  # Get doubly-robust estimation #
  ################################

  w.treat = i.weights * PA4
  w.control = (i.weights * pscore * PAa) / (1 - pscore)
  riesz.treat = w.treat * (deltaY - or_delta)
  riesz.control = w.control * (deltaY - or_delta)
  att.treat = mean(riesz.treat, na.rm = TRUE) / mean(w.treat, na.rm = TRUE)
  att.control = mean(riesz.control, na.rm = TRUE) / mean(w.control, na.rm = TRUE)

  dr.att = att.treat - att.control

  ##########################
  # Get influence function #
  ##########################

  # Influence function related to the estimation of pscores
  cov = stats::model.matrix(xformula, data = condition_data)
  M2 <- base::colMeans(w.control * (deltaY - or_delta - att.control) * cov, na.rm = TRUE) # reg_adjust = deltaY - m_delta(x)

  score_ps <- i.weights * (PA4 - pscore) * cov
  # asymptotic linear representation of logit's beta
  score_ps_no_na <- na.omit(score_ps) # Exclude rows with NA values
  asy.lin.rep.ps <- score_ps_no_na %*% hessian
  inf.control.pscore <- asy.lin.rep.ps %*% as.matrix(M2)

  # Influence function related to the estimation of pscores

  M1 <- base::colMeans(w.treat * cov, na.rm = TRUE)
  M3 <- base::colMeans(w.control * cov, na.rm = TRUE)

  # Influence function related to the estimation of regression model
  or_x <- i.weights * PAa * cov
  or_ex <- i.weights * PAa * (deltaY - or_delta) * cov
  XpX <- crossprod(or_x, cov)/nrow(condition_data)

  #asymptotic linear representation of the beta
  asy.linear.or <- t(solve(XpX, t(or_ex)))

  #or for treat
  inf.treat.or <- -asy.linear.or %*% M1
  #or for control
  inf.cont.or <- -asy.linear.or %*% M3

  # Influence function from did
  inf.control.did <- riesz.control - w.control*att.control
  inf.treat.did <- riesz.treat - w.treat*att.treat

  # Influence function for the ATT
  inf.control <- (inf.control.did + inf.control.pscore + inf.cont.or)/mean(w.control)
  inf.treat <- (inf.treat.did + inf.treat.or)/mean(w.treat)
  # putting all together
  inf.func <- inf.treat - inf.control

  # fill zeros in influence function for observation outside the subgroup analyzed
  data[, "inf.func.result" := numeric(.N)]
  data[data$subgroup %in% c(condition_subgroup, 4), "inf.func.result" := inf.func]

  return(list(dr.att = dr.att, inf.func = data[["inf.func.result"]]))
}

### FUNCTIONS FOR DML ###

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
compute_scores <- function(ml_pa, ml_md){
  #' Function to compute scores based on first step nuisances predictions
  #' @param ml_pa: MLR3 object for propensity score estimation
  #' @param ml_md: MLR3 object for regression adjustment estimation
  #' @return: A list of scores and id for each k fold

  # Warning for overlap condition
  if (any(ml_pa$prediction()$prob[,2] < 0.01) | any(ml_pa$prediction()$prob[,2] > 0.99)) {
    warning(paste("Propensity scores for subgroup", condition_subgroup,
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
  dml_scores <- compute_scores(ml_pa = ml_classif, ml_md = ml_regr)

  # Get a list of dmlddd estimator + scores, for each k-fold
  ddd_dml <- lapply(dml_scores, compute_dml)

  return(ddd_dml)
}

# Function to compute se for DML ATT
compute_se_dml <- function(dml_scores1, dml_scores2, dml_scores3) {
  scores <- c(dml_scores1, dml_scores2, dml_scores3)
  N = length(scores)
  sigma2 = stats::var(scores)
  se = stats::sqrt(sigma2 / N)
  return(se)
}
