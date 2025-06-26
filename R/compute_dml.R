#' Utility functions to compute nuisances parameters based on DML.
#' @importFrom stats update
#' @import parglm
#' @import speedglm
#' @import data.table
#' @import mlr3
#' @import mlr3learners
#' @import mlr3tuning
#' @noRd
#--------------------------------------------------

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

# ---------------------------- #
# FUNCTIONS FOR DML
# ---------------------------- #

# Function to get components of a linear score
#' Compute the score elements for the ATT
#' @param y Outcome variable
#' @param d Treatment variable
#' @param p_hat Propensity score
#' @param m_hat Outcome regression
#' @return A list containing the score elements psi_a and psi_b
#' @keywords internal
get_score_elements <- function(y, d, p_hat, m_hat){

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
#' Function to compute scores based on first step nuisances predictions
#' @param ml_pa MLR3 object for propensity score estimation
#' @param ml_md MLR3 object for regression adjustment estimation
#' @param condition_subgroup Subgroup to analyze
#' @return: A list of scores and id for each k fold
#' @keywords internal
compute_scores <- function(ml_pa, ml_md, condition_subgroup){


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
#' Compute the propensity score for the given condition subgroup
#' @param data A data.table containing the data processed
#' @param condition_subgroup The condition subgroup for which to compute the propensity score
#' @param xformula The formula for the propensity score model
#' @param ml_pa The machine learning algorithm to use for the propensity score model
#' @param ml_md The machine learning algorithm to use for the regression adjustment model
#' @param n_folds The number of folds for cross-fitting
#' @return A list containing ids, estimator for each k fold, and scores
#' @keywords internal
compute_dml_nuisances <- function(data, condition_subgroup, xformula, ml_pa, ml_md, n_folds) {

  set.seed(123)

  formula_pscore <- get_formula_pscore(xformula, weights = TRUE)
  formula_reg <- get_formula_reg(xformula, weights = TRUE)

  # Subset data for subgroup == 4 or the given condition_subgroup
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]

  # drop period column
  condition_data <- condition_data[, period := NULL]

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
