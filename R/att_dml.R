#' DML DDD estimator for ATT, with panel data and 2 periods
#'
#' This function implements a double machine learning estimator for assessing the average
#' treatment effect on the treated (ATT) using a triple differences (DDD) approach
#' in panel data settings across two time periods. The function takes preprocessed
#' data structured specifically for this analysis.
#'
#' @import stats
#' @import mlr3
#' @import mlr3learners
#' @param did_preprocessed A list containing preprocessed data and specifications for the DDD estimation.
#'        Expected elements include:
#'        - `preprocessed_data`: A data table containing the data with variables needed for the analysis.
#'        - `xformula`: The formula for the covariates to be included in the model. It should be of the form \code{~ x1 + x2}.
#'        Default is \code{xformla = ~1} (no covariates).
#'        - `learners`: A list of learners to be used in the estimation. It should be a list of two elements,
#'        the first element being the learner for the propensity score and the second element being the learner for the outcome regression.
#'        Default is \code{NULL}, then OLS and MLE Logit is used to estimate nuisances parameters. If \code{est_method = "dml"}, user have to specify \code{learners}.
#'        - `n_folds`: The number of folds to be used in the cross-fitting. Default is \code{NULL}. If \code{est_method = "dml"}, user have to specify \code{n_folds}.
#'        - `alpha`: The level of significance for the confidence intervals. Default is \code{0.05}.
#'        - `subgroup_counts`: A matrix containing the number of observations in each subgroup.
#'
#' @keywords internal
#' @return A list with the estimated ATT, standard error, upper and lower confidence intervals, and n_folds.
#' @noRd
NULL
# ------------------------------------------------------------------------------
att_dml <- function(did_preprocessed) {

  # get parameters needed
  data <- copy(did_preprocessed$preprocessed_data)
  xformula <- did_preprocessed$xformula
  learners <- did_preprocessed$learners
  n_folds <- did_preprocessed$n_folds
  alpha <- did_preprocessed$alpha
  subgroup_counts <- did_preprocessed$subgroup_counts
  inffunc <- did_preprocessed$inffunc # flag to return influence function
  ## ---------- Learners ---------------------------
  ml_pa = learners$ml_pa # classifier for propensity score
  ml_md = learners$ml_md # regression for outcome model

  ## ---------- create relevant variables -----------------
  data <- get_wide_data(data)

  # Adding treatment variable for model P(1{S=2, Q=1}|X)
  data[, "D" := ifelse(data$subgroup == 4, 1, 0)]
  # Adding outcome variable for model E[deltaY | X]
  data[, deltaY:= y1 - y0]

  ## ---------- Global folds adapted for DDD -----------
  global_folds <- make_stratified_folds(data = data,
                                        strat_col = "subgroup",
                                        n_folds = n_folds,
                                        seed = 1234)

  fold_test <- global_folds$fold_test
  fold_train <- global_folds$fold_train

  # --------------------------------------------------------------------
  # Compute ATT
  # --------------------------------------------------------------------

  # Compute DML Triple Difference Estimator
  dml_att_scores <- compute_dml_nuisances(data = data,
                                          xformula = xformula,
                                          ml_pa = ml_pa,
                                          ml_md = ml_md,
                                          fold_test = fold_test,
                                          fold_train = fold_train,
                                          n_folds = n_folds)

  # Compute ATT
  att_dml = dml_att_scores$att[["3"]] + dml_att_scores$att[["2"]] - dml_att_scores$att[["1"]]

  # Inference
  se_inf_dml_scores <- compute_se_dml(dml_att_scores$influence_matrix)
  se_dml <- se_inf_dml_scores$se
  # Return null if inffunc is FALSE
  if (inffunc == FALSE){
    inf_func <- NULL
  } else {
    inf_func <- se_inf_dml_scores$inf_func
  }


  # estimate upper bound at 1-alpha% confidence level
  ci_upper <- att_dml + qnorm(1-alpha/2) * se_dml
  # estimate lower bound at 95% confidence level
  ci_lower <- att_dml - qnorm(1-alpha/2) * se_dml

  # ------------------------------------------------------------------------------
  # Return results
  # ------------------------------------------------------------------------------

  ret <- (list(ATT = att_dml,
               se = se_dml,
               uci = ci_upper,
               lci = ci_lower,
               inf_func = inf_func,
               subgroup_counts = subgroup_counts
               ))

  return(ret)

}
