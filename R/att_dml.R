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
  data <- did_preprocessed$preprocessed_data
  xformula <- did_preprocessed$xformula
  learners <- did_preprocessed$learners
  n_folds <- did_preprocessed$n_folds
  alpha <- did_preprocessed$alpha
  subgroup_counts <- did_preprocessed$subgroup_counts

  # --------------------------------------------------------------------
  # Compute ATT
  # --------------------------------------------------------------------

  # Compute DML Triple Difference Estimator
  dml_att_scores_3 <- compute_dml_nuisances(data, condition_subgroup = 3,
                                            xformula,
                                            ml_pa = learners$ml_pa,
                                            ml_md = learners$ml_md,
                                            n_folds = n_folds)

  dml_att_scores_2 <- compute_dml_nuisances(data, condition_subgroup = 2,
                                            xformula,
                                            ml_pa = learners$ml_pa,
                                            ml_md = learners$ml_md,
                                            n_folds = n_folds)

  dml_att_scores_1 <- compute_dml_nuisances(data, condition_subgroup = 1,
                                            xformula,
                                            ml_pa = learners$ml_pa,
                                            ml_md = learners$ml_md,
                                            n_folds = n_folds)

  # Compute ATT
  att_dml = mean_ddd_k(dml_att_scores_3) + mean_ddd_k(dml_att_scores_2) - mean_ddd_k(dml_att_scores_1)

  # Inference
  scores_3 <- get_long_scores(dml_att_scores_3)
  scores_2 <- get_long_scores(dml_att_scores_2)
  scores_1 <- get_long_scores(dml_att_scores_1)

  se_dml <- compute_se_dml(scores_3, scores_2, scores_1)

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
               subgroup_counts = subgroup_counts
               ))

  return(ret)

}
