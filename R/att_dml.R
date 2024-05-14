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
#'        - `preprocessed_data`: the data table containing the variables needed for the analysis.
#'        - `xformula`: formula object specifying the model for the nuisance functions.
#'        - `estMethod`: the estimation method to use.
#'        - `learners`: specified machine learning methods for nuisance function estimation.
#'        - `n_folds`: the number of folds for cross-fitting
#'
#' @keywords internal
#' @return A list with the estimated ATT, standard error, upper and lower confidence intervals, and n_folds.
#' @noRd
NULL
# ------------------------------------------------------------------------------

att_dml <- function(did_preprocessed) {

  data <- did_preprocessed$preprocessed_data
  xformula <- did_preprocessed$xformula
  estMethod <- did_preprocessed$estMethod
  learners <- did_preprocessed$learners
  n_folds <- did_preprocessed$n_folds

  # --------------------------------------------------------------------
  # Compute ATT
  # --------------------------------------------------------------------

  # Compute DML Triple Difference Estimator
  dml.att.scores_3 <- compute_dml_nuisances(data, condition_subgroup = 3,
                                            xformula,
                                            ml_pa = learners$ml_pa,
                                            ml_md = learners$ml_md,
                                            n_folds = n_folds)

  dml.att.scores_2 <- compute_dml_nuisances(data, condition_subgroup = 2,
                                            xformula,
                                            ml_pa = learners$ml_pa,
                                            ml_md = learners$ml_md,
                                            n_folds = n_folds)

  dml.att.scores_1 <- compute_dml_nuisances(data, condition_subgroup = 1,
                                            xformula,
                                            ml_pa = learners$ml_pa,
                                            ml_md = learners$ml_md,
                                            n_folds = n_folds)

  # Compute ATT
  att_dml = mean_ddd_k(dml.att.scores_3) - mean_ddd_k(dml.att.scores_2) + mean_ddd_k(dml.att.scores_1)

  # Inference
  scores_3 <- get_long_scores(dml.att.scores_3)
  scores_2 <- get_long_scores(dml.att.scores_2)
  scores_1 <- get_long_scores(dml.att.scores_1)

  se_dml <- compute_se_dml(scores_3, scores_2, scores_1)

  # estimate upper bound at 95% confidence level
  ci_upper <- att_dml + 1.96 * se_dml
  # estimate lower bound at 95% confidence level
  ci_lower <- att_dml - 1.96 * se_dml

  # ------------------------------------------------------------------------------
  # Return results
  # ------------------------------------------------------------------------------

  ret <- (list(ATT = att_dml,
               se = se_dml,
               uci = ci_upper,
               lci = ci_lower
               ))

  return(ret)

}
