#' Doubly robust DDD estimator for ATT, with panel data and 2 periods
#'
#' This function implements a doubly robust estimator for assessing the average
#' treatment effect on the treated (ATT) using a triple differences (DDD) approach
#' in panel data settings across two time periods. The function takes preprocessed
#' data structured specifically for this analysis.
#'
#' @import stats
#' @param did_preprocessed A list containing preprocessed data and specifications for the DDD estimation.
#'        Expected elements include:
#'        - `preprocessed_data`: the data table containing the variables needed for the analysis.
#'        - `xformula`: formula object specifying the model for the nuisance functions.
#'        - `est_method`: the estimation method to use.
#'        - `learners`: specified machine learning methods for nuisance function estimation.
#'        - `boot`: logical indicating whether bootstrap should be used for inference.
#'        - `boot_type`: type of bootstrap to perform.
#'        - `nboot`: number of bootstrap replicates.
#'        - `inffunc`: logical indicating whether return or not influence function.
#'
#' @keywords internal
#' @return A list with the estimated ATT, standard error, upper and lower confidence intervals, and influence function.
#' @noRd
NULL
# ------------------------------------------------------------------------------

att_dr <- function(did_preprocessed) {

  data <- did_preprocessed$preprocessed_data
  xformula <- did_preprocessed$xformula
  est_method <- did_preprocessed$est_method
  learners <- did_preprocessed$learners
  boot <- did_preprocessed$boot
  boot_type <- did_preprocessed$boot_type
  nboot <- did_preprocessed$nboot
  inffunc <- did_preprocessed$inffunc
  subgroup_counts <- did_preprocessed$subgroup_counts

  # --------------------------------------------------------------------
  # Compute ATT
  # --------------------------------------------------------------------

  # Pre-compute propensity scores for each subgroup
  pscores <- lapply(c(3, 2, 1), function(condition_subgroup) {
    compute_pscore(data, condition_subgroup, xformula)
  })


  # Pre-compute the regression adjustment for each subgroup
  reg_adjust <- lapply(c(3, 2, 1), function(condition_subgroup) {
    compute_outcome_regression(data, condition_subgroup, xformula)
  })


  # Compute Doubly Robust Triple Difference Estimator
  dr_att_inf_func_3 <- compute_did(data, condition_subgroup = 3, pscores, reg_adjust, xformula) # S=2, L=A
  dr_att_inf_func_2 <- compute_did(data, condition_subgroup = 2, pscores, reg_adjust, xformula) # S=\infty, L=B
  dr_att_inf_func_1 <- compute_did(data, condition_subgroup = 1, pscores, reg_adjust, xformula) # S=\infty, L=A

  dr_ddd = dr_att_inf_func_3$dr_att - dr_att_inf_func_2$dr_att + dr_att_inf_func_1$dr_att
  inf_func = dr_att_inf_func_3$inf_func - dr_att_inf_func_2$inf_fun + dr_att_inf_func_1$inf_func

  # ---------------------------------------------------------------------
  # Compute Variance
  # ---------------------------------------------------------------------

  if (boot == TRUE){
    if (is.null(nboot)){
      nboot <- 999
    }
    if (boot_type == "multiplier"){
      # perform multiplier bootstrap
      inf_boot <- mboot_did(inf_func, nboot)
      # get bootstrap std errors based on IQR
      se_ddd <- stats::IQR(inf_boot) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmetric critical values
      cv <- stats::quantile(abs(inf_boot/se_ddd), probs = 0.95)
      # Estimate of upper boundary of 95% CI
      ci_upper <- dr_ddd + cv * se_ddd
      # Estimate of lower boundary of 95% CI
      ci_lower <- dr_ddd - cv * se_ddd
    } else {
      stop("Bootstrapping type other than multiplier is currently not supported.")
    }

  } else {
    se_ddd <- stats::sd(inf_func)/sqrt(nrow(data))
    # estimate upper bound at 95% confidence level
    ci_upper <- dr_ddd + 1.96 * se_ddd
    # estimate lower bound at 95% confidence level
    ci_lower <- dr_ddd - 1.96 * se_ddd
  }

  # ------------------------------------------------------------------------------
  # Return results
  # ------------------------------------------------------------------------------

  ret <- (list(ATT = dr_ddd,
              se = se_ddd,
              uci = ci_upper,
              lci = ci_lower,
              nboot = nboot,
              inf_func = inf_func,
              subgroup_counts = subgroup_counts
              ))

  return(ret)

}
