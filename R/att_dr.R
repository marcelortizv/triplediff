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
#'        - `preprocessed_data`: A data table containing the data with variables needed for the analysis.
#'        - `xformula`: The formula for the covariates to be included in the model. It should be of the form \code{~ x1 + x2}.
#'        Default is \code{xformla = ~1} (no covariates).
#'        - `boot`: Logical. If \code{TRUE}, the function use the multiplier bootstrap to compute standard errors. Default is \code{FALSE}.
#'        - `nboot`: The number of bootstrap samples to be used. Default is \code{NULL}. If \code{boot = TRUE}, the default is \code{nboot = 999}.
#'        - `subgroup_counts`: A matrix containing the number of observations in each subgroup.
#'        - `inffunc`: Logical. If \code{TRUE}, the function returns the influence function. Default is \code{FALSE}.
#'
#' @keywords internal
#' @return A list with the estimated ATT, standard error, upper and lower confidence intervals, and influence function.
#' @noRd
NULL
# ------------------------------------------------------------------------------

att_dr <- function(did_preprocessed) {

  data <- did_preprocessed$preprocessed_data
  xformula <- did_preprocessed$xformula
  boot <- did_preprocessed$boot
  nboot <- did_preprocessed$nboot
  cluster <- did_preprocessed$cluster
  alpha <- did_preprocessed$alpha
  cband <- did_preprocessed$cband
  use_parallel <- did_preprocessed$use_parallel
  cores <- did_preprocessed$cores
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

  dr_ddd <- dr_att_inf_func_3$dr_att - dr_att_inf_func_2$dr_att + dr_att_inf_func_1$dr_att
  n <- data[, .N/2]
  n3 <- subgroup_counts$V1[1] + subgroup_counts$V1[2]
  n2 <- subgroup_counts$V1[1] + subgroup_counts$V1[3]
  n1 <- subgroup_counts$V1[1] + subgroup_counts$V1[4]
  w3 <- n/n3
  w2 <- n/n2
  w1 <- n/n1
  # rescaling influence function
  inf_func = w3*dr_att_inf_func_3$inf_func - w2*dr_att_inf_func_2$inf_fun + w3*dr_att_inf_func_1$inf_func

  # ---------------------------------------------------------------------
  # Compute Variance
  # ---------------------------------------------------------------------

  if (boot == TRUE){

    # perform multiplier bootstrap
    boot_result <- mboot(inf_func, did_preprocessed=did_preprocessed, use_parallel=use_parallel, cores=cores)
    se_ddd <- boot_result$se
    cv <- boot_result$crit_val
    if(cv >= 7){
      warning("Simultaneous critical value is arguably `too large' to be realible. This usually happens when number of observations per group is small and/or there is no much variation in outcomes.")
    }
    # Computing uniform confidence band
    # Estimate of upper boundary of 1-alpha% confidence band
    ci_upper <- dr_ddd + cv * se_ddd
    # Estimate of lower boundary of 1-alpha% confidence band
    ci_lower <- dr_ddd - cv * se_ddd

  } else {
    se_ddd <- stats::sd(inf_func)/sqrt(n)
    # estimate upper bound at 1 - alpha% confidence level
    ci_upper <- dr_ddd + qnorm(1-alpha/2) * se_ddd
    # estimate lower bound at 95% confidence level
    ci_lower <- dr_ddd - qnorm(1-alpha/2) * se_ddd
  }

  # ------------------------------------------------------------------------------
  # Return results
  # ------------------------------------------------------------------------------

  # Return null if inffunc is FALSE
  if (inffunc == FALSE){
    inf_func <- NULL
  }

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
