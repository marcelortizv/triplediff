#' Doubly robust DDD estimator for ATT, with repeated cross-section data and 2 periods
#'
#' This function implements a doubly robust estimator for assessing the average
#' treatment effect on the treated (ATT) using a triple differences (DDD) approach
#' in repeated cross-section data settings across two time periods. The function takes preprocessed
#' data structured specifically for this analysis.
#'
#' @import stats
#' @param did_preprocessed A list containing preprocessed data and specifications for the DDD estimation.
#'        Expected elements include:
#'        - `preprocessed_data`: A data table containing the data with variables needed for the analysis.
#'        - `est_method`: The estimation method to be used. Default is \code{est_method = "dr"}.
#'        - `xformula`: The formula for the covariates to be included in the model. It should be of the form \code{~ x1 + x2}.
#'        Default is \code{xformla = ~1} (no covariates).
#'        - `boot`: Logical. If \code{TRUE}, the function use the multiplier bootstrap to compute standard errors. Default is \code{FALSE}.
#'        - `nboot`: The number of bootstrap samples to be used. Default is \code{NULL}. If \code{boot = TRUE}, the default is \code{nboot = 999}.
#'        - `subgroup_counts`: A matrix containing the number of observations in each subgroup.
#'        - `alpha` The level of significance for the confidence intervals.  Default is \code{0.05}.
#'        - `inffunc`: Logical. If \code{TRUE}, the function returns the influence function. Default is \code{FALSE}.
#'        - `use_parallel`: Boolean of whether or not to use parallel processing in the multiplier bootstrap, default is \code{use_parallel=FALSE}
#'        - `cores`: the number of cores to use with parallel processing, default is \code{cores=1}
#'        - `cband`: Boolean of whether or not to compute simultaneous confidence bands, default is \code{cband=FALSE}
#'
#' @keywords internal
#' @return A list with the estimated ATT, standard error, upper and lower confidence intervals, and influence function.
#' @export
att_dr_rc <- function(did_preprocessed) {

  data <- did_preprocessed$preprocessed_data
  est_method <- did_preprocessed$est_method
  xformula <- did_preprocessed$xformula
  boot <- did_preprocessed$boot
  nboot <- did_preprocessed$nboot
  alpha <- did_preprocessed$alpha
  use_parallel <- did_preprocessed$use_parallel # to perform bootstrap
  cores <- did_preprocessed$cores # to perform bootstrap
  cband <- did_preprocessed$cband # to perform bootstrap + simult. conf. band
  inffunc <- did_preprocessed$inffunc # flag to return influence function
  subgroup_counts <- did_preprocessed$subgroup_counts

  # --------------------------------------------------------------------
  # Compute ATT
  # --------------------------------------------------------------------

  # Pre-compute propensity scores for each subgroup
  if (est_method == "reg"){
    # set pscores properly such that weights for control are zero
    pscores <- lapply(c(3, 2, 1), function(condition_subgroup) {
      compute_pscore_null_rc(data, condition_subgroup)
    })
  } else {
    pscores <- lapply(c(3, 2, 1), function(condition_subgroup) {
      compute_pscore_rc(data, condition_subgroup, xformula)
    })
  }

  if (est_method == "ipw"){
    # set or_delta equal to zero
    reg_adjust <- lapply(c(3, 2, 1), function(condition_subgroup) {
      compute_outcome_regression_null_rc(data, condition_subgroup)
    })
  } else {
    # Pre-compute the regression adjustment for each subgroup
    reg_adjust <- lapply(c(3, 2, 1), function(condition_subgroup) {
      compute_outcome_regression_rc(data, condition_subgroup, xformula)
    })
  }

  # Compute Doubly Robust Triple Difference Estimator
  dr_att_inf_func_3 <- compute_did_rc(data, condition_subgroup = 3, pscores = pscores, reg_adjustment = reg_adjust, xformula = xformula, est_method = est_method) # S=g, Q=1 vs. S=g, Q=0
  dr_att_inf_func_2 <- compute_did_rc(data, condition_subgroup = 2, pscores = pscores, reg_adjustment = reg_adjust, xformula = xformula, est_method = est_method) # S=g, Q=1 vs. S=\infty, Q=1
  dr_att_inf_func_1 <- compute_did_rc(data, condition_subgroup = 1, pscores = pscores, reg_adjustment = reg_adjust, xformula = xformula, est_method = est_method) # S=g, Q=1 vs. S=\infty, Q=0

  dr_ddd <- dr_att_inf_func_3$dr_att + dr_att_inf_func_2$dr_att - dr_att_inf_func_1$dr_att
  n <- data[, .N]
  n3 <- subgroup_counts$count[subgroup_counts$subgroup == 3] + subgroup_counts$count[subgroup_counts$subgroup == 4]
  n2 <- subgroup_counts$count[subgroup_counts$subgroup == 2] + subgroup_counts$count[subgroup_counts$subgroup == 4]
  n1 <- subgroup_counts$count[subgroup_counts$subgroup == 1] + subgroup_counts$count[subgroup_counts$subgroup == 4]
  w3 <- n/n3
  w2 <- n/n2
  w1 <- n/n1
  # rescaling influence function
  inf_func = w3*dr_att_inf_func_3$inf_func + w2*dr_att_inf_func_2$inf_func - w1*dr_att_inf_func_1$inf_func

  #-----------------------------------------------------------------------------
  # compute confidence intervals / bands
  #-----------------------------------------------------------------------------

  if (boot){
    # perform multiplier bootstrap
    boot_result <- mboot(inf_func, did_preprocessed=did_preprocessed, use_parallel=use_parallel, cores=cores)
    se_ddd <- boot_result$se # save bootstrap standard error
    bT <- boot_result$bT # save sup-t confidence band
    if (cband){

      # get critical value to compute uniform confidence bands
      cv <- boot_result$unif_crit_val
      if(cv >= 7){
        warning("Simultaneous critical value is arguably `too large' to be reliable. This usually happens when number of observations per group is small and/or there is no much variation in outcomes.")
      }

    } else {
      # use regular critical value
      cv <- qnorm(1-alpha/2)
    }
    # Computing uniform confidence bands
    # Estimate of upper boundary of 1-alpha% confidence band
    ci_upper <- dr_ddd + cv * se_ddd
    # Estimate of lower boundary of 1-alpha% confidence band
    ci_lower <- dr_ddd - cv * se_ddd
  } else {
    # compute point-wise confidence intervals
    se_ddd <- stats::sd(inf_func)/sqrt(n)
    bT <- NULL
    # use regular critical value
    cv <- qnorm(1-alpha/2)
    # estimate upper bound at 1 - alpha% confidence level
    ci_upper <- dr_ddd + cv * se_ddd
    # estimate lower bound at 95% confidence level
    ci_lower <- dr_ddd - cv * se_ddd
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
              bT = bT,
              inf_func = inf_func,
              subgroup_counts = subgroup_counts
              ))

  return(ret)

}

