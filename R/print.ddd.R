#' @title Print
#'
#' @description Prints a ddd Object
#'
#' @param x A ddd object
#' @param alpha The significance level for confidence intervals (optional)
#' @param ... Other params (required as generic function, but not used)
#' @importFrom utils write.table
#' @export
#' @noRd

print.ddd <- function(x, alpha = NULL, ...) {
  # Update confidence interval if alpha is provided and different from the one used in estimation
  if (!is.null(alpha) && alpha != x$argu$alpha) {
    x$argu$alpha <- alpha
    if (x$argu$cband){
      new_cv <- quantile(x$bT, 1 - alpha, type=1, na.rm = T)
    } else {
      new_cv <- qnorm(1 - alpha / 2)
    }

    x$lci <- x$ATT - new_cv * x$se
    x$uci <- x$ATT + new_cv * x$se
  }

  # Methods used
  if (x$argu$multiple_periods == FALSE){
    if (x$argu$est_method[1] == 'dr') {
      est_method1 <- "DR-DDD estimation for the ATT: \n"
      est_method2 <- "Outcome Regression estimated using: OLS"
      est_method3 <- "Propensity score estimated using: Maximum Likelihood"
    } else if (x$argu$est_method[1] == 'ipw') {
      est_method1 <- "IPW-DDD estimation for the ATT: \n"
      est_method2 <- "Outcome Regression estimated using: N/A"
      est_method3 <- "Propensity score estimated using: Maximum Likelihood"
    } else if (x$argu$est_method[1] == 'reg') {
      est_method1 <- "REG-DDD estimation for the ATT: \n"
      est_method2 <- "Outcome Regression estimated using: OLS"
      est_method3 <- "Propensity score estimated using: N/A"
    } #else if (x$argu$est_method[1] == 'dml') {
    #   est_method1 <- "DML-DDD estimation for the ATT: \n"
    #   est_method2 <- paste0("Outcome Regression estimated using: ", x$argu$learners$ml_md$label)
    #   est_method3 <- paste0("Propensity score estimated using: ", x$argu$learners$ml_pa$label)
    # }

    # Front-end Summary Table
    lev_conf <- paste0(round(100 * (1 - x$argu$alpha), digits = 2), "% ")
    band_method <- ifelse(x$argu$cband, "Simult.", "Ptwise.")
    interval_text <- paste0("[", lev_conf, band_method)
    #header <- c("ATT", "Std. Error", "t value", "Pr(>|t|)", interval_text, "Conf. Band]", "")
    header <- c("  ATT", "  Std. Error", "   Pr(>|t|)", interval_text, "Conf. Band]", "")
    sig <- (x$uci < 0) | (x$lci > 0)
    sig[is.na(sig)] <- FALSE
    sig_text <- ifelse(sig, "*", "")

    body <- cbind(
      sprintf("%10.4f", x$ATT),
      sprintf("%10.4f", x$se),
      #sprintf("%10.4f", x$ATT / x$se),
      sprintf("%10.4f", 2 * (stats::pnorm(-abs(x$ATT / x$se)))),
      sprintf("%10.4f", x$lci),
      sprintf("%10.4f", x$uci),
      sig_text
    )
    colnames(body) <- header
  }

  if (x$argu$multiple_periods == TRUE){
    if (x$argu$est_method[1] == 'dr') {
      est_method1 <- "DR-DDD estimation for the ATT(g,t): \n"
      est_method2 <- "Outcome Regression estimated using: OLS"
      est_method3 <- "Propensity score estimated using: Maximum Likelihood"
    } else if (x$argu$est_method[1] == 'ipw') {
      est_method1 <- "IPW-DDD estimation for the ATT(g,t): \n"
      est_method2 <- "Outcome Regression estimated using: N/A"
      est_method3 <- "Propensity score estimated using: Maximum Likelihood"
    } else if (x$argu$est_method[1] == 'reg') {
      est_method1 <- "REG-DDD estimation for the ATT(g,t): \n"
      est_method2 <- "Outcome Regression estimated using: OLS"
      est_method3 <- "Propensity score estimated using: N/A"
    } #else if (x$argu$est_method[1] == 'dml') {
    #   est_method1 <- "DML-DDD estimation for the ATT(g,t): \n"
    #   est_method2 <- paste("Outcome Regression estimated using:", x$argu$learners$ml_md$label)
    #   est_method3 <- paste("Propensity score estimated using:", x$argu$learners$ml_pa$label)
    # }

    # Front-end Summary Table
    lev_conf <- paste0(round(100 * (1 - x$argu$alpha), digits = 2), "% ")
    band_method <- ifelse(x$argu$cband, "Simult. ", "Pointwise ")
    interval_text <- paste0("[", lev_conf, band_method)
    header <- c("Group", "Time", "ATT(g,t)", "Std. Error", interval_text, "Conf. Band]", "")

    sig <- (x$uci < 0) | (x$lci > 0)
    sig[is.na(sig)] <- FALSE
    sig_text <- ifelse(sig, "*", "")

    body <- cbind.data.frame(
      x$groups,
      x$periods,
      sprintf("%10.4f", x$ATT),
      sprintf("%10.4f", x$se),
      sprintf("%10.4f", x$lci),
      sprintf("%10.4f", x$uci),
      sig_text
    )
    colnames(body) <- header
  }

  # Printing results in console
  cat(" Call:\n")
  print(x$call.params)
  cat("=========================== DDD Summary ==============================")
  cat("\n", est_method1)
  utils::write.table(format(rbind(header, body), justify = "centre", digits = 4, nsmall = 4),
                     row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")

  cat("\n")
  cat(" Note: * indicates that the confidence interval does not contain zero.")

  cat("\n --------------------------- Data Info   -----------------------------")
  # Panel data
  cat("\n", "Panel data")
  cat("\n", paste0("Outcome variable: ", x$argu$yname))
  # add partition variable name
  cat("\n", paste0("Qualification variable: ", x$argu$pname))
  # add control group for multiple periods
  if(x$argu$multiple_periods == TRUE){
    ifelse(x$argu$control_group == "nevertreated", control_type <- "Never Treated", control_type <- "Not yet Treated (GMM-based)")
    cat("\n", paste0("Control group: ", control_type))
  }
  if(x$argu$multiple_periods == FALSE){
    cat("\n", "No. of units at each subgroup:")
    cat("\n", paste0("  treated-and-eligible: ", x$subgroup_counts$V1[1]))
    cat("\n", paste0("  treated-but-ineligible: ", x$subgroup_counts$V1[2]))
    cat("\n", paste0("  eligible-but-untreated: ", x$subgroup_counts$V1[3]))
    cat("\n", paste0("  untreated-and-ineligible: ", x$subgroup_counts$V1[4]))
  } else {
    cat("\n", "No. of units per treatment group:")
    for (i in 1:nrow(x$cohort_size)) {

      if (x$cohort_size$first_treat[i] == 0) {
        cat("\n", paste0("  Units never enabling treatment: ", x$cohort_size$V1[i]), sep = "")
      } else {
        cat("\n", paste0("  Units enabling treatment at period ", x$cohort_size$first_treat[i], ": ", x$cohort_size$V1[i]), sep = "")
      }
    }
  }

  # TODO: ADD number of covariates and some examples

  cat("\n --------------------------- Algorithms ------------------------------")
  # Estimation Method
  cat("\n", est_method2)
  cat("\n", est_method3)

  # if (x$argu$est_method[1] == 'dml'){
  #   cat("\n -------------------------- Cross-fitting  ---------------------------")
  #   cat("\n No. of folds: ", x$argu$n_folds)
  #   cat("\n Apply cross-fitting: TRUE")
  # }
  # Analytical vs bootstrapped standard errors
  cat("\n --------------------------- Std. Errors  ----------------------------")
  cat("\n", paste("Level of significance: ", x$argu$alpha))
  if (x$argu$boot == T) {
    boot1 <-
      cat(
        "\n Boostrapped standard error based on",
        x$argu$nboot, "reps.",
        "\n Method: Multiplier Bootstrap."
      )
  } else {
    boot1 <- cat("\n Analytical standard errors.")
  }
  cat("\n", "Type of confidence band: ", ifelse(x$argu$cband, "Uniform Confidence Band ", "Pointwise Confidence Interval"))
  if (!is.null(x$argu$cluster)){
    cat("\n", paste0("Clustering Std. Errors by: ", x$argu$cluster))
  }
  cat("\n =====================================================================")
  cat("\n See Ortiz-Villavicencio and Sant'Anna (2025) for details.")

}
