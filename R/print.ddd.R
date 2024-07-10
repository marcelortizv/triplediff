#' @title Print
#'
#' @description Prints a ddd Object
#'
#' @param x A ddd object
#' @param ... Other params (required as generic function, but not used)
#' @importFrom utils write.table
#' @export
#' @noRd



print.ddd <- function(x, ...) {
  # Methods used
  if (x$argu$multiple_periods == FALSE){
    if (x$argu$est_method[1] == 'trad') {
      est_method1 <- "DRDDD estimator for the ATT: \n"
      est_method2 <- "Outcome Regression estimated using: OLS"
      est_method3 <- "Propensity score estimated using: Maximum Likelihood"
    } else {
      est_method1 <- "DMLDDD estimator for the ATT: \n"
      est_method2 <- paste0("Outcome Regression estimated using: ", x$argu$learners$ml_pa$label)
      est_method3 <- paste0("Propensity score estimated using: ", x$argu$learners$ml_md$label)
    }

    # Front-end Summary Table
    header <-
      c("ATT",
        "Std. Error",
        "t value",
        "Pr(>|t|)",
        "[95% Conf.",
        "Interval]")
    body <- cbind(
      round(x$ATT, digits = 4),
      round(x$se, digits = 4),
      round(x$ATT / x$se, digits = 4),
      round(2 * (stats::pnorm(-abs(
        x$ATT / x$se
      ))), digits = 4),
      round(x$lci, digits = 4),
      round(x$uci, digits = 4)
    )
    colnames(body) <- header
  }

  if (x$argu$multiple_periods == TRUE){
    if (x$argu$est_method[1] == 'trad') {
      est_method1 <- "DRDDD estimator for the ATT(gt): \n"
      est_method2 <- "Outcome Regression estimated using: OLS"
      est_method3 <- "Propensity score estimated using: Maximum Likelihood"
    } else {
      est_method1 <- "DMLDDD estimator for the ATT(gt): \n"
      #TODO: Add the method used to estimate OR and PS (name should be contained in "learners")
      est_method2 <- paste("Outcome Regression estimated using:", x$argu$learners[1]$name)
      est_method3 <- paste("Propensity score estimated using:", x$argu$learners[1]$name)
    }

    # Front-end Summary Table
    header <-
      c("Group",
        "Time",
        "ATT(g,t)",
        "Std. Error",
        "[95% Conf.",
        "Interval]")
    body <- cbind.data.frame(
      x$groups,
      x$periods,
      round(x$ATT, digits = 4),
      round(x$se, digits = 4),
      round(x$lci, digits = 4),
      round(x$uci, digits = 4)
    )
    colnames(body) <- header
  }

  # Printing results in console
  cat(" Call:\n")
  print(x$call.params)
  cat("=========================== DDD Summary ===========================")
  cat("\n", est_method1)
  utils::write.table(format(rbind(header, body), justify= "centre", digits=2, nsmall=2),
                     row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")
  cat("\n --------------------------- Data Info   --------------------------")
  # Panel data
  cat("\n", "Panel data")
  cat("\n", paste("Outcome variable: ", x$argu$yname))
  # add partition variable name
  cat("\n", paste("Partition variable: ", x$argu$partition_name))
  if(x$argu$multiple_periods == FALSE){
    cat("\n", "No. of observations for each partition:")
    cat("\n", paste("(treat = 1, partition = 1): ", x$subgroup_counts$V1[1]))
    cat("\n", paste("(treat = 1, partition = 0): ", x$subgroup_counts$V1[2]))
    cat("\n", paste("(treat = 0, partition = 1): ", x$subgroup_counts$V1[3]))
    cat("\n", paste("(treat = 0, partition = 0): ", x$subgroup_counts$V1[4]))
  } else {
    cat("\n", "No. of observations per cohort:")
    for (i in 1:nrow(x$cohort_size)) {
       cat(paste0("Cohort ", x$cohort_size$first_treat[i], ": ", x$cohort_size$V1[i], "\n"), sep = "")
      }
  }
  # add control group for multiple periods
  if(x$argu$multiple_periods == TRUE){
    cat("\n", paste("Control group: ", x$argu$control_group))
  }
  # TODO: ADD number of observations in each partition
  # TODO: ADD number of covariates and some examples

  cat("\n --------------------------- Algorithm   --------------------------")
  # Estimation Method
  cat("\n", est_method2)
  cat("\n", est_method3)

  if (x$argu$est_method[1] == 'dml'){
    cat("\n -------------------------- Cross-fitting  ------------------------")
    cat("\n No. of folds: ", x$argu$n_folds)
    cat("\n Apply cross-fitting: TRUE")
  }
  # Analytical vs bootstrapped standard errors
  cat("\n --------------------------- Std. Errors  -------------------------")
  if (x$argu$boot == T) {
    boot1 <-
      cat(
        "\n Boostrapped standard error based on",
        x$nboot,
        "bootstrap draws. \n Bootstrap Method:",
        x$argu$boot.type[1],
        "(Mammen,1993) \n"
      )
  } else {
    boot1 <- cat("\n Analytical standard errors.")
  }
  cat("\n ==================================================================")
  cat("\n See Ortiz-Villavicencio and Sant'Anna (2024) for details.")
}
