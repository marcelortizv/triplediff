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
  if (x$argu$multiple.periods == FALSE){
    if (x$argu$estMethod[1] == 'trad') {
      estMethod1 <- "DRDDD estimator for the ATT: \n"
      estMethod2 <- "Outcome Regression estimated using: OLS"
      estMethod3 <- "Propensity score estimated using: Maximum Likelihood"
    } else {
      estMethod1 <- "DMLDDD estimator for the ATT: \n"
      estMethod2 <- paste("Outcome Regression estimated using:", x$argu$learners$ml_pa$learner$label)
      estMethod3 <- paste("Propensity score estimated using:", x$argu$learners$ml_ma$learner$label)
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

  # TODO: Fix this when multiple.periods is done
  # if (x$argu$multiple.periods == TRUE){
  #   if (x$argu$estMethod[1] == 'trad') {
  #     estMethod1 <- "DRDDD estimator for the ATT(gt): \n"
  #     estMethod2 <- "Outcome Regression estimated using: OLS"
  #     estMethod3 <- "Propensity score estimated using: Maximum Likelihood"
  #   } else {
  #     estMethod1 <- "DMLDDD estimator for the ATT(gt): \n"
  #     #TODO: Add the method used to estimate OR and PS (name should be contained in "learners")
  #     estMethod2 <- paste("Outcome Regression estimated using:", x$argu$learners[1]$name)
  #     estMethod3 <- paste("Propensity score estimated using:", x$argu$learners[1]$name)
  #   }
  #
  #   # Front-end Summary Table
  #   header <-
  #     c("Group",
  #       "Time",
  #       "ATT(g,t)",
  #       "Std. Error",
  #       # "t value",
  #       # "Pr(>|t|)",
  #       "[95% Conf.",
  #       "Interval]")
  #   #TODO: Adjust the bod to take into account the group and time variables
  #   body <- cbind(
  #     round(x$ATT, digits = 4),
  #     round(x$se, digits = 4),
  #     round(x$ATT / x$se, digits = 4),
  #     round(2 * (stats::pnorm(-abs(
  #       x$ATT / x$se
  #     ))), digits = 4),
  #     round(x$lci, digits = 4),
  #     round(x$uci, digits = 4)
  #   )
  #   colnames(body) <- header
  # }

  # Printing results in console
  cat(" Call:\n")
  print(x$call.params)
  cat("=========================== DDD Summary ==========================")
  cat("\n", estMethod1, "\n")
  cat("--------------------------- Fit Results --------------------------")
  utils::write.table(format(rbind(header, body), justify= "centre", digits=2, nsmall=2),
                     row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")
  cat("--------------------------- Algorithm   --------------------------")
  # Estimation Method
  cat("\n", estMethod2)
  cat("\n", estMethod3)
  cat("--------------------------- Data Info   --------------------------")
  # Panel data
  cat("\n", "Panel data")
  cat("\n", paste("Outcome variable: ", x$argu$yname))
  # add partition variable name
  cat("\n", paste("Partition variable: ", x$argu$partition.name))

  # add control group for multiple periods
  if(x$argu$multiple.periods == TRUE){
    cat("\n", paste("Control group: ", x$argu$control.group))
  }

  # TODO: ADD number of observations in each partition
  # TODO: ADD number of covariates and some examples

  if (x$argu$estMethod[1] == 'dml'){
    cat("--------------------------- Cross-fitting  -------------------------")
    cat("\n No. of folds: ", x$argu$n_folds)
    cat("\n Apply cross-fitting: TRUE")
  }
  # Analytical vs bootstrapped standard errors
  cat("--------------------------- Std. Errors  -------------------------")
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
    boot1 <- cat("\n Analytical standard errors.\n")
  }
  cat("------------------------------------------------------------------")
  cat("\n See Ortiz-Villavicencio and Sant'Anna (2024) for details.")
}
