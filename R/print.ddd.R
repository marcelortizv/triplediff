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
      #TODO: Add the method used to estimate OR and PS (name should be contained in "learners")
      estMethod2 <- paste("Outcome Regression estimated using:", x$argu$learners[1]$name)
      estMethod3 <- paste("Propensity score estimated using:", x$argu$learners[1]$name)
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
  cat("------------------------------------------------------------------")
  cat("\n", estMethod1, "\n")
  utils::write.table(format(rbind(header, body), justify= "centre", digits=2, nsmall=2),
                     row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")
  cat("------------------------------------------------------------------")
  # Panel data
  cat("\n", "Panel data")
  # Estimation Method
  cat("\n", estMethod2)
  cat("\n", estMethod3)

  # add partition variable name
  cat("\n", paste("Partition variable: ", x$argu$partition.name))

  # add control group for multiple periods
  if(x$argu$multiple.periods == TRUE){
    cat("\n", paste("Control group: ", x$argu$control.group))
  }

  # Analytical vs bootstrapped standard errors
  if (x$argu$boot == T) {
    boot1 <-
      cat(
        "\n Boostrapped standard error based on",
        x$argu$nboot,
        "bootstrap draws. \n Bootstrap method:",
        x$argu$boot.type[1],
        ". \n"
      )
  } else {
    boot1 <- cat("\n Analytical standard errors.\n")
  }
  cat("------------------------------------------------------------------")
  cat("\n See Ortiz-Villavicencio and Sant'Anna (2024) for details.")
}
