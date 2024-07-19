#' @title Print
#'
#' @description Prints a aggddd Object
#'
#' @param x A aggddd object
#' @param ... Other params (required as generic function, but not used)
#' @importFrom utils write.table
#' @export
#' @noRd

print.aggddd <- function(x, ...) {

  pointwise_cval <- qnorm(1-0.05/2)
  overall_cband_upper <- x$aggte$overall.att + pointwise_cval*x$aggte$overall.se
  overall_cband_lower <- x$aggte$overall.att - pointwise_cval*x$aggte$overall.se
  out1 <- cbind.data.frame(x$aggte$overall.att, x$aggte$overall.se, overall_cband_lower, overall_cband_upper)
  out1 <- round(out1, 4)
  overall_sig <- (overall_cband_upper < 0) | (overall_cband_lower > 0)
  overall_sig[is.na(overall_sig)] <- FALSE
  overall_sig_text <- ifelse(overall_sig, "*", "")
  out1 <- cbind.data.frame(out1, overall_sig_text)

  # Printing results in console
  cat(" Call:\n")
  print(x$call.params)
  cat("========================= DDD Aggregation ==========================")
  if (x$aggte$type == "simple") cat("\n", "Overall ATT:")
  if (x$aggte$type == "eventstudy") cat("\n", "Overall summary of ATT\'s based on event-study/dynamic aggregation: ")
  if (x$aggte$type == "group") cat("\n", "Overall summary of ATT\'s based on group/cohort aggregation: ")
  if (x$aggte$type == "calendar") cat("\n", "Overall summary of ATT\'s based on calendar time aggregation: ")

  # Front-end Summary Table
  cat("\n")
  header <-
    c(" ATT",
      "Std. Error",
      "[95% Conf.",
      "Interval]",
      "")
  colnames(out1) <- header

  # utils::write.table(format(rbind(header, out1), justify= "centre", digits=2, nsmall=2),
  #                    row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")
  print(out1, row.names=FALSE)

  # handle cases depending on type
  if (x$aggte$type %in% c("group","eventstudy","calendar")) {

    # header
    if (x$aggte$type=="eventstudy") { c1name <- "Event time"; cat("\n", "Event Study:") }
    if (x$aggte$type=="group") { c1name <- "Group"; cat("\n", "Group Effects:") }
    if (x$aggte$type=="calendar") { c1name <- "Time"; cat("\n", "Calendar Effects:") }

    cat("\n")
    cband_text1a <- paste0(100*(1-0.05),"% ")
    # cband_text1b <- ifelse(x$argu$boot,
    #                        ifelse(x$argu$cband, "Simult. ", "Pointwise "),
    #                        "Pointwise ")
    cband_text1b <- "Pointwise "
    cband_text1 <- paste0("[", cband_text1a, cband_text1b)

    cband_lower <- x$aggte$att.egt - x$aggte$crit.val.egt * x$aggte$se.egt
    cband_upper <- x$aggte$att.egt + x$aggte$crit.val.egt * x$aggte$se.egt

    sig <- (cband_upper < 0) | (cband_lower > 0)
    sig[is.na(sig)] <- FALSE
    sig_text <- ifelse(sig, "*", "")

    out2 <- cbind.data.frame(x$aggte$egt, x$aggte$att.egt, x$aggte$se.egt, cband_lower, cband_upper)
    out2 <- round(out2, 4)
    out2 <- cbind.data.frame(out2, sig_text)

    header2 <- c(c1name, "Estimate","Std. Error", cband_text1, "Conf. Band]", "")

    colnames(out2) <- header2
    # utils::write.table(format(rbind(header2, out2), justify= "centre", digits=2, nsmall=2),
    #                    row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")
    print(out2, row.names=FALSE, justify = "centre")
  }
  cat("\n")
  cat(" Note: * indicates that confidence interval does not contain zero.")

  cat("\n --------------------------- Data Info   --------------------------")
  cat("\n", paste("Outcome variable:", x$aggte$yname))
  # add partition variable name
  cat("\n", paste("Partition variable:", x$aggte$partition_name))
  ifelse(x$aggte$control_group == "nevertreated", control_type <- "Never Treated", control_type <- "Not yet Treated")
  cat("\n", paste("Control group: ", control_type))

  # Analytical vs bootstrapped standard errors
  cat("\n --------------------------- Std. Errors  -------------------------")
  if (x$argu$boot == T) {
    boot1 <-
      cat(
        "\n Boostrapped standard error based on",
        x$argu$nboot,
        "bootstrap draws. \n Bootstrap Method:",
        "Multiplier bootstrap",
        "(Mammen,1993) \n"
      )
  } else {
    boot1 <- cat("\n Analytical standard errors.")
  }
  cat("\n ==================================================================")
  cat("\n See Ortiz-Villavicencio and Sant'Anna (2024) for details.")
}
