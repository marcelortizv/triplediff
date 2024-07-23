#' @title Print
#'
#' @description Prints a agg_ddd Object
#'
#' @param x A agg_ddd object
#' @param ... Other params (required as generic function, but not used)
#' @importFrom utils write.table
#' @export
#' @noRd

print.agg_ddd <- function(x, ...) {

  pointwise_cval <- qnorm(1-x$aggte_ddd$argu$alpha/2)
  overall_cband_upper <- x$aggte_ddd$overall.att + pointwise_cval * x$aggte_ddd$overall.se
  overall_cband_lower <- x$aggte_ddd$overall.att - pointwise_cval * x$aggte_ddd$overall.se
  out1 <- cbind.data.frame(x$aggte_ddd$overall.att, x$aggte_ddd$overall.se, overall_cband_lower, overall_cband_upper)
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
  header <- c("ATT","   Std. Error", paste0("    [ ",100*(1-x$aggte_ddd$argu$alpha),"% "), "Conf. Int.]","")
  colnames(out1) <- header

  # utils::write.table(format(rbind(header, out1), justify= "centre", digits=2, nsmall=2),
  #                    row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")
  print(out1, row.names=FALSE)

  # handle cases depending on type
  if (x$aggte_ddd$type %in% c("group","eventstudy","calendar")) {

    # header
    if (x$aggte_ddd$type=="eventstudy") { c1name <- "Event time"; cat("\n", "Event Study:") }
    if (x$aggte_ddd$type=="group") { c1name <- "Group"; cat("\n", "Group Effects:") }
    if (x$aggte_ddd$type=="calendar") { c1name <- "Time"; cat("\n", "Calendar Effects:") }

    cat("\n")

    # Front-end Summary Table
    lev_conf <- paste0(round(100 * (1 - x$aggte_ddd$argu$alpha), digits = 2), "% ")
    band_method <- ifelse(x$aggte_ddd$argu$cband, "Simult. ", "Ptwise. ")
    interval_text <- paste0("[", lev_conf, band_method)


    # cband_text1a <- paste0(100*(1-x$aggte_ddd$argu$alpha),"% ")
    # # cband_text1b <- ifelse(x$argu$boot,
    # #                        ifelse(x$argu$cband, "Simult. ", "Pointwise "),
    # #                        "Pointwise ")
    # cband_text1b <- "Pointwise "
    # cband_text1 <- paste0("[", cband_text1a, cband_text1b)

    cband_lower <- x$aggte_ddd$att.egt - x$aggte_ddd$crit.val.egt * x$aggte_ddd$se.egt
    cband_upper <- x$aggte_ddd$att.egt + x$aggte_ddd$crit.val.egt * x$aggte_ddd$se.egt

    sig <- (cband_upper < 0) | (cband_lower > 0)
    sig[is.na(sig)] <- FALSE
    sig_text <- ifelse(sig, "*", "")

    out2 <- cbind.data.frame(x$aggte_ddd$egt, x$aggte_ddd$att.egt, x$aggte_ddd$se.egt, cband_lower, cband_upper)
    out2 <- round(out2, 4)
    out2 <- cbind.data.frame(out2, sig_text)

    header2 <- c(c1name, "Estimate","Std. Error", interval_text, "Conf. Band]", "")

    colnames(out2) <- header2
    # utils::write.table(format(rbind(header2, out2), justify= "centre", digits=2, nsmall=2),
    #                    row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")
    print(out2, row.names=FALSE, justify = "centre")
  }
  cat("\n")
  cat(" Note: * indicates that confidence interval does not contain zero.")

  cat("\n --------------------------- Data Info   --------------------------")
  cat("\n", paste("Outcome variable:", x$aggte_ddd$yname))
  # add partition variable name
  cat("\n", paste("Partition variable:", x$aggte_ddd$partition_name))
  ifelse(x$aggte_ddd$control_group == "nevertreated", control_type <- "Never Treated", control_type <- "Not yet Treated")
  cat("\n", paste("Control group: ", control_type))
  cat("\n", paste("Level of significance: ", x$aggte_ddd$argu$alpha))

  # Analytical vs bootstrapped standard errors
  cat("\n --------------------------- Std. Errors  -------------------------")
  if (x$aggte_ddd$argu$boot == T) {
    boot1 <-
      cat(
        "\n Boostrapped standard error based on",
        x$aggte_ddd$argu$nboot, "reps.",
        "\n Bootstrap Method: Multiplier bootstrap."
      )
  } else {
    boot1 <- cat("\n Analytical standard errors.")
  }
  if (!is.null(x$aggte_ddd$argu$cluster)){
    cat("\n", paste0("Clustering Std. Errors by: ", x$aggte_ddd$argu$cluster))
  }
  cat("\n ==================================================================")
  cat("\n See Ortiz-Villavicencio and Sant'Anna (2024) for details.")
}
