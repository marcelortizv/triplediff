#' Doubly robust DDD estimator for ATT, with panel data and multiple periods
#'
#' This function implements a doubly robust estimator for assessing the average
#' treatment effect on the treated (ATT) using a triple differences (DDD) approach
#' in panel data settings with staggered treatment adoption. The function takes preprocessed
#' data structured specifically for this analysis.
#'
#' @importFrom Matrix Matrix
#' @param did_preprocessed A list containing preprocessed data and specifications for the DDD estimation.
#'        Expected elements include:
#'        - `preprocessed_data`: A data table containing the data with variables needed for the analysis.
#'        - `xformula`: The formula for the covariates to be included in the model. It should be of the form \code{~ x1 + x2}.
#'        Default is \code{xformla = ~1} (no covariates).
#'        - `control_group`: A character string indicating the control group. Default is \code{"notyettreated"}.
#'        - `n`: The number of unique ID observations in the data.
#'        - `cohorts`: A vector containing the number of treated cohorts
#'        - `time_periods`: A vector containing the number of time periods
#'        - `tlist`: A vector containing the time periods
#'        - `glist`: A vector containing the groups
#'        - `base_period`: A character string indicating the base period. Default is \code{"universal"}.
#'        - `cohort_size`: A vector containing the size of each cohort
#'        - `boot`: Logical. If \code{TRUE}, the function computes the bootstrap standard errors. Default is \code{FALSE}.
#'        - `nboot`: The number of bootstrap samples to be used. Default is \code{NULL}. If \code{boot = TRUE}, the default is \code{nboot = 999}.
#'        - `use_parallel`: Boolean of whether or not to use parallel processing in the multiplier bootstrap, default is \code{use_parallel=FALSE}
#'        - `cores`: the number of cores to use with parallel processing, default is \code{cores=1}
#'        - `cband`: Boolean of whether or not to compute simultaneous confidence bands, default is \code{cband=FALSE}
#'        - `alpha` The level of significance for the confidence intervals.  Default is \code{0.05}.
#'
#' @keywords internal
#' @return A list with the estimated ATT, standard error, upper and lower confidence intervals, and influence function.
#' @noRd
NULL
# ------------------------------------------------------------------------------

att_gt <- function(did_preprocessed){
  data <- did_preprocessed$preprocessed_data
  xformla <- did_preprocessed$xformula
  control_group <- did_preprocessed$control_group
  n <- did_preprocessed$n
  cohorts <- did_preprocessed$cohorts
  time_periods <- did_preprocessed$time_periods
  tlist <- did_preprocessed$tlist
  glist <- did_preprocessed$glist
  base_period <- did_preprocessed$base_period
  boot <- did_preprocessed$boot
  alpha <- did_preprocessed$alpha
  cohort_size <- did_preprocessed$cohort_size
  use_parallel <- did_preprocessed$use_parallel # to perform bootstrap
  cores <- did_preprocessed$cores # to perform bootstrap
  cband <- did_preprocessed$cband # to perform bootstrap + simult. conf. band

  orig_data <- copy(data)

  attgt_list <- list()

  counter <- 1
  tlist_length <- length(tlist)
  tfac <- 0
  if (base_period != "universal") {
    tlist_length <- tlist_length - 1
    tfac <- 1
  }

  # influence function matrix to be populated
  inf_func_mat <- Matrix::Matrix(data=0,nrow=n, ncol=cohorts*(time_periods - tfac), sparse=TRUE)
  se_gt_ddd_nyt <- rep(NA, cohorts * (time_periods - tfac)) # standard errors for not yet treated control group
  # If "nevertreated", create the control column that gonna indicate the comparison group for each iteration
  if (control_group == "nevertreated") {
    data[, control := as.integer(first_treat == 0)]
  }


  for (g in 1:cohorts){
    data[, treat := as.integer(first_treat == glist[g])]

    for (t in 1:tlist_length){
      pret <- t
      if (base_period == "universal"){
        pret <- tail(which(tlist < glist[g]), 1)
      }

      # creating control in case of "not yet treated option"
      if (control_group == "notyettreated") {
        # max_period <- tlist[max(t, pret) + tfac] + anticipation
        max_period <- tlist[max(t, pret) + tfac]
        data[, control := as.integer((first_treat == 0) | (first_treat > max_period & first_treat != glist[g]))]
      }


      # check if in post-treatment period
      if ((glist[g] <= tlist[(t + tfac)])) {

        # update pre-period if in post-treatment period to be period (g-delta-1)
        pret <- tail(which(tlist < glist[g]), 1)

        # print a warning message if there are no pre-treatment period
        if (length(pret) == 0) {
          warning(paste0("There are no pre-treatment periods for the group first treated at ", glist[g], "\nUnits from this group are dropped"))

          # if there are not pre-treatment periods, move on to next iteration
          break
        }
      }

      # normalize in period (g-1) to be equal to 0 when based period is universal and move on to next iteration
      if (base_period == "universal") {
        if (tlist[pret] == tlist[t + tfac]) {
          attgt_list[[counter]] <- list(att = 0, group = glist[g], year = tlist[t + tfac], post = 0)
          inf_func_mat[, counter] <- rep(0, n)
          counter <- counter + 1
          next
        }
      }

      # post treatment dummy variable
      post_treat <- 1*(glist[g] <= tlist[t + tfac])

      # filter data with only 2 periods
      cohort_data <- data[period %in% c(tlist[t + tfac], tlist[pret])]

      # -------------------------------------
      # PANEL DATA ONLY
      # -------------------------------------

      # filter data for treated and control groups in each (g,t) cell. Save index

      # get total number of units
      n_size = uniqueN(cohort_data[, id])

      # index of unit in current cell (g,t) when treat = 1 and control = 1
      index_units_in_gt <- cohort_data[, treat == 1 | control == 1]

      # filter by units selected in current cell (g,t)
      cohort_data <- cohort_data[index_units_in_gt]

      # save number of units after filtering
      size_gt = uniqueN(cohort_data[, id])

      # Identify all available control states in this cohort.
      available_controls <- sort(unique(cohort_data$first_treat[cohort_data$control == 1 ]))
      available_controls <- available_controls[available_controls != glist[g]]

      if (length(available_controls) == 1) {
        # creating subgroup variable
        # 4 if (eligible ==1 & enabled == 1); 3 if (eligible ==0 & enabled == 1); 2 if (eligible ==1 & enabled == 0); 1 if (eligible ==0 & enabled == 0)
        # Create subgroup variable based on the control_group option
        cohort_data[, subgroup := NA_integer_]
        cohort_data[(treat == 1 & partition == 1), subgroup := 4]
        cohort_data[(treat == 1 & partition == 0), subgroup := 3]
        cohort_data[(treat == 0 & partition == 1), subgroup := 2]
        cohort_data[(treat == 0 & partition == 0), subgroup := 1]

        # add post treatment dummy variable if period is equal to t + tfac
        cohort_data[, post := 0]

        # Ypre is always the baseline, even in pre-treatment periods were baseline could be in later periods with respect to t.
        pseudo_post = ifelse(post_treat == 1, tlist[max(t, pret) + tfac], tlist[t + tfac])
        cohort_data[(period == pseudo_post), post := 1]

        # cohort_data[(period == tlist[max(t, pret) + tfac]), post := 1]
        # Calculate the size of each subgroup in the 'subgroup' column
        subgroup_counts <- cohort_data[, .N/2, by = subgroup][order(-subgroup)]

        # Reassign dp object and run att_dr
        did_preprocessed$preprocessed_data <- cohort_data
        did_preprocessed$subgroup_counts <- subgroup_counts
        did_preprocessed$boot <- FALSE # forcing false to avoid bootstrap inside att_dr() function.
        did_preprocessed$inffunc <- TRUE # forcing true to recover influence function

        # run att_dr based on 2 time periods subdata
        attgt_inf_func <- att_dr(did_preprocessed)

        # adjust influence function
        attgt_inf_func$inf_func <- (n_size/size_gt) * attgt_inf_func$inf_func
        # save results in a list
        attgt_list[[counter]] <- list(att = attgt_inf_func$ATT, group = glist[g], year = tlist[t + tfac], post = post_treat)

        # recover influence function
        inff <- rep(0, n_size)
        # avoid repetition in index
        inff[index_units_in_gt[seq(1, length(index_units_in_gt), by = 2)]] <- attgt_inf_func$inf_func
        # save in influence function matrix
        inf_func_mat[, counter] <- inff
        # update counter
        counter  <- counter + 1
      } else {
        # Here we have more than one control group available -> notyettreated

        # Initialize local containers.
        ddd_over_controls_res <- list()
        inf_mat_local <- NULL

        # Loop over each available control state.
        for (ctrl in available_controls) {
          # Subset data: keep treated units (first_treat == glist[g]) and control units (first_treat == ctrl)

          index_units_in_gt_ctrl <- cohort_data[, first_treat == glist[g] | first_treat == ctrl]
          subset_data <- cohort_data[index_units_in_gt_ctrl]
          size_gt_ctrl <- uniqueN(subset_data[, id])

          # Create subgroup variable based on the control_group option
          subset_data[, subgroup := NA_integer_]
          subset_data[(treat == 1 & partition == 1), subgroup := 4]
          subset_data[(treat == 1 & partition == 0), subgroup := 3]
          subset_data[(treat == 0 & partition == 1), subgroup := 2]
          subset_data[(treat == 0 & partition == 0), subgroup := 1]

          # add post treatment dummy variable if period is equal to t + tfac
          subset_data[, post := 0]

          # Ypre is always the baseline, even in pre-treatment periods were baseline could be in later periods with respect to t.
          pseudo_post = ifelse(post_treat == 1, tlist[max(t, pret) + tfac], tlist[t + tfac])
          subset_data[(period == pseudo_post), post := 1]

          # subset_data[(period == tlist[max(t, pret) + tfac]), post := 1]
          # Calculate the size of each subgroup in the 'subgroup' column
          subgroup_counts <- subset_data[, .N/2, by = subgroup][order(-subgroup)]

          # Reassign dp object and run att_dr
          did_preprocessed$preprocessed_data <- subset_data
          did_preprocessed$subgroup_counts <- subgroup_counts
          did_preprocessed$boot <- FALSE # forcing false to avoid bootstrap inside att_dr() function.
          did_preprocessed$inffunc <- TRUE # forcing true to recover influence function


          # run att_dr based on 2 time periods subdata
          ddd_out <- triplediff::att_dr(did_preprocessed)

          # computing the aggregate influence function for every g'>t: RIF_{g,g',t}
          # adjust influence function (rescale + recentered)
          ddd_out$inf_func <- ((size_gt/size_gt_ctrl)*ddd_out$inf_func) #+ 1*ddd_out$ATT

          # Save the ddd output in the local list using the control state as name.
          ddd_over_controls_res[[as.character(ctrl)]] <- ddd_out

          # recover influence function
          inff <- rep(0, n_size)
          # avoid repetition in index
          inff[index_units_in_gt_ctrl[seq(1, length(index_units_in_gt_ctrl), by = 2)]] <- ddd_out$inf_func

          # Append this vector as a new column in the local matrix.
          if (is.null(inf_mat_local)) {
            inf_mat_local <- matrix(inff, ncol = 1)
            colnames(inf_mat_local) <- as.character(ctrl)
          } else {
            inf_mat_local <- cbind(inf_mat_local, inff)
            colnames(inf_mat_local)[ncol(inf_mat_local)] <- as.character(ctrl)
          }
        } # end loop over available_controls

        # If no ddd estimations were obtained, skip this (g,t) cell.
        if (length(ddd_over_controls_res) == 0) next

        # Compute the naive ATT as the average of all ddd ATT estimates.
        att_vals <- sapply(ddd_over_controls_res, function(res) res$ATT)

        # -----------------------------------------
        # Compute the GMM-based estimator for ATT(g,t)
        # -----------------------------------------

        # computing the aggregate influence function for every g'>t: RIF_{g,g',t}
        OMEGA <- cov(inf_mat_local, use = "complete.obs")

        ones <- rep(1, length(ddd_over_controls_res))
        inv_OMEGA <- solve(OMEGA)
        w <- colSums(inv_OMEGA) / sum(inv_OMEGA)

        # Compute the GMM estimator for ATT(g,t)
        ATT_gmm <- sum(w * att_vals) / sum(w)
        # IF_gmm
        IF_gmm <- (inf_mat_local) %*% w
        # gmm_se
        gmm_se <- sqrt(1/ (n * sum(inv_OMEGA)))

        # save the results for this (g,t) cell
        attgt_list[[counter]] <- list(att = ATT_gmm, group = glist[g], year = tlist[t + tfac], post = post_treat)

        # save in influence function matrix
        inf_func_mat[, counter] <- IF_gmm

        # save the standard errors from GMM
        se_gt_ddd_nyt[counter] <- gmm_se

        # update counter
        counter  <- counter + 1

      } # end of GMM-based procedure

    } # end of tlist loop
  } # end of glist loop

  # recover original arguments in did_preprocessed
  did_preprocessed$preprocessed_data <- orig_data
  did_preprocessed$boot <- boot # restoring boot argument

  # PREPROCESS attgt_list AND inf_func_mat
  attgt_res <- process_attgt(attgt_list)
  groups <- attgt_res$group
  periods <- attgt_res$periods
  att_gt_ddd <- attgt_res$att
  # NEW PROCEDURE

  # Get analytical errors: This is analogous to cluster robust standard errors at the unit level
  n <- did_preprocessed$n
  V <- Matrix::t(inf_func_mat)%*%inf_func_mat/n
  se_gt_ddd <- sqrt(Matrix::diag(V)/n)

  # Zero standard error replaced by NA
  se_gt_ddd[se_gt_ddd <= sqrt(.Machine$double.eps)*10] <- NA

  # if control group is "notyettreated", we need to adjust the standard errors
  if (control_group == "notyettreated") {
    # adjust standard errors for not yet treated control group
    se_gt_ddd[!is.na(se_gt_ddd_nyt)] <- se_gt_ddd_nyt[!is.na(se_gt_ddd_nyt)]
  }

  # Identify entries of main diagonal V that are zero or NA
  zero_na_sd_entry <- unique(which(is.na(se_gt_ddd)))

  # compute bootstrap standard errors
  bT <- NULL
  if (boot){
    # perform multiplier bootstrap
    boot_result <- mboot(inf_func_mat, did_preprocessed=did_preprocessed, use_parallel=use_parallel, cores=cores)
    bT <- boot_result$bT # sup-t confidence band
    # save bootstrap standard errors
    if(length(zero_na_sd_entry)>0) {
      se_gt_ddd[-zero_na_sd_entry] <- boot_result$se[-zero_na_sd_entry]
    } else {
      se_gt_ddd <- boot_result$se
    }
  }

  # Zero standard error replaced by NA
  se_gt_ddd[se_gt_ddd <= sqrt(.Machine$double.eps)*10] <- NA

  #-----------------------------------------------------------------------------
  # compute confidence intervals / bands
  #-----------------------------------------------------------------------------

  # compute critical values for point-wise interval
  cv <- qnorm(1 - alpha/2)

  # in case uniform confidence bands are requested
  if (boot){
    if (cband){
      # get critical value to compute uniform confidence bands
      cv <- boot_result$unif_crit_val
      if(cv >= 7){
        warning("Simultaneous critical value is arguably `too large' to be realible. This usually happens when number of observations per group is small and/or there is no much variation in outcomes.")
      }
    }
  }

  # compute confidence intervals/bands
  ci_upper <- att_gt_ddd + cv * se_gt_ddd
  ci_lower <- att_gt_ddd - cv * se_gt_ddd

  # ------------------------------------------------------------------------------
  # Return results
  # ------------------------------------------------------------------------------

  # we need this for aggregation
  first_period_dta <- orig_data[period == tlist[1]]
  ret <- (list(ATT = att_gt_ddd,
               se = se_gt_ddd,
               uci = ci_upper,
               lci = ci_lower,
               groups = groups,
               periods = periods,
               tlist = tlist,
               glist = glist,
               cohort_size = cohort_size,
               n = n,
               bT = bT,
               inf_func_mat = inf_func_mat,
               first_period_dta = first_period_dta
              ))

  return(ret)
}
