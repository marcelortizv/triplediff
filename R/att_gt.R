att_gt_dr <- function(did_preprocessed){
  data <- did_preprocessed$preprocessed_data
  xformla <- did_preprocessed$xformula
  control_group <- did_preprocessed$control_group
  n <- did_preprocessed$n
  cohorts <- did_preprocessed$cohorts
  time_periods <- did_preprocessed$time_periods
  tlist <- did_preprocessed$tlist
  glist <- did_preprocessed$glist
  base_period <- did_preprocessed$base_period

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
          inf_func_mat[, counter] <- rep(0,n)
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
      size = uniqueN(cohort_data[, id])

      # index of unit in current cell (g,t) when treat = 1 and control = 1
      index_units_in_gt <- cohort_data[, treat == 1 | control == 1]

      # filter by units selected in current cell (g,t)
      cohort_data <- cohort_data[index_units_in_gt]

      # save number of units after filtering
      size_gt = uniqueN(cohort_data[, id])

      # creating subgroup variable
      # 4 if (partition ==1 & treat == 1); 3 if (partition ==0 & treat == 1); 2 if (partition ==1 & treat == 0); 1 if (partition ==0 & treat == 0)
      # Create subgroup variable based on the control_group option
      cohort_data[, subgroup := NA_integer_]
      cohort_data[(treat == 1 & partition == 1), subgroup := 4]
      cohort_data[(treat == 1 & partition == 0), subgroup := 3]
      cohort_data[(treat == 0 & partition == 1), subgroup := 2]
      cohort_data[(treat == 0 & partition == 0), subgroup := 1]

      # add post treatment dummy variable if period is equal to t + tfac
      cohort_data[, post := 0]
      cohort_data[(period == tlist[max(t, pret) + tfac]), post := 1]
      # Calculate the size of each subgroup in the 'subgroup' column
      subgroup_counts <- cohort_data[, .N/2, by = subgroup][order(-subgroup)]

      # Reassign dp object and run att_dr
      dp$preprocessed_data <- cohort_data
      dp$subgroup_counts <- subgroup_counts
      dp$boot <- FALSE # forcing false to avoid bootstrapping inside att_dr() function.

      att_gt <- att_dr(dp)

      # adjust influence function
      att_gt$inf_func <- (size/size_gt) * att_gt$inf_func
      # save results in a list
      attgt_list[[counter]] <- list(att = att_gt$ATT, group = glist[g], year = tlist[t + tfac], post = post_treat)

      # recover influence function
      inff <- rep(0, size)
      # avoid repetition in index
      inff[index_units_in_gt[seq(1, length(index_units_in_gt), by = 2)]] <- att_gt$inf_func
      # save in influence function matrix
      inf_func_mat[, counter] <- inff
      # update counter
      counter  <- counter + 1

    } # end of tlist loop
  } # end of glist loop

  # TODO; PREPROCESS attgt_list AND inf_func_mat

  # TODO; COMPUTE STD ERRORS; EITHER ANALYTICA

  # TODO; COMPUTE CONFIDENCE BANDS

  # TODO; RETURN RESULTS

  ret <- list()
  return(ret)
}
