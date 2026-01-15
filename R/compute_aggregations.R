#' Utility functions to compute aggregation procedures.
#' @importFrom stats update
#' @import data.table
#' @noRd
#--------------------------------------------------

# --------------------------------------
# FUNCTIONS FOR AGGREGATION PROCEDURES
# --------------------------------------

#' @title Compute extra term in influence function due to estimating weights
#'
#' @description A function to compute the extra term that shows up in the
#'  influence function for aggregated treatment effect parameters
#'  due to estimating the weights
#'
#' @param keepers a vector of indices for which group-time average
#'  treatment effects are used to compute a particular aggregated parameter
#' @param pg a vector with same length as total number of group-time average
#'  treatment effects that contains the probability of being in particular group
#' @param weights additional sampling weights (nx1)
#' @param G vector containing which group a unit belongs to (nx1)
#' @param group vector of groups
#'
#' @return nxk influence function matrix
#'
#' @keywords internal
get_weight_influence <- function(keepers, pg, weights, G, group) {
  # note: weights are all of the form P(G=g|cond)/sum_cond(P(G=g|cond))
  # this is equal to P(G=g)/sum_cond(P(G=g)) which simplifies things here

  # effect of estimating weights in the numerator
  if1 <- sapply(keepers, function(k) {
    (weights * BMisc::TorF(G==group[k]) - pg[k]) /
      sum(pg[keepers])
  })
  # effect of estimating weights in the denominator
  if2 <- base::rowSums( sapply( keepers, function(k) {
    weights * BMisc::TorF(G==group[k]) - pg[k]
  })) %*%
    t(pg[keepers]/(sum(pg[keepers])^2))

  # return the influence function for the weights
  return(if1 - if2)
}

#' @title Get an influence function for particular aggregate parameters
#'
#' @title This is a generic internal function for combining influence
#'  functions across ATT(g,t)'s to return an influence function for
#'  various aggregated treatment effect parameters.
#'
#' @param att vector of group-time average treatment effects
#' @param inf_func influence function for all group-time average treatment effects
#'  (matrix)
#' @param whichones which elements of att will be used to compute the aggregated
#'  treatment effect parameter
#' @param weights_agg the weights to apply to each element of att(whichones);
#'  should have the same dimension as att(whichones)
#' @param wif extra influence function term coming from estimating the weights;
#'  should be n x k matrix where k is dimension of whichones
#'
#' @return nx1 influence function
#'
#' @keywords internal
get_agg_inf_func <- function(att, inf_func, whichones, weights_agg, wif=NULL) {
  # enforce weights are in matrix form
  weights_agg <- as.matrix(weights_agg)

  # multiplies influence function times weights and sums to get vector of weighted IF (of length n)
  agg_inf_func <- inf_func[,whichones] %*% weights_agg

  # Incorporate influence function of the weights
  if (!is.null(wif)) {
    agg_inf_func <- agg_inf_func + wif%*%as.matrix(att[whichones])
  }

  # return influence function
  return(agg_inf_func)
}

#' @title Take influence function and compute standard errors
#'
#' @description Function to take an nx1 influence function and return
#'  a standard error
#'
#' @param influence_function An influence function
#' @param boot a boolean indicating whether bootstrapping was performed
#' @param boot_std_errors a vector of bootstrapped standard errors
#'
#' @return scalar standard error
#'
#' @keywords internal
compute_se_agg <- function(influence_function, boot=FALSE, boot_std_errors = NA) {

  n <- length(influence_function)

  # if we performed boot, boot_std_errors will be a vector of bootstrapped standard errors no null
  if (boot) {
    return(boot_std_errors)
  } else {
    return(sqrt( mean((influence_function)^2)/n ))
  }
}

#' @title Compute Aggregated Treatment Effect Parameters
#'
#' @description Does the heavy lifting on computing aggregated group-time
#'  average treatment effects

#' @param ddd_obj a ddd object (i.e., the results of the [ddd()] function)
#' @param type Which type of aggregated treatment effect parameter to compute.
#'   \code{"simple"} just computes a weighted average of all
#'   group-time average treatment effects with weights proportional to group
#'   size.
#'   \code{"eventstudy"} computes average effects across
#'   different lengths of exposure to the treatment (event times). Here the overall effect averages the effect of the
#'   treatment across the positive lengths of exposure. This is the default option;
#'   \code{"group"} computes average treatment effects across different groups/cohorts; here
#'   the overall effect averages the effect across different groups using group size as weights;
#'   \code{"calendar"} computes average treatment effects across different
#'   time periods, with weights proportional to the group size; here the overall effect averages the effect across each
#'   time period.
#' @param cluster The name of the variable to be used for clustering. The maximum number of cluster variables is 1. Default is \code{NULL}.
#' @param balance_e If set (and if one computes event study), it balances
#'  the sample with respect to event time.  For example, if `balance_e=2`,
#'  `agg_ddd` will drop groups that are not exposed to treatment for
#'  at least three periods, the initial period `e=0` as well as the
#'  next two periods, `e=1` and `e=2`.  This ensures that
#'  the composition of groups does not change when event time changes.
#' @param min_e For event studies, this is the smallest event time to compute
#'  dynamic effects for.  By default, `min_e = -Inf` so that effects at
#'  all lengths of exposure are computed.
#' @param max_e For event studies, this is the largest event time to compute
#'  dynamic effects for.  By default, `max_e = Inf` so that effects at
#'  all lengths of exposure are computed.
#' @param na.rm Logical value if we are to remove missing Values from analyses. Defaults is FALSE.
#' @param boot Boolean for whether or not to compute standard errors using
#'  the multiplier bootstrap.  If standard errors are clustered, then one
#'  must set `boot=TRUE`. Default is value set in the ddd object.  If `boot = FALSE`, then analytical
#'  standard errors are reported.
#' @param nboot The number of bootstrap iterations to use.  The default is the value set in the ddd object,
#'  and this is only applicable if `boot=TRUE`.
#' @param alpha The level of confidence for the confidence intervals.  The default is 0.05. Otherwise, it will
#' use the value set in the ddd object.
#' @param cband Boolean for whether or not to compute a uniform confidence
#'  band that covers all of the group-time average treatment effects
#'  with fixed probability `1 - alpha`.  In order to compute uniform confidence
#'  bands, `boot` must also be set to `TRUE`.  The default is
#'  the value set in the ddd object
#'
#' @return Aggregation object (list) of class [`agg_ddd`]
#'
#' @keywords internal
#'
#' @export
compute_aggregation <- function(ddd_obj,
                                type = "simple",
                                cluster = NULL,
                                balance_e = NULL,
                                min_e = -Inf,
                                max_e = Inf,
                                na.rm = FALSE,
                                boot = FALSE,
                                nboot = NULL,
                                cband = NULL,
                                alpha = 0.05){
  # check if the object is of class `ddd`
  if (!inherits(ddd_obj, "ddd")) {
    stop("Object must be of class `ddd`")
  }

  # check if the object is a multiple period object
  if (!(ddd_obj$argu$multiple_periods)) {
    stop("Object must be a multiple period object")
  }

  # validate types of aggregation
  if (!(type %in% c("simple", "eventstudy", "group", "calendar"))){
    stop("type must be one of 'simple', 'eventstudy', 'group' or 'calendar'")
  }

  # get parameters
  groups <- ddd_obj$groups
  periods <- ddd_obj$periods
  ATT <- ddd_obj$ATT
  inf_func_mat <- ddd_obj$inf_func_mat
  n <- ddd_obj$n
  #data <- ddd_obj$data
  tlist <- ddd_obj$tlist
  glist <- ddd_obj$glist
  dta <- ddd_obj$first_period_dta # data only for first period
  yname <- ddd_obj$argu$yname
  pname <- ddd_obj$argu$pname
  control_group <- ddd_obj$argu$control_group

  # overwriting parameters for multiplier bootstrap if needed:

  # for cluster variable
  if (is.null(cluster)){
    cluster <- ddd_obj$argu$cluster
  }

  # for bootstrap Boolean
  if (is.null(boot)){
    boot <- ddd_obj$argu$boot
  }

  # for number of bootstrap iterations
  if (is.null(nboot)){
    nboot <- ddd_obj$argu$nboot
  }

  # for uniform confidence bands Boolean
  if (is.null(cband)){
    cband <- ddd_obj$argu$cband
  }

  # for alpha level of significance
  if (is.null(alpha)){
    alpha <- ddd_obj$argu$alpha
  }

  # this is useful for summary tables only
  new_argu <- list(cluster = cluster,
                   boot = boot,
                   nboot = nboot,
                   cband = cband,
                   alpha = alpha)

  # flag for boot and cband
  if ((!boot) && (cband)){
    stop("cband is only available when boot = TRUE")
  }

  # recreating a `did_preprocessed` object only with needed parameters to compute bootstrapped standard errors
  did_preprocessed <- list()
  did_preprocessed$preprocessed_data <- dta # we only require data from the first period
  did_preprocessed$cluster <- cluster # cluster variable
  did_preprocessed$nboot <- nboot # number of bootstrap iterations
  did_preprocessed$alpha <- alpha # level of significance
  did_preprocessed$panel <- ddd_obj$argu$panel # panel indicator
  did_preprocessed$allow_unbalanced_panel <- ddd_obj$argu$allow_unbalanced_panel # unbalanced panel indicator


  if((na.rm == FALSE) && base::anyNA(ATT)) stop("Missing values at att_gt found. If you want to remove these, set `na.rm = TRUE'.")

  # removing NA values from ATT(g,t), groups and periods objects
  if(na.rm){
    notna <- !is.na(ATT)
    groups <- groups[notna]
    periods <- periods[notna]
    ATT <- ATT[notna]
    inf_func_mat <- inf_func_mat[, notna]
    glist <- sort(unique(groups))

    # Ensure we have non-missing post-treatment ATTs for each group if the type is "group"
    if(type == "group"){
      # Get the groups that have some non-missing ATT(g,t) in post-treatmemt periods
      gnotna <- sapply(glist, function(g) {
        # look at post-treatment periods for group g
        whichg <- which((groups == g) & (g <= periods))
        attg <- ATT[whichg]
        group_select <- !is.na(mean(attg))
        return(group_select)
      })
      gnotna <- glist[gnotna]
      # indicator for not all post-treatment ATT(g,t) missing
      not_all_na <- groups %in% gnotna
      # Re-do the na.rm thing to update the groups
      groups <- groups[not_all_na]
      periods <- periods[not_all_na]
      ATT <- ATT[not_all_na]
      inf_func_mat <- inf_func_mat[, not_all_na]
      glist <- sort(unique(groups))
    }
  }

  #-----------------------------------------------------------------------------
  # data manipulation
  #-----------------------------------------------------------------------------

  orig_periods <- periods
  orig_group <- groups
  orig_glist <- glist
  orig_tlist <- tlist
  # In case g's are not part of tlist
  orig_gtlist <- sort(unique(c(orig_tlist, orig_glist)))
  uniquet <- seq(1,length(unique(orig_gtlist)))

  # function to switch from "new" t values to  original t values
  t2orig <- function(t) {
    origt <- unique(c(orig_gtlist,0))[which(c(uniquet,0)==t)]
    return(origt)
  }

  # function to switch between "original" t values and new t values
  orig2t <- function(orig) {
    new_t <- c(uniquet,0)[which(unique(c(orig_gtlist,0))==orig)]
    out <- ifelse(length(new_t) == 0, NA, new_t)
    return(out)
  }

  t <- sapply(orig_periods, orig2t)
  group <- sapply(orig_group, orig2t)
  glist <- sapply(orig_glist, orig2t)
  tlist <- unique(t)
  maxT <- max(t)

  weights = dta[["weights"]]

  # compute probability of belonging of each group in glist
  pg <- sapply(orig_glist, function(g) mean(weights*(dta[["first_treat"]]==g)))
  pgg <- pg

  # making pg have the same length of ATT(g,t)
  pg <- pg[match(group, glist)]

  # which group time average treatment effects are post-treatment
  keepers <- which(group <= t & t <= (group + max_e))

  # n x 1 vector of group variable
  G <-  unlist(lapply(dta[["first_treat"]], orig2t))

  #-----------------------------------------------------------------------------
  # Compute the simple ATT summary
  #-----------------------------------------------------------------------------

  if (type == "simple") {

    # simple att
    # averages all post-treatment ATT(g,t) with weights
    # given by group size
    simple.att <- sum(ATT[keepers]*pg[keepers])/(sum(pg[keepers]))
    if(is.nan(simple.att)) simple.att <- NA

    # get the part of the influence function coming from estimated weights
    simple.wif <- get_weight_influence(keepers, pg, weights, G, group)

    # get the overall influence function
    simple.if <- get_agg_inf_func(att=ATT,
                                  inf_func=inf_func_mat,
                                  whichones=keepers,
                                  weights_agg=pg[keepers]/sum(pg[keepers]),
                                  wif=simple.wif)
    # Make it as vector
    simple.if <- as.numeric(simple.if)
    # RUN MULTIPLIER BOOTSTRAP
    if (boot){
      simple_bres <- mboot(inf_func = simple.if, did_preprocessed = did_preprocessed)
      simple_boot_see <- simple_bres$se
    } else {
      simple_boot_see <- NA
    }

    # get standard errors from overall influence function
    simple.se <- compute_se_agg(simple.if, boot, simple_boot_see)
    if(!is.na(simple.se)){
      if(simple.se <= sqrt(.Machine$double.eps)*10) simple.se <- NA
    }


    return(list(overall.att = simple.att,
                overall.se = simple.se,
                type = type,
                yname = yname,
                argu = new_argu,
                pname = pname,
                control_group = control_group,
                inf.function = list(simple.att = simple.if)))
  }

  #-----------------------------------------------------------------------------
  # Compute the group (i.e., selective) treatment timing estimators
  #-----------------------------------------------------------------------------

  if (type == "group") {

    # get group specific ATTs
    # note: there are no estimated weights here
    selective.att.g <- sapply(glist, function(g) {
      # look at post-treatment periods for group g
      whichg <- which( (group == g) & (g <= t) & (t <= (group + max_e))) # added last condition to allow for limit on longest period included in att
      attg <- ATT[whichg]
      mean(attg)
    })
    selective.att.g[is.nan(selective.att.g)] <- NA


    # get standard errors for each group specific ATT
    selective.se.inner <- lapply(glist, function(g) {
      whichg <- which( (group == g) & (g <= t) & (t <= (group + max_e)))  # added last condition to allow for limit on longest period included in att
      inf.func.g <- as.numeric(get_agg_inf_func(att=ATT,
                                                inf_func=inf_func_mat,
                                                whichones=whichg,
                                                weights_agg=pg[whichg]/sum(pg[whichg]),
                                                wif=NULL))
      if (boot){
        g_bres <- mboot(inf_func = inf.func.g, did_preprocessed = did_preprocessed)
        g_boot_se <- g_bres$se
      } else{
        g_boot_se <- NA
      }
      se.g <- compute_se_agg(inf.func.g, boot, g_boot_se)

      list(inf.func=inf.func.g, se=se.g)
    })

    # recover standard errors separately by group
    selective.se.g <- unlist(BMisc::getListElement(selective.se.inner, "se"))
    selective.se.g[selective.se.g <= sqrt(.Machine$double.eps)*10] <- NA

    # recover influence function separately by group
    selective.inf.func.g <- simplify2array(BMisc::getListElement(selective.se.inner, "inf.func"))

    # use multiplier bootstrap (across groups) to get critical value
    # for constructing uniform confidence bands
    selective.crit.val <- stats::qnorm(1 - alpha/2)

    # GET CRITICAL VALUES FOR UNIFORM CONFIDENCE BANDS
    if (cband){
      # if we enter here, it's because we already perform bootstrap. Cannot allow cband without doing bootstrap
      selective.crit.val <- mboot(inf_func = selective.inf.func.g, did_preprocessed = did_preprocessed)$unif_crit_val

      if(is.na(selective.crit.val) | is.infinite(selective.crit.val)){
        warning('Simultaneous critical value is NA. This probably happened because we cannot compute t-statistic (std errors are NA). We then report pointwise conf. intervals.')
        selective.crit.val <- stats::qnorm(1 - alpha/2)
      }

      if(selective.crit.val < stats::qnorm(1 - alpha/2)){
        warning('Simultaneous conf. band is somehow smaller than pointwise one using normal approximation. Since this is unusual, we are reporting pointwise confidence intervals')
        selective.crit.val <- stats::qnorm(1 - alpha/2)
      }

      if(selective.crit.val >= 7){
        warning("Simultaneous critical value is arguably `too large' to be reliable. This usually happens when number of observations per group is small and/or there is no much variation in outcomes.")
      }
    }

    # get overall att under selective treatment timing
    # (here use pgg instead of pg because we can just look at each group)
    selective.att <- sum(selective.att.g * pgg)/sum(pgg)

    # account for having to estimate pgg in the influence function
    selective.wif <- get_weight_influence(keepers=1:length(glist),
                                          pg=pgg,
                                          weights=weights,
                                          G=G,
                                          group=group)

    # get overall influence function
    selective.inf.func <- get_agg_inf_func(att=selective.att.g,
                                           inf_func=selective.inf.func.g,
                                           whichones=(1:length(glist)),
                                           weights_agg=pgg/sum(pgg),
                                           wif=selective.wif)


    selective.inf.func <- as.numeric(selective.inf.func)
    # RUN MULTIPLIER BOOTSTRAP
    if (boot){
      selective_bres <- mboot(inf_func = selective.inf.func, did_preprocessed = did_preprocessed)
      selective_boot_se <- selective_bres$se
    } else{
      selective_boot_se <- NA
    }
    # get overall standard error
    selective.se <- compute_se_agg(selective.inf.func, boot, selective_boot_se)
    if(!is.na(selective.se)){
      if((selective.se <= sqrt(.Machine$double.eps)*10)) selective.se <- NA
    }

    return(list(overall.att=selective.att,
                overall.se=selective.se,
                type=type,
                yname = yname,
                pname = pname,
                control_group = control_group,
                argu = new_argu,
                egt=orig_glist,
                att.egt=selective.att.g,
                se.egt=selective.se.g,
                crit.val.egt=selective.crit.val,
                inf.function = list(selective.inf.func.g = selective.inf.func.g,
                                    selective.inf.func = selective.inf.func)
    ))

  }

  #-----------------------------------------------------------------------------
  # calendar time effects
  #-----------------------------------------------------------------------------

  if (type == "calendar") {

    # drop time periods where no one is treated yet
    # (can't get treatment effects in those periods)
    minG <- min(group)
    calendar.tlist <- tlist[tlist>=minG]

    # calendar time specific atts
    calendar.att.t <- sapply(calendar.tlist, function(t1) {
      # look at post-treatment periods for group g
      whicht <- which( (t == t1) & (group <= t))
      attt <- ATT[whicht]
      pgt <- pg[whicht]/(sum(pg[whicht]))
      sum(pgt * attt)
    })

    # get standard errors and influence functions
    # for each time specific att
    calendar.se.inner <- lapply(calendar.tlist, function(t1) {
      whicht <- which( (t == t1) & (group <= t))
      pgt <- pg[whicht]/(sum(pg[whicht]))
      wif.t <- get_weight_influence(keepers=whicht,
                                    pg=pg,
                                    weights=weights,
                                    G=G,
                                    group=group)
      inf.func.t <- as.numeric(get_agg_inf_func(att=ATT,
                                                inf_func=inf_func_mat,
                                                whichones=whicht,
                                                weights_agg=pgt,
                                                wif=wif.t))
      if (boot){
        t_bres <- mboot(inf_func = inf.func.t, did_preprocessed = did_preprocessed)
        t_boot_se <- t_bres$se
      } else{
        t_boot_se <- NA
      }
      se.t <- compute_se_agg(inf.func.t, boot, t_boot_se)
      list(inf.func=inf.func.t, se=se.t)
    })

    # recover standard errors separately by time
    calendar.se.t <- unlist(BMisc::getListElement(calendar.se.inner, "se"))
    calendar.se.t[calendar.se.t <= sqrt(.Machine$double.eps)*10] <- NA
    # recover influence function separately by time
    calendar.inf.func.t <- simplify2array(BMisc::getListElement(calendar.se.inner, "inf.func"))

    # use multiplier boostrap (across groups) to get critical value
    # for constructing uniform confidence bands
    calendar.crit.val <-  stats::qnorm(1-alpha/2)

    # GET CRITICAL VALUES FOR UNIFORM CONFIDENCE BANDS
    if (cband){
      # if we enter here, it's because we already perform bootstrap. Cannot allow cband without doing bootstrap
      calendar.crit.val <- mboot(inf_func = calendar.inf.func.t, did_preprocessed = did_preprocessed)$unif_crit_val

      if(is.na(calendar.crit.val) | is.infinite(calendar.crit.val)){
        warning('Simultaneous critical value is NA. This probably happened because we cannot compute t-statistic (std errors are NA). We then report pointwise conf. intervals.')
        calendar.crit.val <- stats::qnorm(1 - alpha/2)
      }

      if(calendar.crit.val < stats::qnorm(1 - alpha/2)){
        warning('Simultaneous conf. band is somehow smaller than pointwise one using normal approximation. Since this is unusual, we are reporting pointwise confidence intervals')
        calendar.crit.val <- stats::qnorm(1 - alpha/2)
      }

      if(calendar.crit.val >= 7){
        warning("Simultaneous critical value is arguably `too large' to be reliable. This usually happens when number of observations per group is small and/or there is no much variation in outcomes.")
      }
    }

    # get overall att under calendar time effects
    # this is just average over all time periods
    calendar.att <- mean(calendar.att.t)

    # get overall influence function
    calendar.inf.func <- get_agg_inf_func(att=calendar.att.t,
                                          inf_func=calendar.inf.func.t,
                                          whichones=(1:length(calendar.tlist)),
                                          weights_agg=rep(1/length(calendar.tlist), length(calendar.tlist)),
                                          wif=NULL)
    calendar.inf.func <- as.numeric(calendar.inf.func)
    # RUN MULTIPLIER BOOTSTRAP
    if (boot){
      calendar_bres <- mboot(inf_func = calendar.inf.func, did_preprocessed = did_preprocessed)
      calendar_boot_se <- calendar_bres$se
    } else {
      calendar_boot_se <- NA
    }

    # get overall standard error
    calendar.se <- compute_se_agg(calendar.inf.func, boot, calendar_boot_se)
    if(!is.na(calendar.se)){
      if (calendar.se <= sqrt(.Machine$double.eps)*10) calendar.se <- NA
    }
    return(list(overall.att=calendar.att,
                overall.se=calendar.se,
                type=type,
                yname = yname,
                pname = pname,
                control_group = control_group,
                argu = new_argu,
                egt=sapply(calendar.tlist,t2orig),
                att.egt=calendar.att.t,
                se.egt=calendar.se.t,
                crit.val.egt=calendar.crit.val,
                inf.function = list(calendar.inf.func.t = calendar.inf.func.t,
                                    calendar.inf.func = calendar.inf.func)
    ))
  }

  #-----------------------------------------------------------------------------
  # Compute the event-study estimators
  #-----------------------------------------------------------------------------

  if (type == "eventstudy") {

    # note: event times can be negative here.
    # note: event time = 0 corresponds to "on impact"
    # event times
    eseq <- unique(orig_periods - orig_group)
    eseq <- eseq[order(eseq)]

    # if the user specifies balance_e, then we are going to
    # drop some event times and some groups; if not, we just
    # keep everything (that is what this variable is for)
    include.balanced.gt <- rep(TRUE, length(orig_group))

    # if we balance the sample with respect to event time
    if (!is.null(balance_e)) {
      include.balanced.gt <- (t2orig(maxT) - orig_group >= balance_e)

      eseq <- unique(orig_periods[include.balanced.gt] - orig_group[include.balanced.gt])
      eseq <- eseq[order(eseq)]

      eseq <- eseq[ (eseq <= balance_e) & (eseq >= balance_e - t2orig(maxT) + t2orig(1))]

    }

    # only looks at some event times
    eseq <- eseq[ (eseq >= min_e) & (eseq <= max_e) ]

    # compute atts that are specific to each event time
    dynamic.att.e <- sapply(eseq, function(e) {
      # keep att(g,t) for the right g&t as well as ones that
      # are not trimmed out from balancing the sample
      whiche <- which( (orig_periods - orig_group == e) & (include.balanced.gt) )
      atte <- ATT[whiche]
      pge <- pg[whiche]/(sum(pg[whiche]))
      sum(atte*pge)
    })

    # compute standard errors for dynamic effects
    dynamic.se.inner <- lapply(eseq, function(e) {
      whiche <- which( (orig_periods - orig_group == e) & (include.balanced.gt) )
      pge <- pg[whiche]/(sum(pg[whiche]))
      wif.e <- get_weight_influence(whiche, pg, weights, G, group)
      inf.func.e <- as.numeric(get_agg_inf_func(att=ATT,
                                                inf_func=inf_func_mat,
                                                whichones=whiche,
                                                weights_agg=pge,
                                                wif=wif.e))
      if (boot){
        e_bres <- mboot(inf_func = inf.func.e, did_preprocessed = did_preprocessed)
        e_boot_se <- e_bres$se
      } else{
        e_boot_se <- NA
      }
      se.e <- compute_se_agg(inf.func.e, boot, e_boot_se)
      list(inf.func=inf.func.e, se=se.e)
    })

    dynamic.se.e <- unlist(BMisc::getListElement(dynamic.se.inner, "se"))
    dynamic.se.e[dynamic.se.e <= sqrt(.Machine$double.eps)*10] <- NA

    dynamic.inf.func.e <- simplify2array(BMisc::getListElement(dynamic.se.inner, "inf.func"))

    dynamic.crit.val <- stats::qnorm(1 - alpha/2)

    # GET CRITICAL VALUES FOR UNIFORM CONFIDENCE BANDS
    if (cband){
      # if we enter here, it's because we already perform bootstrap. Cannot allow cband without doing bootstrap
      dynamic.crit.val <- mboot(inf_func = dynamic.inf.func.e, did_preprocessed = did_preprocessed)$unif_crit_val

      if(is.na(dynamic.crit.val) | is.infinite(dynamic.crit.val)){
        warning('Simultaneous critical value is NA. This probably happened because we cannot compute t-statistic (std errors are NA). We then report pointwise conf. intervals.')
        dynamic.crit.val <- stats::qnorm(1 - alpha/2)
      }

      if(dynamic.crit.val < stats::qnorm(1 - alpha/2)){
        warning('Simultaneous conf. band is somehow smaller than pointwise one using normal approximation. Since this is unusual, we are reporting pointwise confidence intervals')
        dynamic.crit.val <- stats::qnorm(1 - alpha/2)
      }

      if(dynamic.crit.val >= 7){
        warning("Simultaneous critical value is arguably `too large' to be reliable. This usually happens when number of observations per group is small and/or there is no much variation in outcomes.")
      }
    }


    # get overall average treatment effect
    # by averaging over positive dynamics
    epos <- eseq >= 0
    dynamic.att <- mean(dynamic.att.e[epos])
    dynamic.inf.func <- get_agg_inf_func(att=dynamic.att.e[epos],
                                         inf_func=as.matrix(dynamic.inf.func.e[,epos]),
                                         whichones=(1:sum(epos)),
                                         weights_agg=(rep(1/sum(epos), sum(epos))),
                                         wif=NULL)

    dynamic.inf.func <- as.numeric(dynamic.inf.func)
    # RUN MULTIPLIER BOOTSTRAP
    if (boot){
      dynamic_bres <- mboot(inf_func = dynamic.inf.func, did_preprocessed = did_preprocessed)
      dynamic_boot_se <- dynamic_bres$se
    } else {
      dynamic_boot_se <- NA
    }
    dynamic.se <- compute_se_agg(dynamic.inf.func, boot, dynamic_boot_se)
    if(!is.na(dynamic.se)){
      if (dynamic.se <= sqrt(.Machine$double.eps)*10) dynamic.se <- NA
    }

    return(list(overall.att=dynamic.att,
                overall.se=dynamic.se,
                type=type,
                yname = yname,
                pname = pname,
                control_group = control_group,
                argu = new_argu,
                egt=eseq,
                att.egt=dynamic.att.e,
                se.egt=dynamic.se.e,
                crit.val.egt=dynamic.crit.val,
                inf.function = list(dynamic.inf.func.e = dynamic.inf.func.e,
                                    dynamic.inf.func = dynamic.inf.func),
                min_e=min_e,
                max_e=max_e,
                balance_e=balance_e
    ))
  }

}
