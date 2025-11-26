#' Utility functions to compute nuisances parameters.
#' @importFrom stats update
#' @import parglm
#' @import speedglm
#' @import data.table
#' @noRd
#--------------------------------------------------

# utility function to generate equation to estimate pscore
get_formula_pscore <- function(xformula, weights = TRUE){
  # Get a formula object for pscore estimation
  if (weights) {
    formula_obj <- stats::update(xformula, PA4 ~ . + weights)
  } else {
    formula_obj <- stats::update(xformula, PA4 ~ .)
  }
  return(formula_obj)
}

# utility function to generate equation to estimate outcome regression
get_formula_reg <- function(xformula, weights = TRUE){
  # Get a formula object for outcome regression
  if (weights) {
    formula_obj <- stats::update(xformula, deltaY ~ . + weights)
  } else {
    formula_obj <- stats::update(xformula, deltaY ~ .)
  }
  return(formula_obj)
}

# Function to compute propensity scores using parglm for multiple subgroups
compute_pscore <- function(data, condition_subgroup, xformula) {
  # get formula for pscore estimation using covariates
  formula_pscore <- get_formula_pscore(xformula, weights = FALSE)
  # Subset data for condition_subgroup and subgroup == 4 or the given condition_subgroup
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]
  # Adding treatment variable P(1{PA4 = 4}|X)
  condition_data[, "PA4" := ifelse(condition_data$subgroup == 4, 1, 0)]
  uid_condition_data <- unique(condition_data, by = "id")

  # Fit logistic regression model using parglm
  model <- withCallingHandlers(
    parglm::parglm(
      formula_pscore,
      data     = uid_condition_data,
      family   = stats::binomial(),
      weights  = weights,
      control  = parglm.control(nthreads = data.table::getDTthreads()),
      intercept = FALSE
    ),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("^Too few observation compared to the number of threads", msg)) {
        invokeRestart("muffleWarning")
      }
    }
  )

  # Flag for convergence of glm
  if (!model$converged) {
    warning(paste("Logistic regression model for subgroup", condition_subgroup,
                  "did not converge."))
  }

  # Flag for weird coefficients
  if(anyNA(model$coefficients)) {
    stop(paste("Pscore model coefficients for subgroup", condition_subgroup,
               "have NA components. Multicollinearity of covariates is a likely reason"))
  }

  # Compute propensity scores
  propensity_scores <- predict(model, newdata = uid_condition_data, type = "response")

  # Warning for overlap condition
  if (any(propensity_scores < 5e-4)) {
    warning(paste("Propensity scores for comparison subgroup", condition_subgroup,
                  "have poor overlap."))
  }

  # Avoid divide by zero
  propensity_scores <- pmin(propensity_scores, 1 - 1e-16)

  # Save Hessian matrix for influence function
  hessian_matrix <- stats::vcov(model) * nrow(uid_condition_data)

  return(list(propensity_scores = propensity_scores, hessian_matrix = hessian_matrix))
}

compute_pscore_null <- function(data, condition_subgroup) {
  # This is valid when REG is used
  # Subset data for condition_subgroup and subgroup == 4 or the given condition_subgroup
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]
  uid_condition_data <- unique(condition_data, by = "id")

  propensity_scores <- rep(1, nrow(uid_condition_data))
  hessian_matrix <- NA

  return(list(propensity_scores = propensity_scores, hessian_matrix = hessian_matrix))
}

# Function to compute outcome regression for multiple subgroups
compute_outcome_regression <- function(data, condition_subgroup, xformula){
  # Subset data for condition_subgroup and subgroup == 4 or the given condition_subgroup
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]
  # Subset data for the control group (subgroup != 4)
  control_data <- condition_data[condition_data$subgroup == condition_subgroup]
  y1_control = control_data[control_data$post == 1, y]
  y0_control = control_data[control_data$post == 0, y]
  # generate deltaY for control group
  deltaY_control = y1_control - y0_control
  # get covariates including the intercept (conditioning in pre-treatment period since covariates are invariant)
  cov_control <- stats::model.matrix(as.formula(xformula), data = control_data[control_data$post == 0])
  # get weights (conditioning in pre-treatment period since weights are invariant)
  i_weights = control_data[control_data$post == 0, weights]

  # get coefficients of outcome regression model for the control group (subgroup != 4)
  # Attempt to compute reg.coeff using speedglm
  try_speedglm <- tryCatch({
    reg.coeff <- stats::coef(speedglm::speedglm.wfit(y = deltaY_control,
                                                     X = as.matrix(cov_control),
                                                     intercept = FALSE,
                                                     weights = i_weights))
  }, error = function(e) {
    # If error, attempt to compute reg.coeff using lm.wfit
    reg.coeff <- tryCatch({
      stats::coef(stats::lm.wfit(x = as.matrix(cov_control),
                                 y = deltaY_control,
                                 w = i_weights))
    }, error = function(e2) {
      # If error with lm.wfit, stop the program and send an error message
      stop("Error in computing regression coefficients: Subgroup ", condition_subgroup, " may have insufficient data.")
    })
  })

  # Flag for NA coefficients
  if(anyNA(reg.coeff)) {
    stop(paste("Outcome regression model coefficients for subgroup", condition_subgroup,
               "have NA components. Multicollinearity of covariates is a likely reason"))
  }

  # compute regression adjustment
  y1 = condition_data[condition_data$post == 1, y]
  y0 = condition_data[condition_data$post == 0, y]
  deltaY = y1 - y0
  # get covariates including the intercept (conditioning in pre-treatment period since covariates are invariant)
  covX <- stats::model.matrix(xformula, data = condition_data[condition_data$post == 0])
  # compute OR delta
  or_delta = as.vector(tcrossprod(reg.coeff, as.matrix(covX)))
  # compute regression adjustment
  # reg_adjust = deltaY - or_delta
  return(list(deltaY = deltaY, or_delta = or_delta))
}

compute_outcome_regression_null <- function(data, condition_subgroup){
  # This is valid when IPW is used
  # Subset data for condition_subgroup and subgroup == 4 or the given condition_subgroup
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]

  # compute regression adjustment
  y1 = condition_data[condition_data$post == 1, y]
  y0 = condition_data[condition_data$post == 0, y]
  deltaY = y1 - y0

  # Filter the data for the pre-treatment period
  pre_treatment_data <- condition_data[condition_data$post == 0]
  # Create a zero vector with the same length as the number of rows in pre_treatment_data
  or_delta <- numeric(nrow(pre_treatment_data))

  return(list(deltaY = deltaY, or_delta = or_delta))
}

# Function to compute the average treatment effect for multiple subgroups
compute_did <- function(data, condition_subgroup, pscores, reg_adjustment, xformula, est_method){

  data <- unique(data, by = "id")
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]
  PA4 = ifelse(condition_data$subgroup == 4, 1, 0)
  PAa = ifelse(condition_data$subgroup == condition_subgroup, 1, 0)

  # Compute propensity scores
  if (condition_subgroup == 3) {
    pscore <- pscores[[1]]$propensity_scores
    hessian <- pscores[[1]]$hessian_matrix
    deltaY <- reg_adjustment[[1]]$deltaY
    or_delta <- reg_adjustment[[1]]$or_delta

  } else if (condition_subgroup == 2) {
    pscore <- pscores[[2]]$propensity_scores
    hessian <- pscores[[2]]$hessian_matrix
    deltaY <- reg_adjustment[[2]]$deltaY
    or_delta <- reg_adjustment[[2]]$or_delta

  } else if (condition_subgroup == 1) {
    pscore <- pscores[[3]]$propensity_scores
    hessian <- pscores[[3]]$hessian_matrix
    deltaY <- reg_adjustment[[3]]$deltaY
    or_delta <- reg_adjustment[[3]]$or_delta

  } else {
    stop("Invalid condition_subgroup")
  }

  # get weights
  # i_weights = condition_data[condition_data$post == 0, weights]
  i_weights = condition_data$weights

  ################################
  # Get doubly-robust estimation #
  ################################
  w_treat = i_weights * PA4
  if (est_method == "reg") {
    # Compute doubly-robust estimation
    w_control = (i_weights * PAa)
  } else {
    w_control = (i_weights * pscore * PAa) / (1 - pscore)
  }
  riesz_treat = w_treat * (deltaY - or_delta)
  riesz_control = w_control * (deltaY - or_delta)
  att_treat = mean(riesz_treat, na.rm = TRUE) / mean(w_treat, na.rm = TRUE)
  att_control = mean(riesz_control, na.rm = TRUE) / mean(w_control, na.rm = TRUE)

  dr_att = att_treat - att_control

  ##########################
  # Get influence function #
  ##########################

  # Influence function related to the estimation of pscores
  covX = stats::model.matrix(xformula, data = condition_data)
  if (est_method == "reg") {
    inf_control_pscore <- 0
  } else {
    M2 <- base::colMeans(w_control * (deltaY - or_delta - att_control) * covX, na.rm = TRUE) # reg_adjust = deltaY - m_delta(x)
    score_ps <- i_weights * (PA4 - pscore) * covX
    # asymptotic linear representation of logit's beta
    score_ps_no_na <- na.omit(score_ps) # Exclude rows with NA values
    asy_lin_rep_ps <- score_ps_no_na %*% hessian
    inf_control_pscore <- asy_lin_rep_ps %*% as.matrix(M2)
  }

  if (est_method == "ipw") {
    inf_treat_or <- 0
    inf_cont_or <- 0
  } else {
    # Influence function related to the estimation of outcome model
    M1 <- base::colMeans(w_treat * covX, na.rm = TRUE)
    M3 <- base::colMeans(w_control * covX, na.rm = TRUE)

    # Influence function related to the estimation of regression model
    or_x <- i_weights * PAa * covX
    or_ex <- i_weights * PAa * (deltaY - or_delta) * covX
    XpX <- crossprod(or_x, covX)/nrow(condition_data)

    #asymptotic linear representation of the beta
    asy_linear_or <- t(solve(XpX, t(or_ex)))

    #or for treat
    inf_treat_or <- -asy_linear_or %*% M1
    #or for control
    inf_cont_or <- -asy_linear_or %*% M3
  }

  # Influence function from did
  inf_control_did <- riesz_control - w_control*att_control
  inf_treat_did <- riesz_treat - w_treat*att_treat

  # Influence function for the ATT
  inf_control <- (inf_control_did + inf_control_pscore + inf_cont_or)/mean(w_control)
  inf_treat <- (inf_treat_did + inf_treat_or)/mean(w_treat)
  # putting all together
  inf_func <- inf_treat - inf_control

  # fill zeros in influence function for observation outside the subgroup analyzed
  data[, "inf_func_result" := numeric(.N)]
  data[data$subgroup %in% c(condition_subgroup, 4), "inf_func_result" := inf_func]

  return(list(dr_att = dr_att, inf_func = data[["inf_func_result"]]))
}


#--------------------------------------------------
# Repeated Cross-Section Functions
#--------------------------------------------------

compute_pscore_rc <- function(data, condition_subgroup, xformula) {
  # Similar to compute_pscore but no unique(by="id")
  formula_pscore <- get_formula_pscore(xformula, weights = FALSE)
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]
  condition_data[, "PA4" := ifelse(condition_data$subgroup == 4, 1, 0)]
  # Use all data for RC propensity score (assumes pooling)
  
  model <- withCallingHandlers(
    parglm::parglm(
      formula_pscore,
      data     = condition_data,
      family   = stats::binomial(),
      weights  = weights,
      control  = parglm.control(nthreads = data.table::getDTthreads()),
      intercept = FALSE
    ),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("^Too few observation compared to the number of threads", msg)) {
        invokeRestart("muffleWarning")
      }
    }
  )

  if (!model$converged) warning(paste("RC pscore model for subgroup", condition_subgroup, "did not converge."))
  if(anyNA(model$coefficients)) stop(paste("RC pscore model coefficients for subgroup", condition_subgroup, "have NA components. Multicollinearity of covariates is a likely reason. "))

  propensity_scores <- predict(model, newdata = condition_data, type = "response")
  
  if (any(propensity_scores < 5e-8)) warning(paste("Propensity scores for comparison subgroup", condition_subgroup, "have poor overlap."))
  propensity_scores <- pmin(propensity_scores, 1 - 1e-16)

  hessian_matrix <- stats::vcov(model) * nrow(condition_data)

  return(list(propensity_scores = propensity_scores, hessian_matrix = hessian_matrix))
}

compute_pscore_null_rc <- function(data, condition_subgroup) {
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]
  propensity_scores <- rep(1, nrow(condition_data))
  hessian_matrix <- NA
  return(list(propensity_scores = propensity_scores, hessian_matrix = hessian_matrix))
}

compute_outcome_regression_rc <- function(data, condition_subgroup, xformula){
  # Subset data for condition_subgroup and subgroup == 4 or the given condition_subgroup
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]
  # compute outcome regressioin model for each of the four cells: (subgroup=4, post=1), (subgroup=4, post=0), (subgroup=condition_subgroup, post=1), (subgroup=condition_subgroup, post=0)
  fit_predict_mu <- function(subg, pst) {
     subset_data <- condition_data[condition_data$subgroup == subg & condition_data$post == pst]
     if(nrow(subset_data) == 0) return(rep(NA, nrow(condition_data)))
     
     # Covariates including intercept for fitting
     cov_fit <- stats::model.matrix(as.formula(xformula), data = subset_data)
     y_fit <- subset_data$y
     w_fit <- subset_data$weights
     
     # Fit
     reg.coeff <- tryCatch({
        stats::coef(speedglm::speedglm.wfit(y = y_fit, X = as.matrix(cov_fit), intercept = FALSE, weights = w_fit))
     }, error = function(e) {
        stats::coef(stats::lm.wfit(x = as.matrix(cov_fit), y = y_fit, w = w_fit))
     })
     
     if(anyNA(reg.coeff)) stop(paste("Outcome regression coefficients NA for subgroup", subg, "post", pst))
     
     # Predict for ALL
     cov_all <- stats::model.matrix(as.formula(xformula), data = condition_data)
     as.vector(tcrossprod(reg.coeff, as.matrix(cov_all)))
  }

  or_trt_post <- fit_predict_mu(4, 1) # OR for: S=g, Q=1, T=post
  or_trt_pre <- fit_predict_mu(4, 0) # OR for: S=g, Q=1, T=pre
  or_ctrl_post <- fit_predict_mu(condition_subgroup, 1) # OR for: S=g', Q=q', T=post
  or_ctrl_pre <- fit_predict_mu(condition_subgroup, 0) # OR for: S=g', Q=q', T=pre


  return(list(or_trt_post=or_trt_post, or_trt_pre=or_trt_pre, or_ctrl_post=or_ctrl_post, or_ctrl_pre=or_ctrl_pre))
}

compute_outcome_regression_null_rc <- function(data, condition_subgroup){
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]
  n <- nrow(condition_data)
  return(list(or_trt_post=rep(0,n), or_trt_pre=rep(0,n), or_ctrl_post=rep(0,n), or_ctrl_pre=rep(0,n)))
}

compute_did_rc <- function(data, condition_subgroup, pscores, reg_adjustment, xformula, est_method){
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]

  if (condition_subgroup == 3) idx <- 1 else if (condition_subgroup == 2) idx <- 2 else idx <- 3

  pscore <- pscores[[idx]]$propensity_scores
  hessian <- pscores[[idx]]$hessian_matrix
  mu <- reg_adjustment[[idx]]

  # Basic variables
  n <- nrow(condition_data)
  post <- condition_data$post
  PA4 <- ifelse(condition_data$subgroup == 4, 1, 0) # PA4=1 for treated (subgroup 4), PA4=0 for control
  PAa <- ifelse(condition_data$subgroup == condition_subgroup, 1, 0) # PAa=1 for control subgroup
  y <- condition_data$y
  i_weights <- condition_data$weights

  # Get covariate matrix
  covX <- stats::model.matrix(xformula, data = condition_data)

  #############################################
  # Method-specific estimation
  #############################################

  if (est_method == "ipw") {
    #-----------------------------------------
    # IPW Estimator
    #-----------------------------------------
    # Riesz representers
    riesz_treat_pre <- i_weights * PA4 * (1 - post)
    riesz_treat_post <- i_weights * PA4 * post
    riesz_control_pre <- i_weights * pscore * PAa * (1 - post) / (1 - pscore)
    riesz_control_post <- i_weights * pscore * PAa * post / (1 - pscore)

    # Elements of the influence function (summands)
    eta_treat_pre <- riesz_treat_pre * y / mean(riesz_treat_pre)
    eta_treat_post <- riesz_treat_post * y / mean(riesz_treat_post)
    eta_control_pre <- riesz_control_pre * y / mean(riesz_control_pre)
    eta_control_post <- riesz_control_post * y / mean(riesz_control_post)

    # Estimator of each component
    att_treat_pre <- mean(eta_treat_pre)
    att_treat_post <- mean(eta_treat_post)
    att_control_pre <- mean(eta_control_pre)
    att_control_post <- mean(eta_control_post)

    # ATT estimator
    att <- (att_treat_post - att_treat_pre) - (att_control_post - att_control_pre)

    # Influence function
    # Asymptotic linear representation of logit's beta's
    W <- pscore * (1 - pscore) * i_weights
    score_ps <- i_weights * (PA4 - pscore) * covX
    asy_lin_rep_ps <- score_ps %*% hessian

    # Leading term of the influence function
    inf_treat_pre <- eta_treat_pre - riesz_treat_pre * att_treat_pre / mean(riesz_treat_pre)
    inf_treat_post <- eta_treat_post - riesz_treat_post * att_treat_post / mean(riesz_treat_post)
    inf_treat <- inf_treat_post - inf_treat_pre

    inf_control_pre <- eta_control_pre - riesz_control_pre * att_control_pre / mean(riesz_control_pre)
    inf_control_post <- eta_control_post - riesz_control_post * att_control_post / mean(riesz_control_post)
    inf_control <- inf_control_post - inf_control_pre

    # Estimation effect from gamma hat (pscore)
    M2_pre <- base::colMeans(riesz_control_pre * (y - att_control_pre) * covX) / mean(riesz_control_pre)
    M2_post <- base::colMeans(riesz_control_post * (y - att_control_post) * covX) / mean(riesz_control_post)
    inf_control_ps <- asy_lin_rep_ps %*% (M2_post - M2_pre)

    # Final influence function
    inf_control <- inf_control + inf_control_ps
    att_inf_func <- inf_treat - inf_control

  } else if (est_method == "reg") {
    #-----------------------------------------
    # REG Estimator
    #-----------------------------------------
    or_ctrl_pre <- mu$or_ctrl_pre
    or_ctrl_post <- mu$or_ctrl_post

    # Riesz representers
    riesz_treat_pre <- i_weights * PA4 * (1 - post)
    riesz_treat_post <- i_weights * PA4 * post
    riesz_control <- i_weights * PA4

    reg_att_treat_pre <- riesz_treat_pre * y
    reg_att_treat_post <- riesz_treat_post * y
    reg_att_control <- riesz_control * (or_ctrl_post - or_ctrl_pre)

    eta_treat_pre <- mean(reg_att_treat_pre) / mean(riesz_treat_pre)
    eta_treat_post <- mean(reg_att_treat_post) / mean(riesz_treat_post)
    eta_control <- mean(reg_att_control) / mean(riesz_control)

    att <- (eta_treat_post - eta_treat_pre) - eta_control

    # Influence function
    # Asymptotic linear representation of OLS parameters
    # Pre-period
    weights_ols_pre <- i_weights * PAa * (1 - post)
    wols_x_pre <- weights_ols_pre * covX
    wols_eX_pre <- weights_ols_pre * (y - or_ctrl_pre) * covX
    XpX_pre <- base::crossprod(wols_x_pre, covX) / n
    XpX_inv_pre <- solve(XpX_pre)
    asy_lin_rep_ols_pre <- wols_eX_pre %*% XpX_inv_pre

    # Post-period
    weights_ols_post <- i_weights * PAa * post
    wols_x_post <- weights_ols_post * covX
    wols_eX_post <- weights_ols_post * (y - or_ctrl_post) * covX
    XpX_post <- base::crossprod(wols_x_post, covX) / n
    XpX_inv_post <- solve(XpX_post)
    asy_lin_rep_ols_post <- wols_eX_post %*% XpX_inv_post

    # Leading term
    inf_treat_pre <- (reg_att_treat_pre - riesz_treat_pre * eta_treat_pre) / mean(riesz_treat_pre)
    inf_treat_post <- (reg_att_treat_post - riesz_treat_post * eta_treat_post) / mean(riesz_treat_post)
    inf_treat <- inf_treat_post - inf_treat_pre

    # Control component
    inf_control_1 <- (reg_att_control - riesz_control * eta_control)
    M1 <- base::colMeans(riesz_control * covX)
    inf_control_2_post <- asy_lin_rep_ols_post %*% M1
    inf_control_2_pre <- asy_lin_rep_ols_pre %*% M1
    inf_control <- (inf_control_1 + inf_control_2_post - inf_control_2_pre) / mean(riesz_control)

    att_inf_func <- inf_treat - inf_control

  } else if (est_method == "dr") {
    #-----------------------------------------
    # DR Estimator
    #-----------------------------------------
    or_ctrl_pre <- mu$or_ctrl_pre
    or_ctrl_post <- mu$or_ctrl_post
    or_ctrl <- post * or_ctrl_post + (1 - post) * or_ctrl_pre

    or_trt_pre <- mu$or_trt_pre
    or_trt_post <- mu$or_trt_post

    # Riesz representers
    riesz_treat_pre <- i_weights * PA4 * (1 - post)
    riesz_treat_post <- i_weights * PA4 * post
    riesz_control_pre <- i_weights * pscore * PAa * (1 - post) / (1 - pscore)
    riesz_control_post <- i_weights * pscore * PAa * post / (1 - pscore)

    riesz_d <- i_weights * PA4
    riesz_dt1 <- i_weights * PA4 * post
    riesz_dt0 <- i_weights * PA4 * (1 - post)

    # Elements of the influence function (summands)
    eta_treat_pre <- riesz_treat_pre * (y - or_ctrl) / mean(riesz_treat_pre)
    eta_treat_post <- riesz_treat_post * (y - or_ctrl) / mean(riesz_treat_post)
    eta_control_pre <- riesz_control_pre * (y - or_ctrl) / mean(riesz_control_pre)
    eta_control_post <- riesz_control_post * (y - or_ctrl) / mean(riesz_control_post)

    # Extra elements for the locally efficient DRDID
    eta_d_post <- riesz_d * (or_trt_post - or_ctrl_post) / mean(riesz_d)
    eta_dt1_post <- riesz_dt1 * (or_trt_post - or_ctrl_post) / mean(riesz_dt1)
    eta_d_pre <- riesz_d * (or_trt_pre - or_ctrl_pre) / mean(riesz_d)
    eta_dt0_pre <- riesz_dt0 * (or_trt_pre - or_ctrl_pre) / mean(riesz_dt0)

    # Estimator of each component
    att_treat_pre <- mean(eta_treat_pre)
    att_treat_post <- mean(eta_treat_post)
    att_control_pre <- mean(eta_control_pre)
    att_control_post <- mean(eta_control_post)

    att_d_post <- mean(eta_d_post)
    att_dt1_post <- mean(eta_dt1_post)
    att_d_pre <- mean(eta_d_pre)
    att_dt0_pre <- mean(eta_dt0_pre)

    # ATT estimator
    att <- (att_treat_post - att_treat_pre) - (att_control_post - att_control_pre) +
      (att_d_post - att_dt1_post) - (att_d_pre - att_dt0_pre)

    # Influence function
    # Asymptotic linear representation of OLS parameters
    # Control pre-period
    weights_ols_pre <- i_weights * PAa * (1 - post)
    wols_x_pre <- weights_ols_pre * covX
    wols_eX_pre <- weights_ols_pre * (y - or_ctrl_pre) * covX
    XpX_pre <- base::crossprod(wols_x_pre, covX) / n
    XpX_inv_pre <- solve(XpX_pre)
    asy_lin_rep_ols_pre <- wols_eX_pre %*% XpX_inv_pre

    # Control post-period
    weights_ols_post <- i_weights * PAa * post
    wols_x_post <- weights_ols_post * covX
    wols_eX_post <- weights_ols_post * (y - or_ctrl_post) * covX
    XpX_post <- base::crossprod(wols_x_post, covX) / n
    XpX_inv_post <- solve(XpX_post)
    asy_lin_rep_ols_post <- wols_eX_post %*% XpX_inv_post

    # Treated pre-period
    weights_ols_pre_treat <- i_weights * PA4 * (1 - post)
    wols_x_pre_treat <- weights_ols_pre_treat * covX
    wols_eX_pre_treat <- weights_ols_pre_treat * (y - or_trt_pre) * covX
    XpX_pre_treat <- base::crossprod(wols_x_pre_treat, covX) / n
    XpX_inv_pre_treat <- solve(XpX_pre_treat)
    asy_lin_rep_ols_pre_treat <- wols_eX_pre_treat %*% XpX_inv_pre_treat

    # Treated post-period
    weights_ols_post_treat <- i_weights * PA4 * post
    wols_x_post_treat <- weights_ols_post_treat * covX
    wols_eX_post_treat <- weights_ols_post_treat * (y - or_trt_post) * covX
    XpX_post_treat <- base::crossprod(wols_x_post_treat, covX) / n
    XpX_inv_post_treat <- solve(XpX_post_treat)
    asy_lin_rep_ols_post_treat <- wols_eX_post_treat %*% XpX_inv_post_treat

    # Asymptotic linear representation of logit's beta's
    W <- pscore * (1 - pscore) * i_weights
    score_ps <- i_weights * (PA4 - pscore) * covX
    asy_lin_rep_ps <- score_ps %*% hessian

    # Influence function of the "treat" component
    inf_treat_pre <- eta_treat_pre - riesz_treat_pre * att_treat_pre / mean(riesz_treat_pre)
    inf_treat_post <- eta_treat_post - riesz_treat_post * att_treat_post / mean(riesz_treat_post)

    # Estimation effect from beta hat
    M1_post <- -base::colMeans(riesz_treat_post * post * covX) / mean(riesz_treat_post)
    M1_pre <- -base::colMeans(riesz_treat_pre * (1 - post) * covX) / mean(riesz_treat_pre)

    inf_treat_or_post <- asy_lin_rep_ols_post %*% M1_post
    inf_treat_or_pre <- asy_lin_rep_ols_pre %*% M1_pre

    # Influence function of control component
    inf_control_pre <- eta_control_pre - riesz_control_pre * att_control_pre / mean(riesz_control_pre)
    inf_control_post <- eta_control_post - riesz_control_post * att_control_post / mean(riesz_control_post)

    # Estimation effect from gamma hat (pscore)
    M2_pre <- base::colMeans(riesz_control_pre * (y - or_ctrl - att_control_pre) * covX) / mean(riesz_control_pre)
    M2_post <- base::colMeans(riesz_control_post * (y - or_ctrl - att_control_post) * covX) / mean(riesz_control_post)
    inf_control_ps <- asy_lin_rep_ps %*% (M2_post - M2_pre)

    # Estimation effect from beta hat
    M3_post <- -base::colMeans(riesz_control_post * post * covX) / mean(riesz_control_post)
    M3_pre <- -base::colMeans(riesz_control_pre * (1 - post) * covX) / mean(riesz_control_pre)

    inf_control_or_post <- asy_lin_rep_ols_post %*% M3_post
    inf_control_or_pre <- asy_lin_rep_ols_pre %*% M3_pre

    # Adjustment terms
    inf_eff1 <- eta_d_post - riesz_d * att_d_post / mean(riesz_d)
    inf_eff2 <- eta_dt1_post - riesz_dt1 * att_dt1_post / mean(riesz_dt1)
    inf_eff3 <- eta_d_pre - riesz_d * att_d_pre / mean(riesz_d)
    inf_eff4 <- eta_dt0_pre - riesz_dt0 * att_dt0_pre / mean(riesz_dt0)
    inf_eff <- (inf_eff1 - inf_eff2) - (inf_eff3 - inf_eff4)

    # Estimation effect of OR coefficients
    mom_post <- base::colMeans((riesz_d / mean(riesz_d) - riesz_dt1 / mean(riesz_dt1)) * covX)
    mom_pre <- base::colMeans((riesz_d / mean(riesz_d) - riesz_dt0 / mean(riesz_dt0)) * covX)
    inf_or_post <- (asy_lin_rep_ols_post_treat - asy_lin_rep_ols_post) %*% mom_post
    inf_or_pre <- (asy_lin_rep_ols_pre_treat - asy_lin_rep_ols_pre) %*% mom_pre

    # Put all pieces together
    inf_treat_or <- inf_treat_or_post + inf_treat_or_pre
    inf_control_or <- inf_control_or_post + inf_control_or_pre
    inf_or <- inf_or_post - inf_or_pre

    inf_treat <- inf_treat_post - inf_treat_pre + inf_treat_or
    inf_control <- inf_control_post - inf_control_pre + inf_control_ps + inf_control_or

    # Final influence function
    att_inf_func <- inf_treat - inf_control + inf_eff + inf_or

  } else {
    stop("Invalid est_method. Must be 'ipw', 'reg', or 'dr'")
  }

  # Fill zeros in influence function for observations outside the subgroup analyzed
  data[, "inf_func_result" := numeric(.N)]
  data[data$subgroup %in% c(condition_subgroup, 4), "inf_func_result" := att_inf_func]

  return(list(dr_att = att, inf_func = data[["inf_func_result"]]))
}
