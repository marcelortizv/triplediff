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

  if (!model$converged) warning(paste("RC Pscore model for subgroup", condition_subgroup, "did not converge."))
  if(anyNA(model$coefficients)) stop(paste("RC Pscore model coefficients for subgroup", condition_subgroup, "have NA components. Multicollinearity of covariates is a likely reason. "))

  propensity_scores <- predict(model, newdata = condition_data, type = "response")
  
  if (any(propensity_scores < 5e-4)) warning(paste("Propensity scores for comparison subgroup", condition_subgroup, "have poor overlap."))
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

  mu_11 <- fit_predict_mu(4, 1)
  mu_10 <- fit_predict_mu(4, 0)
  mu_01 <- fit_predict_mu(condition_subgroup, 1)
  mu_00 <- fit_predict_mu(condition_subgroup, 0)

  return(list(mu_11=mu_11, mu_10=mu_10, mu_01=mu_01, mu_00=mu_00))
}

compute_outcome_regression_null_rc <- function(data, condition_subgroup){
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]
  n <- nrow(condition_data)
  return(list(mu_11=rep(0,n), mu_10=rep(0,n), mu_01=rep(0,n), mu_00=rep(0,n)))
}

compute_did_rc <- function(data, condition_subgroup, pscores, reg_adjustment, xformula, est_method){
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]
  
  if (condition_subgroup == 3) idx <- 1 else if (condition_subgroup == 2) idx <- 2 else idx <- 3
  
  pscore <- pscores[[idx]]$propensity_scores
  mu <- reg_adjustment[[idx]]
  
  post <- condition_data$post
  D <- ifelse(condition_data$subgroup == 4, 1, 0) # Treat=1, Control=0
  y <- condition_data$y
  weights <- condition_data$weights
  
  # Simplified variables
  mu_11 <- mu$mu_11
  mu_10 <- mu$mu_10
  mu_01 <- mu$mu_01
  mu_00 <- mu$mu_00
  
  # DR Estimator for RC
  # Weights
  w_treat_post = weights * D * post
  w_treat_pre  = weights * D * (1 - post)
  w_cont_post  = weights * (1 - D) * post * pscore / (1 - pscore)
  w_cont_pre   = weights * (1 - D) * (1 - post) * pscore / (1 - pscore)
  
  mean_w_tr_post = mean(w_treat_post)
  mean_w_tr_pre  = mean(w_treat_pre)
  mean_w_co_post = mean(w_cont_post)
  mean_w_co_pre  = mean(w_cont_pre)
  
  if(mean_w_tr_post==0 || mean_w_tr_pre==0 || mean_w_co_post==0 || mean_w_co_pre==0) {
     warning("Weights sum to zero in one group/period.")
     return(list(dr_att=NA, inf_func=rep(0, nrow(data))))
  }

  # Influence Functions
  
  E_D = mean(weights * D)
  E_D_Post = mean(weights * D * post)
  E_D_Pre = mean(weights * D * (1-post))
  E_C_Post_w = mean(w_cont_post)
  E_C_Pre_w = mean(w_cont_pre)
  
  psi_11 = (weights * D * post * (y - mu_11)) / E_D_Post + (weights * D * mu_11) / E_D
  att_11 = mean(psi_11)
  
  psi_10 = (weights * D * (1-post) * (y - mu_10)) / E_D_Pre + (weights * D * mu_10) / E_D
  att_10 = mean(psi_10)
  
  psi_01 = (w_cont_post * (y - mu_01)) / E_C_Post_w + (weights * D * mu_01) / E_D 
  att_01 = mean(psi_01)
  
  psi_00 = (w_cont_pre * (y - mu_00)) / E_C_Pre_w + (weights * D * mu_00) / E_D 
  att_00 = mean(psi_00)
  
  dr_att = (att_11 - att_10) - (att_01 - att_00)
  
  inf_func = (psi_11 - att_11) - (psi_10 - att_10) - ((psi_01 - att_01) - (psi_00 - att_00))
  
  # fill zeros in influence function for observation outside the subgroup analyzed
  data[, "inf_func_result" := numeric(.N)]
  data[data$subgroup %in% c(condition_subgroup, 4), "inf_func_result" := inf_func]
  
  return(list(dr_att = dr_att, inf_func = data[["inf_func_result"]]))
}
