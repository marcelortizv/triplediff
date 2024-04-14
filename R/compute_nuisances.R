#' Function to compute nuisances parameters and did
#' @import stats
#' @import parglm
#' @import speedglm
#' @import data.table
#' @noRd
#--------------------------------------------------

# utility function to generate equation to estimate pscore
get_formula_pscore <- function(xformula){
  # Get a formula object for pscore estimation
  formula_obj <- stats::update(xformula, PA4 ~ .)

  return(formula_obj)
}


# Function to compute propensity scores using parglm for multiple subgroups
compute_pscore <- function(data, condition_subgroup, xformula) {
  # get formula for pscore estimation using covariates
  formula_pscore <- get_formula_pscore(xformula)
  # Subset data for condition_subgroup and subgroup == 4 or the given condition_subgroup
  condition_data <- data[data$subgroup %in% c(condition_subgroup, 4)]
  # Adding treatment variable P(1{PA4 = 4}|X)
  condition_data[, "PA4" := ifelse(condition_data$subgroup == 4, 1, 0)]
  uid_condition_data <- unique(condition_data, by = "id")

  # Fit logistic regression model using parglm
  model <- parglm::parglm(formula_pscore, data = uid_condition_data,
                          family = stats::binomial(),
                          weights = weights,
                          control = parglm.control(nthreads = getDTthreads()),
                          intercept = FALSE)

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
  if (any(propensity_scores < 0.01) | any(propensity_scores > 0.99)) {
    warning(paste("Propensity scores for subgroup", condition_subgroup,
                  "have poor overlap."))
  }

  # Avoid divide by zero
  propensity_scores <- pmin(propensity_scores, 1 - 1e-16)

  # Save Hessian matrix for influence function
  hessian_matrix <- stats::vcov(model) * nrow(uid_condition_data)

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
  i.weights = control_data[control_data$post == 0, weights]

  # get coefficients of outcome regression model for the control group (subgroup != 4)
  # Attempt to compute reg.coeff using speedglm
  try_speedglm <- tryCatch({
    reg.coeff <- stats::coef(speedglm::speedglm.wfit(y = deltaY_control,
                                                     X = as.matrix(cov_control),
                                                     intercept = FALSE,
                                                     weights = i.weights))
  }, error = function(e) {
    # If error, attempt to compute reg.coeff using lm.wfit
    reg.coeff <- tryCatch({
      stats::coef(stats::lm.wfit(x = as.matrix(cov_control),
                                 y = deltaY_control,
                                 w = i.weights))
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
  cov <- stats::model.matrix(xformula, data = condition_data[condition_data$post == 0])
  # compute OR delta
  or_delta = as.vector(tcrossprod(reg.coeff, as.matrix(cov)))
  # compute regression adjustment
  # reg_adjust = deltaY - or_delta
  return(list(deltaY = deltaY, or_delta = or_delta))
}

# Function to compute the average treatment effect for multiple subgroups
compute_did <- function(data, condition_subgroup, pscores, reg_adjustment, xformula){

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
  i.weights = condition_data[condition_data$post == 0, weights]

  ################################
  # Get doubly-robust estimation #
  ################################

  w.treat = i.weights * PA4
  w.control = (i.weights * pscore * PAa) / (1 - pscore)
  riesz.treat = w.treat * (deltaY - or_delta)
  riesz.control = w.control * (deltaY - or_delta)
  att.treat = mean(riesz.treat, na.rm = TRUE) / mean(w.treat, na.rm = TRUE)
  att.control = mean(riesz.control, na.rm = TRUE) / mean(w.control, na.rm = TRUE)

  dr.att = att.treat - att.control

  ##########################
  # Get influence function #
  ##########################

  # Influence function related to the estimation of pscores
  cov = stats::model.matrix(xformula, data = condition_data)
  M2 <- base::colMeans(w.control * (deltaY - or_delta - att.control) * cov, na.rm = TRUE) # reg_adjust = deltaY - m_delta(x)

  score_ps <- i.weights * (PA4 - pscore) * cov
  # asymptotic linear representation of logit's beta
  score_ps_no_na <- na.omit(score_ps) # Exclude rows with NA values
  asy.lin.rep.ps <- score_ps_no_na %*% hessian
  inf.control.pscore <- asy.lin.rep.ps %*% as.matrix(M2)

  # Influence function related to the estimation of pscores

  M1 <- base::colMeans(w.treat * cov, na.rm = TRUE)
  M3 <- base::colMeans(w.control * cov, na.rm = TRUE)

  # Influence function related to the estimation of regression model
  or_x <- i.weights * PAa * cov
  or_ex <- i.weights * PAa * (deltaY - or_delta) * cov
  XpX <- crossprod(or_x, cov)/nrow(condition_data)

  #asymptotic linear representation of the beta
  asy.linear.or <- t(solve(XpX, t(or_ex)))

  #or for treat
  inf.treat.or <- -asy.linear.or %*% M1
  #or for control
  inf.cont.or <- -asy.linear.or %*% M3

  # Influence function from did
  inf.control.did <- riesz.control - w.control*att.control
  inf.treat.did <- riesz.treat - w.treat*att.treat

  # Influence function for the ATT
  inf.control <- (inf.control.did + inf.control.pscore + inf.cont.or)/mean(w.control)
  inf.treat <- (inf.treat.did + inf.treat.or)/mean(w.treat)
  # putting all together
  inf.func <- inf.treat - inf.control

  # fill zeros in influence function for observation outside the subgroup analyzed
  data[, "inf.func.result" := numeric(.N)]
  data[data$subgroup %in% c(condition_subgroup, 4), "inf.func.result" := inf.func]

  return(list(dr.att = dr.att, inf.func = data[["inf.func.result"]]))
}
