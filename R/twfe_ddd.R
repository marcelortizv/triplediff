# Complementary function to run TWFE DDD estimator
NULL
#' TWFE DDD estimator, with panel data and 2 periods
#'
#' @description
#' \code{twfe_ddd} function implements a two-way fixed effect regression for assessing the average
#' treatment effect on the treated (ATT) using a triple differences (DDD) specification
#' in panel data settings across two time periods.
#'
#' @param yname A character string specifying the name of the outcome variable in the data.
#' @param tname A character string specifying the name of the time variable in the data.
#' @param dname A character string specifying the name of the state variable in the data.
#' @param pname A character string specifying the name of the partition variable in the data.
#' @param xformla A formula specifying the covariates to be included in the regression.
#' @param data A data frame containing the preprocessed data.
#'
#' @return A list with the estimated ATT, standard error, upper and lower confidence intervals, and influence function.
# #' @noRd
#' @import stats
#' @export

twfe_ddd <- function(yname, tname, dname,
                    pname, xformla = ~1, data) {

  # Extract data
  y <- data[[yname]]
  state <- data[[dname]]
  partition <- data[[pname]]

  # Creating a post dummy variable based on tlist[2] (second period = post treatment)
  tlist <- unique(data[[tname]])[base::order(unique(data[[tname]]))]
  post <- as.numeric(data[[tname]] == tlist[2])

  # get matrix of covariates based on xformla
  X <- model.matrix(xformla, data)

  # Estimate TWFE regression
  twfe <- stats::lm(y ~ X + state + partition + post + state:partition + state:post + partition:post + state:partition:post - 1, x = TRUE)

  twfe_att <- twfe$coefficients["state:partition:post"][[1]]

  # compute influence function
  inf.twfe <- (twfe$x * twfe$residuals) %*%
    base::qr.solve(crossprod(twfe$x, twfe$x) / dim(twfe$x)[1])

  sel.theta <- matrix(c(rep(0, dim(inf.twfe)[2])))

  index.theta <- which(dimnames(twfe$x)[[2]]=="state:partition:post",
                       arr.ind = TRUE)

  sel.theta[index.theta, ] <- 1

  #get the influence function of the TWFE regression
  twfe_inf_func <- as.vector(inf.twfe %*% sel.theta)

  # Compute standard errors
  # Estimate of standard error
  se_twfe <- stats::sd(twfe_inf_func)/sqrt(length(twfe_inf_func))
  # Estimate of upper boudary of 95% CI
  ci_upper <- twfe_att + 1.96 * se_twfe
  # Estimate of lower doundary of 95% CI
  ci_lower <- twfe_att - 1.96 * se_twfe



  return(list(ATT = twfe_att,
              se = se_twfe,
              uci = ci_upper,
              lci = ci_lower,
              inf_func = twfe_inf_func))
}
