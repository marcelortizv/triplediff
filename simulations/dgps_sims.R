# ------------------------------------------------------------------------
# Purpose: This script contains the functions to generate the data for
#          Monte Carlo simulations
# Date: 05/29/2024
# ------------------------------------------------------------------------


fps <- function(psi_index, coefs, xvars){
  elem_wise_mult <- xvars %*% coefs
  return(psi_index * elem_wise_mult)
}

freg <- function(coefs, xvars){
  return(210 + xvars %*% coefs)
}

get_long_data <- function(y0, y1, state, partition, covs){
  # ensure that y0 and y1 have the same length
  if(length(y0) != length(y1)){
    stop("y0 and y1 must have the same length")
  }

  # check if covs has 4 columns
  if(ncol(covs) != 4){
    stop("covs must have 4 columns")
  }

  n <- length(y0)

  dt <- data.table(
    id = rep(1:n, 2),
    state = rep(state, 2),
    partition = rep(partition, 2),
    time = rep(1:2, each = n),
    y = c(y0, y1),
    cov1 = rep(covs[, 1], 2),
    cov2 = rep(covs[, 2], 2),
    cov3 = rep(covs[, 3], 2),
    cov4 = rep(covs[, 4], 2)
  )

  return(dt)
}

# Mean and Std deviation of Z's without truncation
mean.z1 <- exp(0.25/2)
sd.z1 <- sqrt((exp(0.25) - 1) * exp(0.25))

mean.z2 <- 10
sd.z2 <- 0.54164

mean.z3 <- 0.21887
sd.z3 <-   0.04453

mean.z4 <- 402
sd.z4 <-  56.63891

# Generate efficiency bounds for MC simulations
gen_dgp <- function(size, dgp_type){
  # Generate data for MC simulations
  # size: number of observations
  # dgp_type: type of DGP
  # Returns a data table in long format

  att <- 0
  # Gen covariates
  x1 <- stats::rnorm(size, mean = 0, sd = 1)
  x2 <- stats::rnorm(size, mean = 0, sd = 1)
  x3 <- stats::rnorm(size, mean = 0, sd = 1)
  x4 <- stats::rnorm(size, mean = 0, sd = 1)

  z1_tilde <- exp(x1/2)
  z2_tilde <- x2/(1 + exp(x1)) + 10
  z3_tilde <- (x1 * x3/25 + 0.6)^3
  z4_tilde <- (x1 + x4 + 20)^2


  z1 <- (z1_tilde - mean.z1) / sd.z1
  z2 <- (z2_tilde - mean.z2) / sd.z2
  z3 <- (z3_tilde - mean.z3) / sd.z3
  z4 <- (z4_tilde - mean.z4) / sd.z4

  x <- cbind(x1, x2, x3, x4)
  z <- cbind(z1, z2, z3, z4)

  w1 <- c(-1, .5, -.25, -.1)
  w2 <- c(-0.5, 2, .5, -.2)
  w3 <- c(3, -1.5, +.75, -.3)

  b1 <- c(27.4, 13.7, 13.7, 13.7)
  b2 <- c(6.85, 3.43, 3.43, 3.43)

  if(dgp_type == 1){
    # both models are functions of Z
    # index for pscores
    fps1 <- fps(0.3, w1, z)
    fps2 <- fps(0.2, w2, z)
    fps3 <- fps(0.05, w3, z)

    # index for outcome regression model
    freg1 <- freg(b1, z)
    freg0 <- freg(b2, z)
    eff <- 34.76

  } else if(dgp_type == 2){
    # Pscore depends on X but Regressions depend on Z
    # index for pscores
    fps1 <- fps(0.3, w1, x)
    fps2 <- fps(0.2, w2, x)
    fps3 <- fps(0.05, w3, x)

    # index for outcome regression model
    freg1 <- freg(b1, z)
    freg0 <- freg(b2, z)
    eff <- 34.72

  } else if(dgp_type == 3){
    # Pscore depends on X but Regressions depend on Z
    # index for pscores
    fps1 <- fps(0.3, w1, z)
    fps2 <- fps(0.2, w2, z)
    fps3 <- fps(0.05, w3, z)

    # index for outcome regression model
    freg1 <- freg(b1, x)
    freg0 <- freg(b2, x)
    eff <- 34.78
  } else if (dgp_type == 4){
    # Both models are functions of X (misspecified)
    # index for pscores
    fps1 <- fps(0.3, w1, x)
    fps2 <- fps(0.2, w2, x)
    fps3 <- fps(0.05, w3, x)

    # index for outcome regression model
    freg1 <- freg(b1, x)
    freg0 <- freg(b2, x)
    eff <- 34.73
  } else {
    stop("Invalid DGP type")
  }

  # Generate groups
  exp_f1 <- exp(fps1)
  exp_f2 <- exp(fps2)
  exp_f3 <- exp(fps3)

  sum_exp_f <- exp_f1 + exp_f2 + exp_f3

  p1 <- exp_f1 / (1 + sum_exp_f)
  p2 <- exp_f2 / (1 + sum_exp_f)
  p3 <- exp_f3 / (1 + sum_exp_f)
  p4 <- 1 - p1 - p2 - p3

  # determine treatment group PA based on U and cumulative probabilities

  U <- runif(size)

  PA <- sapply(1:size, function(i) {
    if (U[i] <= p1[i]) {
      return(1)
    } else if (U[i] <= p1[i] + p2[i]) {
      return(2)
    } else if (U[i] <= 1 - p4[i]) {
      return(3)
    } else {
      return(4)
    }
  })

  # Generate state variable
  state <- ifelse(PA == 3 | PA == 4, 1, 0)

  # Generate partition
  partition <- ifelse(PA == 2 | PA == 4, 1, 0)

  # Generate potential outcomes

  unobs_het <- state*partition*freg1 + (1 - state)*partition*freg0
  or_lin <- state*freg1 + (1 - state)*freg0
  # v is the unobserved heterogeneity
  v <- stats::rnorm(size, mean = unobs_het, sd = 1)

  #Gen realized outcome at time 0
  y0 <- or_lin + v + stats::rnorm(size, mean = 0, sd = 1)

  # gen outcomes at time 1
  # First let's generate potential outcomes: y_1_potential
  y10 <- or_lin + v + stats::rnorm(size, mean = 0, sd = 1) + # This is the baseline
    or_lin # this is for the trend based on X

  y11 <- or_lin + v + stats::rnorm(size, mean = 0, sd = 1) +#This is the baseline
    or_lin + #this is for the trend based on X
    att # This is the treatment effects

  # Get unfeasible att
  att.unf <- (mean(state*partition*y11) - mean(state*partition*y10))/mean(state*partition)

  # Now we generate the realized outcomes
  y1 <- state*partition*y11 + (1 - state*partition)*y10

  # Generate the data assuming we always observe Z
  dt <- get_long_data(y0, y1, state, partition, z)

  return(list(data = dt, att = 0, att.unf = att.unf, eff = eff))

}
