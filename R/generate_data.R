#' Function to generate a fake dataset for testing purposes only.
#' @description
#' Function to generate fake dataset to test internal procedures.
#'
#' @param seed Seed for reproducibility
#' @param num_ids Number of IDs
#' @param time Number of time periods
#' @param initial.year Initial year
#' @param treatment.year Treatment year
#'
#' @return A data.table with the following columns:
#'
#' - id: ID
#' - state: State variable
#' - year: Time variable
#' - partition: Partition variable
#' - x1: Covariate 1
#' - x2: Covariate 2
#' - treat: Treatment variable
#' - outcome: Outcome variable
#'
#' @keywords internal
generate_test_panel <- function(seed = 123,
                                num_ids = 100,
                                time = 2,
                                initial.year = 2019,
                                treatment.year = 2020) {
  # Flag for initial year and treatment year
  if (initial.year >= treatment.year) {
    stop("The initial year must be less than the treatment year")
  }

  if (time < 2) {
    stop("The time variable must be greater than 1")
  }

  if (num_ids < 1) {
    stop("The number of IDs must be greater than 0")
  }

  if ((treatment.year - initial.year +1 ) > time){
    stop("The difference between the treatment year and the initial year must be less or equal than the time variable")
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Number of IDs
  num_ids <- num_ids

  # Number of observations per ID
  num_obs <- time

  # Generate IDs
  id = rep(1:num_ids, each = num_obs)

  # Generate x1 and x2 for each state
  x1 <- rep(rnorm(num_ids), each = num_obs)
  x2 <- rep(rnorm(num_ids), each = num_obs)

  # Generate state variable
  state = rep(sample(c("A", "B", "C", "D"), num_ids, replace = TRUE), each = num_obs)

  # generate time variable
  year = rep(seq(initial.year, initial.year+(time-1), 1), times = num_ids)

  # Generate partition variable
  partition = rep(rbinom(num_ids, 1, 0.5), each = num_obs)

  # Generate target group
  s_2 = as.numeric(state %in% c("A", "B"))

  # Generate post treatment variable
  t_2 = as.numeric(year >= treatment.year)

  # TWFE outcome variable
  y = 2 + 5*s_2 - 2 * partition + 3 * t_2 +
    4*s_2*partition + 2*s_2*t_2 + 3*partition*t_2 +
    1*s_2*partition*t_2 + rnorm(num_ids*num_obs)

  # Create a data.table
  dt <- data.table::data.table(
    id = id,
    state = state,
    year = year,
    partition = partition,
    x1 = x1,
    x2 = x2,
    treat = s_2,
    #post = t_2,
    outcome = y
  )
  return(dt)
}

# ---------------------------------------------------------------------
# Functions to generate data for 2 time periods
# ---------------------------------------------------------------------

fps <- function(psi_index, coefs, xvars){
  elem_wise_mult <- xvars %*% coefs
  return(psi_index * elem_wise_mult)
}

fps2 <- function(psi_index, coefs, xvars, c){
  elem_wise_mult <- c + xvars %*% coefs
  return(psi_index * elem_wise_mult)
}

freg <- function(coefs, xvars){
  return(210 + xvars %*% coefs)
}

get_long_data <- function(y_list, state, partition, covs){
  # ensure that y_list is a list of vectors
  if (!is.list(y_list) || any(!sapply(y_list, is.vector))) {
    stop("y_list must be a list of vectors")
  }

  # ensure that all y vectors have the same length
  y_length <- sapply(y_list, length)
  if (length(unique(y_length)) != 1) {
    stop("All y vectors must have the same length")
  }

  # check if covs has 4 columns
  if(ncol(covs) != 4){
    stop("covs must have 4 columns")
  }

  n <- y_length[1]
  t <- length(y_list)

  dt <- data.table(
    id = rep(1:n, t),
    state = rep(state, t),
    partition = rep(partition, t),
    time = rep(1:t, each = n),
    y = unlist(y_list),
    cov1 = rep(covs[, 1], t),
    cov2 = rep(covs[, 2], t),
    cov3 = rep(covs[, 3], t),
    cov4 = rep(covs[, 4], t)
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

# ---------------------------------------------------------------------
# Generate data for 2 time treatment periods
# ---------------------------------------------------------------------

#' Function that generates panel data with single treatment date assignment and two time periods.
#' @description
#' Generate panel data with a single treatment date and two periods
#'
#' @param size Integer. Number of units.
#' @param dgp_type Integer in \{1,2,3,4\}.
#'   1 = both nuisance functions correct;
#'   2 = only the outcome model correct;
#'   3 = only the propensity score correct;
#'   4 = both nuisance functions incorrect.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{data}{A \code{data.table} in long format with columns:
#'     \itemize{
#'       \item \code{id}: unit identifier
#'       \item \code{state}: state variable
#'       \item \code{time}: time variable
#'       \item \code{partition}: partition assignment
#'       \item \code{x1}, \code{x2}, \code{x3}, \code{x4}: covariates
#'       \item \code{y}: outcome variable
#'       \item \code{cluster}: cluster ID (no within-cluster correlation)
#'     }
#'   }
#'   \item{att}{True average treatment effect on the treated (ATT), set to 0.}
#'   \item{att.unf}{Oracle ATT computed under the unfeasible specification.}
#'   \item{eff}{Theoretical efficiency bound for the estimator.}
#' }
#'
#' @export
gen_dgp_2periods <- function(size, dgp_type){
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
    fps1 <- fps(0.2, w1, z)
    fps2 <- fps(0.2, w2, z)
    fps3 <- fps(0.05, w3, z)

    # index for outcome regression model
    freg1 <- freg(b1, z)
    freg0 <- freg(b2, z)
    eff <- 32.82

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
    dt <- get_long_data(list(as.vector(y0), as.vector(y1)), state, partition, z)
    setorder(dt, "id", "time")
    # generate clusters (there's no actual within-cluster correlation) with sample with replacement by id
    # get unique ids
    unique_ids <- unique(dt$id)
    # sample cluster values for each unique id
    clusters <- data.table(id = unique_ids, cluster = sample(1:50, size = length(unique_ids), replace = TRUE))
    # merge the sampled clusters back into the original data table
    dt <- merge(dt, clusters, by = "id", all.x = TRUE)

    return(list(data = dt, att = 0, att.unf = att.unf, eff = eff))


  } else if(dgp_type == 2){
    # Pscore depends on X but Regressions depend on Z
    # index for pscores
    fps1 <- fps(0.2, w1, x)
    fps2 <- fps(0.2, w2, x)
    fps3 <- fps(0.05, w3, x)

    # index for outcome regression model
    freg1 <- freg(b1, z)
    freg0 <- freg(b2, z)
    eff <- 32.52

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
    dt <- get_long_data(list(as.vector(y0), as.vector(y1)), state, partition, z)
    setorder(dt, "id", "time")
    # generate clusters (there's no actual within-cluster correlation) with sample with replacement by id
    # get unique ids
    unique_ids <- unique(dt$id)
    # sample cluster values for each unique id
    clusters <- data.table(id = unique_ids, cluster = sample(1:50, size = length(unique_ids), replace = TRUE))
    # merge the sampled clusters back into the original data table
    dt <- merge(dt, clusters, by = "id", all.x = TRUE)

    return(list(data = dt, att = 0, att.unf = att.unf, eff = eff))


  } else if(dgp_type == 3){
    # Pscore depends on Z but Regressions depend on X
    # index for pscores
    fps1 <- fps(0.2, w1, z)
    fps2 <- fps(0.2, w2, z)
    fps3 <- fps(0.05, w3, z)

    # index for outcome regression model
    freg1 <- freg(b1, x)
    freg0 <- freg(b2, x)
    eff <- 32.82

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
    dt <- get_long_data(list(as.vector(y0), as.vector(y1)), state, partition, z)
    setorder(dt, "id", "time")
    # generate clusters (there's no actual within-cluster correlation) with sample with replacement by id
    # get unique ids
    unique_ids <- unique(dt$id)
    # sample cluster values for each unique id
    clusters <- data.table(id = unique_ids, cluster = sample(1:50, size = length(unique_ids), replace = TRUE))
    # merge the sampled clusters back into the original data table
    dt <- merge(dt, clusters, by = "id", all.x = TRUE)

    return(list(data = dt, att = 0, att.unf = att.unf, eff = eff))

  } else if (dgp_type == 4){
    # Both models are functions of X (misspecified)
    # index for pscores
    fps1 <- fps(0.2, w1, x)
    fps2 <- fps(0.2, w2, x)
    fps3 <- fps(0.05, w3, x)

    # index for outcome regression model
    freg1 <- freg(b1, x)
    freg0 <- freg(b2, x)
    eff <- 32.52

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
    dt <- get_long_data(list(as.vector(y0), as.vector(y1)), state, partition, z)
    setorder(dt, "id", "time")
    # generate clusters (there's no actual within-cluster correlation) with sample with replacement by id
    # get unique ids
    unique_ids <- unique(dt$id)
    # sample cluster values for each unique id
    clusters <- data.table(id = unique_ids, cluster = sample(1:50, size = length(unique_ids), replace = TRUE))
    # merge the sampled clusters back into the original data table
    dt <- merge(dt, clusters, by = "id", all.x = TRUE)

    return(list(data = dt, att = 0, att.unf = att.unf, eff = eff))

  } else {
    stop("Invalid DGP type")
  }

}

# ---------------------------------------------------------------------
# Functions to generate data for multiple time periods
# ---------------------------------------------------------------------

#' Generate panel data with staggered treatment adoption (three periods)
#' @description
#' Generate panel data where units adopt treatment at different times across three periods.
#'
#' @param size Integer. Number of units to simulate.
#' @param dgp_type Integer in \{1,2,3,4\}. Only used when \code{include_covariates = TRUE}.
#'   1 = both nuisance functions correct (default);
#'   2 = only the outcome model correct;
#'   3 = only the propensity-score model correct;
#'   4 = both nuisance functions misspecified.
#' @param include_covariates Logical. If \code{TRUE} (default), generates covariates with
#'   transformations and uses \code{dgp_type} specification. If \code{FALSE}, uses constant
#'   covariates and fixed propensity score probabilities for a simpler DGP.
#'
#' @return A named list with components:
#' \describe{
#'   \item{data}{A \code{data.table} in long format with columns:
#'     \itemize{
#'       \item \code{id}: unit identifier
#'       \item \code{cohort}: first period when treatment is assigned
#'       \item \code{partition}: partition indicator
#'       \item \code{x1}, \code{x2}, \code{x3}, \code{x4}: covariates
#'       \item \code{cluster}: cluster identifier (no within-cluster correlation)
#'       \item \code{time}: time period index
#'       \item \code{y}: observed outcome
#'     }
#'   }
#'   \item{data_wide}{A \code{data.table} in wide format (one row per \code{id}) with columns:
#'     \itemize{
#'       \item \code{id}, \code{cohort}, \code{partition}, \code{x1}, \code{x2}, \code{x3}, \code{x4}, \code{cluster}
#'       \item \code{y_t0}, \code{y_t1}, \code{y_t2}: outcomes in periods 0, 1, and 2
#'     }
#'   }
#'   \item{ES_0_unf}{(Only if \code{include_covariates = TRUE}) Unfeasible (oracle) event-study parameter at time 0.}
#'   \item{prob_g2_p1}{(Only if \code{include_covariates = TRUE}) Proportion of units with \code{cohort == 2} and eligibility in period 1.}
#'   \item{prob_g3_p1}{(Only if \code{include_covariates = TRUE}) Proportion of units with \code{cohort == 3} and eligibility in period 1.}
#' }
#'
#' @export
gen_dgp_mult_periods <- function(size, dgp_type = 1, include_covariates = TRUE){

  # Validate parameters
  if (!include_covariates && dgp_type != 1) {
    warning("dgp_type is ignored when include_covariates = FALSE. Using simple DGP specification.")
  }

  # Generates covariates
  Xsi.ps = 0.4 # This scale things in pscore

  if (include_covariates) {
    x1 <- stats::rnorm(size, mean = 0, sd = 1)
    x2 <- stats::rnorm(size, mean = 0, sd = 1)
    x3 <- stats::rnorm(size, mean = 0, sd = 1)
    x4 <- stats::rnorm(size, mean = 0, sd = 1)
  } else {
    x1 <- rep(1, size)
    x2 <- rep(1, size)
    x3 <- rep(1, size)
    x4 <- rep(1, size)
  }

  x <- cbind(x1, x2, x3, x4)

  # Transform to Z only if covariates are included
  if (include_covariates) {
    z1_tilde <- exp(x1/2)
    z2_tilde <- x2/(1 + exp(x1)) + 10
    z3_tilde <- (x1 * x3/25 + 0.6)^3
    z4_tilde <- (x1 + x4 + 20)^21

    z1 <- (z1_tilde - mean.z1) / sd.z1
    z2 <- (z2_tilde - mean.z2) / sd.z2
    z3 <- (z3_tilde - mean.z3) / sd.z3
    z4 <- (z4_tilde - mean.z4) / sd.z4

    z <- cbind(z1, z2, z3, z4)
  } else {
    # For no covariates, z = x (constants)
    z <- x
  }

  w1 <- c(-1, 0.5, -0.25, -0.1) # for (2,A), (2,B)
  w2 <- c(-0.5, 1, -0.1, -0.25) # for (3,A), (3,B)
  w3 <- c(-0.25, 0.1, -1, -0.1) # for (0,A), (0,B)

  b1 <- c(27.4, 13.7, 13.7, 13.7)

  #-----------------------------------------------------------------------------
  # Generate groups and partitions
  # DGP has 3 timing groups, G=0,2,3, three timing periods, t=1,2,3,
  # and two partitions for each group, P=0,1 (P=0 is not eligible for treatment)
  # Thus, in total, we have 6 Group-Partition combinations

  if (!include_covariates) {
    #---------------------------------------------------------------------------
    # SIMPLE DGP: No covariates (constants only)
    #---------------------------------------------------------------------------
    # Use fixed probability matrix for propensity scores
    probs_pscore <- matrix(cbind(0.20, 0.15, 0.30, 0.20, 0.05, 0.10),
                           nrow=size, ncol=6, byrow = TRUE)
    colnames(probs_pscore) <- c("pi_2A", "pi_2B", "pi_3A", "pi_3B", "pi_0A", "pi_0B")

    # Sample group types
    group_types <- apply(probs_pscore, 1,
                         function(pvec) sample(x = seq(1,6), size=1, prob=pvec))

    # Translate group types into partitions and cohorts
    partition <- ifelse(group_types %in% c(1,3,5), 1, 0)
    cohort <- ifelse(group_types %in% c(1,2),
                     2,
                     ifelse(group_types %in% c(3,4),
                            3,
                            0))

    #---------------------------------------------------------------------------
    # Create indexes for potential outcomes (simple specification)
    #---------------------------------------------------------------------------
    # Index for PO (using simpler linear index without Z transformation)
    index_lin <- 210 + 27.4*x1 + 13.7*(x2 + x3 + x4)
    index_partition <- partition * index_lin
    # Index for unobserved heterogeneity
    index_unobs_het <- cohort * index_lin + index_partition
    # Index for Conditional Parallel Trends
    index_trend <- index_lin
    # Generate unobserved heterogeneity
    v <- stats::rnorm(size, mean = index_unobs_het, sd = 1)
    # Index to violates the parallel trends assumption
    index_pt_violation <- v/10

    # Index for per-period ATT(g,t)
    index_att_g2 <- 10
    index_att_g3 <- 25

    #---------------------------------------------------------------------------
    # Generate the potential outcomes
    #---------------------------------------------------------------------------
    # Baseline index
    baseline_t1 <- index_lin + index_partition + v

    # Generate untreated outcome at time 1
    y_t1 <- baseline_t1 + stats::rnorm(size)

    # Generate potential outcomes at time 2
    baseline_t2 <- baseline_t1 +
      index_pt_violation + # Violate the CPT
      index_trend # trend based on X

    y_t2_never <- baseline_t2 + stats::rnorm(size)
    y_t2_g2 <- baseline_t2 + stats::rnorm(size) +
      index_att_g2 * partition

    # Generate potential outcomes at time 3
    baseline_t3 <- baseline_t1 +
      2*index_trend +
      2*index_pt_violation

    y_t3_never <- baseline_t3 + stats::rnorm(size)
    y_t3_g2 <- baseline_t3 + stats::rnorm(size) +
      2*index_att_g2 * partition
    y_t3_g3 <- baseline_t3 + stats::rnorm(size) +
      index_att_g3 * partition

    #---------------------------------------------------------------------------
    # Get the realized outcomes
    #---------------------------------------------------------------------------
    y_t2 <- ifelse(cohort == 2 & partition==1, y_t2_g2, y_t2_never)
    y_t3 <- ifelse(cohort == 2 & partition==1, y_t3_g2,
                   ifelse(cohort == 3 & partition==1, y_t3_g3,
                          y_t3_never))

    #---------------------------------------------------------------------------
    # Put all data in a data frame
    #---------------------------------------------------------------------------
    dta_wide <- data.frame(cbind(id = 1:size,
                                 y_t1, y_t2 , y_t3,
                                 cohort , partition,
                                 x))
    # Transform data into long format
    dta <- get_long_data(list(as.vector(y_t1), as.vector(y_t2), as.vector(y_t3)), cohort, partition, x)
    setorder(dta, "id", "time")
    # generate clusters
    unique_ids <- unique(dta$id)
    clusters <- data.table(id = unique_ids, cluster = sample(1:50, size = length(unique_ids), replace = TRUE))
    dta <- merge(dta, clusters, by = "id", all.x = TRUE)

    #---------------------------------------------------------------------------
    # Return the data (no unfeasible parameters for simple DGP)
    #---------------------------------------------------------------------------
    return(list(data = as.data.table(dta),
                data_wide = as.data.table(dta_wide)))

  } else if(dgp_type == 1){
    # both models are functions of Z
    # index for pscores

    # Generalize Propensity score (3 groups and two partitions)
    pi_2A <- exp(fps2(Xsi.ps, w1, z, 1.25))
    pi_2B <- exp(fps2(-Xsi.ps, w1, z, -0.5))
    pi_3A <- exp(fps2(Xsi.ps, w2, z, 2))
    pi_3B <- exp(fps2(-Xsi.ps, w2, z, -1.25))
    pi_0A <- exp(fps2(Xsi.ps, w3, z, -0.5))
    # Get the sum of the exponentials
    sum_pi <- 1 + pi_2A + pi_2B + pi_3A + pi_3B + pi_0A
    # Construct the generalized pscore
    pi_2A <- pi_2A/sum_pi
    pi_2B <- pi_2B/sum_pi
    pi_3A <- pi_3A/sum_pi
    pi_3B <- pi_3B/sum_pi
    pi_0A <- pi_0A/sum_pi
    pi_0B <- 1 - (pi_2A + pi_2B + pi_3A + pi_3B + pi_0A)

    # Bind them into a matrix
    probs_pscore <- cbind(pi_2A, pi_2B, pi_3A, pi_3B, pi_0A, pi_0B)

    # Sample group types (1 and 3 are the group of interest)
    group_types <- apply(probs_pscore, 1,
                         function(pvec) sample(seq(1,6), size=1, prob=pvec))

    # Translate group types into partitions and cohorts
    partition <- ifelse(group_types %in% c(1,3,5), 1, 0)
    cohort <- ifelse(group_types %in% c(1,2),
                     2,
                     ifelse(group_types %in% c(3,4),
                            3,
                            0))
    #-----------------------------------------------------------------------------
    # Now create the indexes for the potential outcomes
    #-----------------------------------------------------------------------------
    # Index for PO
    index_lin <- freg(b1, z)
    index_partition <- partition * index_lin
    # Index for unobserved heterogeneity
    index_unobs_het <- cohort * index_lin + index_partition
    #Index for Conditional Parallel Trends (this does not vary with treatment group)
    index_trend <- index_lin
    # Generate unobserved heterogeneity
    # (fixed effects that are correlated with X, partition, and cohort)
    v <- stats::rnorm(size, mean = index_unobs_het, sd = 1)
    # Index to violates the parallel trends assumption
    index_pt_violation <- v/10

    # Index for per-period ATT(g,t)
    index_att_g2 <- 10
    index_att_g3 <- 25
    #-----------------------------------------------------------------------------
    # Generate the potential outcomes
    # Baseline index
    baseline_t1 <- index_lin + index_partition + v

    #Generate untreated outcome at time 1 (that is observed for all units)
    y_t1 <- baseline_t1 + stats::rnorm(size)

    # Generate potential outcomes at time 2
    baseline_t2 <- baseline_t1 +
      index_pt_violation + #Violate the CPT
      index_trend #this is for the trend based on X

    y_t2_never <- baseline_t2 + stats::rnorm(size)
    y_t2_g2 <- baseline_t2 + stats::rnorm(size) +
      index_att_g2 * partition # This is the treatment effects for cohort 2
    y_t2_g3 <-  baseline_t2 + stats::rnorm(size)

    # Generate potential outcomes at time 3
    baseline_t3 <- baseline_t1 +
      2*index_trend + #this is for the trend based on X
      2*index_pt_violation #Violate the CPT

    y_t3_never <- baseline_t3 + stats::rnorm(size)
    y_t3_g2 <- baseline_t3 + stats::rnorm(size) +
      2*index_att_g2 * partition # This is the treatment effects for cohort 2
    y_t3_g3 <- baseline_t3 + stats::rnorm(size) +
      index_att_g3 * partition # This is the treatment effects for cohort 2
    #-----------------------------------------------------------------------------
    # Get the realized outcomes now
    # y_t2<- y_t2_g2*(cohort == 2) + y_t2_g3*(cohort==3) + y_t2_never*(cohort==0)
    # y_t3<- y_t3_g2*(cohort == 2) + y_t3_g3*(cohort==3) + y_t3_never*(cohort==0)

    y_t2<-  ifelse(cohort == 2 & partition==1, y_t2_g2, y_t2_never)
    y_t3<- ifelse(cohort == 2 & partition==1, y_t3_g2,
                  ifelse(cohort == 3 & partition==1, y_t3_g3,
                         y_t3_never))
    #-----------------------------------------------------------------------------
    # Get the event study parameter
    # Get unfeasible ATT(g,t) (assuming we see the Potential outcomes)
    att_g2_t2_unf <- (mean((group_types==1)*y_t2_g2) - mean((group_types==1)*y_t2_never))/
      mean((group_types==1))
    att_g2_t3_unf <- (mean((group_types==1)*y_t3_g2) - mean((group_types==1)*y_t3_never))/
      mean((group_types==1))
    att_g3_t2_unf <- (mean((group_types==3)*y_t2_g3) - mean((group_types==3)*y_t2_never))/
      mean((group_types==3))
    att_g3_t3_unf <- (mean((group_types==3)*y_t3_g3) - mean((group_types==3)*y_t3_never))/
      mean((group_types==3))

    # Get unfeasible Event-Study (ES) (assuming we see the Potential outcomes and pscores)
    # get weights
    prob_g2_p1 <- mean(pi_2A/(pi_2A + pi_3A))
    prob_g3_p1 <- mean(pi_3A/(pi_2A + pi_3A))
    # Get ES
    ES_0_unf <- att_g2_t2_unf * prob_g2_p1 + att_g3_t3_unf * prob_g3_p1
    ES_1_unf <- att_g2_t3_unf
    # ------------------------------------------------------------------------------
    # Put all data in a data frame
    dta_wide <- data.frame(cbind(id = 1:size,
                                 y_t1 = y_t1,
                                 y_t2 = y_t2,
                                 y_t3 = y_t3,
                                 cohort = cohort,
                                 partition = partition,
                                 z = z)) # we always observe z
    # Transform data into long format
    # Generate the data assuming we always observe Z
    dta <- get_long_data(list(as.vector(y_t1), as.vector(y_t2), as.vector(y_t3)), cohort, partition, z)
    setorder(dta, "id", "time")
    # generate clusters (there's no actual within-cluster correlation) with sample with replacement by id
    # get unique ids
    unique_ids <- unique(dta$id)
    # sample cluster values for each unique id
    clusters <- data.table(id = unique_ids, cluster = sample(1:50, size = length(unique_ids), replace = TRUE))
    # merge the sampled clusters back into the original data table
    dta <- merge(dta, clusters, by = "id", all.x = TRUE)
    #-----------------------------------------------------------------------------
    # Return the data
    return(list(data = as.data.table(dta),
                data_wide = as.data.table(dta_wide),
                ES_0_unf = ES_0_unf,
                prob_g2_p1 = prob_g2_p1,
                prob_g3_p1 = prob_g3_p1))

  } else if(dgp_type == 2){
    # Pscore depends on X but Regressions depend on Z
    # index for pscores

    # Generalize Propensity score (3 groups and two partitions)
    pi_2A <- exp(fps2(Xsi.ps, w1, x, 1.25))
    pi_2B <- exp(fps2(-Xsi.ps, w1, x, -0.5))
    pi_3A <- exp(fps2(Xsi.ps, w2, x, 2))
    pi_3B <- exp(fps2(-Xsi.ps, w2, x, -1.25))
    pi_0A <- exp(fps2(Xsi.ps, w3, x, -0.5))
    # Get the sum of the exponentials
    sum_pi <- 1 + pi_2A + pi_2B + pi_3A + pi_3B + pi_0A
    # Construct the generalized pscore
    pi_2A <- pi_2A/sum_pi
    pi_2B <- pi_2B/sum_pi
    pi_3A <- pi_3A/sum_pi
    pi_3B <- pi_3B/sum_pi
    pi_0A <- pi_0A/sum_pi
    pi_0B <- 1 - (pi_2A + pi_2B + pi_3A + pi_3B + pi_0A)

    # Bind them into a matrix
    probs_pscore <- cbind(pi_2A, pi_2B, pi_3A, pi_3B, pi_0A, pi_0B)

    # Sample group types (1 and 3 are the group of interest)
    group_types <- apply(probs_pscore, 1,
                         function(pvec) sample(seq(1,6), size=1, prob=pvec))

    # Translate group types into partitions and cohorts
    partition <- ifelse(group_types %in% c(1,3,5), 1, 0)
    cohort <- ifelse(group_types %in% c(1,2),
                     2,
                     ifelse(group_types %in% c(3,4),
                            3,
                            0))
    #-----------------------------------------------------------------------------
    # Now create the indexes for the potential outcomes
    #-----------------------------------------------------------------------------
    # Index for PO
    index_lin <- freg(b1, z)
    index_partition <- partition * index_lin
    # Index for unobserved heterogeneity
    index_unobs_het <- cohort * index_lin + index_partition
    #Index for Conditional Parallel Trends (this does not vary with treatment group)
    index_trend <- index_lin
    # Generate unobserved heterogeneity
    # (fixed effects that are correlated with X, partition, and cohort)
    v <- stats::rnorm(size, mean = index_unobs_het, sd = 1)
    # Index to violates the parallel trends assumption
    index_pt_violation <- v/10

    # Index for per-period ATT(g,t)
    index_att_g2 <- 10
    index_att_g3 <- 25
    #-----------------------------------------------------------------------------
    # Generate the potential outcomes
    # Baseline index
    baseline_t1 <- index_lin + index_partition + v

    #Generate untreated outcome at time 1 (that is observed for all units)
    y_t1 <- baseline_t1 + stats::rnorm(size)

    # Generate potential outcomes at time 2
    baseline_t2 <- baseline_t1 +
      index_pt_violation + #Violate the CPT
      index_trend #this is for the trend based on X

    y_t2_never <- baseline_t2 + stats::rnorm(size)
    y_t2_g2 <- baseline_t2 + stats::rnorm(size) +
      index_att_g2 * partition # This is the treatment effects for cohort 2
    y_t2_g3 <-  baseline_t2 + stats::rnorm(size)

    # Generate potential outcomes at time 3
    baseline_t3 <- baseline_t1 +
      2*index_trend + #this is for the trend based on X
      2*index_pt_violation #Violate the CPT

    y_t3_never <- baseline_t3 + stats::rnorm(size)
    y_t3_g2 <- baseline_t3 + stats::rnorm(size) +
      2*index_att_g2 * partition # This is the treatment effects for cohort 2
    y_t3_g3 <- baseline_t3 + stats::rnorm(size) +
      index_att_g3 * partition # This is the treatment effects for cohort 2
    #-----------------------------------------------------------------------------
    # Get the realized outcomes now
    # y_t2<- y_t2_g2*(cohort == 2) + y_t2_g3*(cohort==3) + y_t2_never*(cohort==0)
    # y_t3<- y_t3_g2*(cohort == 2) + y_t3_g3*(cohort==3) + y_t3_never*(cohort==0)

    y_t2<-  ifelse(cohort == 2 & partition==1, y_t2_g2, y_t2_never)
    y_t3<- ifelse(cohort == 2 & partition==1, y_t3_g2,
                  ifelse(cohort == 3 & partition==1, y_t3_g3,
                         y_t3_never))
    #-----------------------------------------------------------------------------
    # Get the event study parameter
    # Get unfeasible ATT(g,t) (assuming we see the Potential outcomes)
    att_g2_t2_unf <- (mean((group_types==1)*y_t2_g2) - mean((group_types==1)*y_t2_never))/
      mean((group_types==1))
    att_g2_t3_unf <- (mean((group_types==1)*y_t3_g2) - mean((group_types==1)*y_t3_never))/
      mean((group_types==1))
    att_g3_t2_unf <- (mean((group_types==3)*y_t2_g3) - mean((group_types==3)*y_t2_never))/
      mean((group_types==3))
    att_g3_t3_unf <- (mean((group_types==3)*y_t3_g3) - mean((group_types==3)*y_t3_never))/
      mean((group_types==3))

    # Get unfeasible Event-Study (ES) (assuming we see the Potential outcomes and pscores)
    # get weights
    prob_g2_p1 <- mean(pi_2A/(pi_2A + pi_3A))
    prob_g3_p1 <- mean(pi_3A/(pi_2A + pi_3A))
    # Get ES
    ES_0_unf <- att_g2_t2_unf * prob_g2_p1 + att_g3_t3_unf * prob_g3_p1
    ES_1_unf <- att_g2_t3_unf
    # ------------------------------------------------------------------------------
    # Put all data in a data frame
    dta_wide <- data.frame(cbind(id = 1:size,
                                 y_t1, y_t2 , y_t3,
                                 cohort , partition,
                                 z)) # we always observe z
    # Transform data into long format
    dta <- get_long_data(list(as.vector(y_t1), as.vector(y_t2), as.vector(y_t3)), cohort, partition, z)
    setorder(dta, "id", "time")
    # generate clusters (there's no actual within-cluster correlation) with sample with replacement by id
    # get unique ids
    unique_ids <- unique(dta$id)
    # sample cluster values for each unique id
    clusters <- data.table(id = unique_ids, cluster = sample(1:50, size = length(unique_ids), replace = TRUE))
    # merge the sampled clusters back into the original data table
    dta <- merge(dta, clusters, by = "id", all.x = TRUE)
    #-----------------------------------------------------------------------------
    # Return the data
    return(list(data = as.data.table(dta),
                data_wide = as.data.table(dta_wide),
                ES_0_unf = ES_0_unf,
                prob_g2_p1 = prob_g2_p1,
                prob_g3_p1 = prob_g3_p1))

  } else if(dgp_type == 3){
    # Pscore depends on Z but Regressions depend on X
    # index for pscores

    # Generalize Propensity score (3 groups and two partitions)
    pi_2A <- exp(fps2(Xsi.ps, w1, z, 1.25))
    pi_2B <- exp(fps2(-Xsi.ps, w1, z, -0.5))
    pi_3A <- exp(fps2(Xsi.ps, w2, z, 2))
    pi_3B <- exp(fps2(-Xsi.ps, w2, z, -1.25))
    pi_0A <- exp(fps2(Xsi.ps, w3, z, -0.5))
    # Get the sum of the exponentials
    sum_pi <- 1 + pi_2A + pi_2B + pi_3A + pi_3B + pi_0A
    # Construct the generalized pscore
    pi_2A <- pi_2A/sum_pi
    pi_2B <- pi_2B/sum_pi
    pi_3A <- pi_3A/sum_pi
    pi_3B <- pi_3B/sum_pi
    pi_0A <- pi_0A/sum_pi
    pi_0B <- 1 - (pi_2A + pi_2B + pi_3A + pi_3B + pi_0A)

    # Bind them into a matrix
    probs_pscore <- cbind(pi_2A, pi_2B, pi_3A, pi_3B, pi_0A, pi_0B)

    # Sample group types (1 and 3 are the group of interest)
    group_types <- apply(probs_pscore, 1,
                         function(pvec) sample(seq(1,6), size=1, prob=pvec))

    # Translate group types into partitions and cohorts
    partition <- ifelse(group_types %in% c(1,3,5), 1, 0)
    cohort <- ifelse(group_types %in% c(1,2),
                     2,
                     ifelse(group_types %in% c(3,4),
                            3,
                            0))
    #-----------------------------------------------------------------------------
    # Now create the indexes for the potential outcomes
    #-----------------------------------------------------------------------------
    # Index for PO
    index_lin <- freg(b1, x)
    index_partition <- partition * index_lin
    # Index for unobserved heterogeneity
    index_unobs_het <- cohort * index_lin + index_partition
    #Index for Conditional Parallel Trends (this does not vary with treatment group)
    index_trend <- index_lin
    # Generate unobserved heterogeneity
    # (fixed effects that are correlated with X, partition, and cohort)
    v <- stats::rnorm(size, mean = index_unobs_het, sd = 1)
    # Index to violates the parallel trends assumption
    index_pt_violation <- v/10

    # Index for per-period ATT(g,t)
    index_att_g2 <- 10
    index_att_g3 <- 25
    #-----------------------------------------------------------------------------
    # Generate the potential outcomes
    # Baseline index
    baseline_t1 <- index_lin + index_partition + v

    #Generate untreated outcome at time 1 (that is observed for all units)
    y_t1 <- baseline_t1 + stats::rnorm(size)

    # Generate potential outcomes at time 2
    baseline_t2 <- baseline_t1 +
      index_pt_violation + #Violate the CPT
      index_trend #this is for the trend based on X

    y_t2_never <- baseline_t2 + stats::rnorm(size)
    y_t2_g2 <- baseline_t2 + stats::rnorm(size) +
      index_att_g2 * partition # This is the treatment effects for cohort 2
    y_t2_g3 <-  baseline_t2 + stats::rnorm(size)

    # Generate potential outcomes at time 3
    baseline_t3 <- baseline_t1 +
      2*index_trend + #this is for the trend based on X
      2*index_pt_violation #Violate the CPT

    y_t3_never <- baseline_t3 + stats::rnorm(size)
    y_t3_g2 <- baseline_t3 + stats::rnorm(size) +
      2*index_att_g2 * partition # This is the treatment effects for cohort 2
    y_t3_g3 <- baseline_t3 + stats::rnorm(size) +
      index_att_g3 * partition # This is the treatment effects for cohort 2
    #-----------------------------------------------------------------------------
    # Get the realized outcomes now
    # y_t2<- y_t2_g2*(cohort == 2) + y_t2_g3*(cohort==3) + y_t2_never*(cohort==0)
    # y_t3<- y_t3_g2*(cohort == 2) + y_t3_g3*(cohort==3) + y_t3_never*(cohort==0)

    y_t2<-  ifelse(cohort == 2 & partition==1, y_t2_g2, y_t2_never)
    y_t3<- ifelse(cohort == 2 & partition==1, y_t3_g2,
                  ifelse(cohort == 3 & partition==1, y_t3_g3,
                         y_t3_never))
    #-----------------------------------------------------------------------------
    # Get the event study parameter
    # Get unfeasible ATT(g,t) (assuming we see the Potential outcomes)
    att_g2_t2_unf <- (mean((group_types==1)*y_t2_g2) - mean((group_types==1)*y_t2_never))/
      mean((group_types==1))
    att_g2_t3_unf <- (mean((group_types==1)*y_t3_g2) - mean((group_types==1)*y_t3_never))/
      mean((group_types==1))
    att_g3_t2_unf <- (mean((group_types==3)*y_t2_g3) - mean((group_types==3)*y_t2_never))/
      mean((group_types==3))
    att_g3_t3_unf <- (mean((group_types==3)*y_t3_g3) - mean((group_types==3)*y_t3_never))/
      mean((group_types==3))

    # Get unfeasible Event-Study (ES) (assuming we see the Potential outcomes and pscores)
    # get weights
    prob_g2_p1 <- mean(pi_2A/(pi_2A + pi_3A))
    prob_g3_p1 <- mean(pi_3A/(pi_2A + pi_3A))
    # Get ES
    ES_0_unf <- att_g2_t2_unf * prob_g2_p1 + att_g3_t3_unf * prob_g3_p1
    ES_1_unf <- att_g2_t3_unf
    # ------------------------------------------------------------------------------
    # Put all data in a data frame
    dta_wide <- data.frame(cbind(id = 1:size,
                                 y_t1, y_t2 , y_t3,
                                 cohort , partition,
                                 z)) # we always observe z
    # Transform data into long format
    dta <- get_long_data(list(as.vector(y_t1), as.vector(y_t2), as.vector(y_t3)), cohort, partition, z)
    setorder(dta, "id", "time")
    # generate clusters (there's no actual within-cluster correlation) with sample with replacement by id
    # get unique ids
    unique_ids <- unique(dta$id)
    # sample cluster values for each unique id
    clusters <- data.table(id = unique_ids, cluster = sample(1:50, size = length(unique_ids), replace = TRUE))
    # merge the sampled clusters back into the original data table
    dta <- merge(dta, clusters, by = "id", all.x = TRUE)
    #-----------------------------------------------------------------------------
    # Return the data
    return(list(data = as.data.table(dta),
                data_wide = as.data.table(dta_wide),
                ES_0_unf = ES_0_unf,
                prob_g2_p1 = prob_g2_p1,
                prob_g3_p1 = prob_g3_p1))

  } else if (dgp_type == 4){
    # Both models are functions of X (misspecified)
    # index for pscores

    # Generalize Propensity score (3 groups and two partitions)
    pi_2A <- exp(fps2(Xsi.ps, w1, x, 1.25))
    pi_2B <- exp(fps2(-Xsi.ps, w1, x, -0.5))
    pi_3A <- exp(fps2(Xsi.ps, w2, x, 2))
    pi_3B <- exp(fps2(-Xsi.ps, w2, x, -1.25))
    pi_0A <- exp(fps2(Xsi.ps, w3, x, -0.5))
    # Get the sum of the exponentials
    sum_pi <- 1 + pi_2A + pi_2B + pi_3A + pi_3B + pi_0A
    # Construct the generalized pscore
    pi_2A <- pi_2A/sum_pi
    pi_2B <- pi_2B/sum_pi
    pi_3A <- pi_3A/sum_pi
    pi_3B <- pi_3B/sum_pi
    pi_0A <- pi_0A/sum_pi
    pi_0B <- 1 - (pi_2A + pi_2B + pi_3A + pi_3B + pi_0A)

    # Bind them into a matrix
    probs_pscore <- cbind(pi_2A, pi_2B, pi_3A, pi_3B, pi_0A, pi_0B)

    # Sample group types (1 and 3 are the group of interest)
    group_types <- apply(probs_pscore, 1,
                         function(pvec) sample(seq(1,6), size=1, prob=pvec))

    # Translate group types into partitions and cohorts
    partition <- ifelse(group_types %in% c(1,3,5), 1, 0)
    cohort <- ifelse(group_types %in% c(1,2),
                     2,
                     ifelse(group_types %in% c(3,4),
                            3,
                            0))
    #-----------------------------------------------------------------------------
    # Now create the indexes for the potential outcomes
    #-----------------------------------------------------------------------------
    # Index for PO
    index_lin <- freg(b1, x)
    index_partition <- partition * index_lin
    # Index for unobserved heterogeneity
    index_unobs_het <- cohort * index_lin + index_partition
    #Index for Conditional Parallel Trends (this does not vary with treatment group)
    index_trend <- index_lin
    # Generate unobserved heterogeneity
    # (fixed effects that are correlated with X, partition, and cohort)
    v <- stats::rnorm(size, mean = index_unobs_het, sd = 1)
    # Index to violates the parallel trends assumption
    index_pt_violation <- v/10

    # Index for per-period ATT(g,t)
    index_att_g2 <- 10
    index_att_g3 <- 25
    #-----------------------------------------------------------------------------
    # Generate the potential outcomes
    # Baseline index
    baseline_t1 <- index_lin + index_partition + v

    #Generate untreated outcome at time 1 (that is observed for all units)
    y_t1 <- baseline_t1 + stats::rnorm(size)

    # Generate potential outcomes at time 2
    baseline_t2 <- baseline_t1 +
      index_pt_violation + #Violate the CPT
      index_trend #this is for the trend based on X

    y_t2_never <- baseline_t2 + stats::rnorm(size)
    y_t2_g2 <- baseline_t2 + stats::rnorm(size) +
      index_att_g2 * partition # This is the treatment effects for cohort 2
    y_t2_g3 <-  baseline_t2 + stats::rnorm(size)

    # Generate potential outcomes at time 3
    baseline_t3 <- baseline_t1 +
      2*index_trend + #this is for the trend based on X
      2*index_pt_violation #Violate the CPT

    y_t3_never <- baseline_t3 + stats::rnorm(size)
    y_t3_g2 <- baseline_t3 + stats::rnorm(size) +
      2*index_att_g2 * partition # This is the treatment effects for cohort 2
    y_t3_g3 <- baseline_t3 + stats::rnorm(size) +
      index_att_g3 * partition # This is the treatment effects for cohort 2
    #-----------------------------------------------------------------------------
    # Get the realized outcomes now
    # y_t2<- y_t2_g2*(cohort == 2) + y_t2_g3*(cohort==3) + y_t2_never*(cohort==0)
    # y_t3<- y_t3_g2*(cohort == 2) + y_t3_g3*(cohort==3) + y_t3_never*(cohort==0)

    y_t2<-  ifelse(cohort == 2 & partition==1, y_t2_g2, y_t2_never)
    y_t3<- ifelse(cohort == 2 & partition==1, y_t3_g2,
                  ifelse(cohort == 3 & partition==1, y_t3_g3,
                         y_t3_never))
    #-----------------------------------------------------------------------------
    # Get the event study parameter
    # Get unfeasible ATT(g,t) (assuming we see the Potential outcomes)
    att_g2_t2_unf <- (mean((group_types==1)*y_t2_g2) - mean((group_types==1)*y_t2_never))/
      mean((group_types==1))
    att_g2_t3_unf <- (mean((group_types==1)*y_t3_g2) - mean((group_types==1)*y_t3_never))/
      mean((group_types==1))
    att_g3_t2_unf <- (mean((group_types==3)*y_t2_g3) - mean((group_types==3)*y_t2_never))/
      mean((group_types==3))
    att_g3_t3_unf <- (mean((group_types==3)*y_t3_g3) - mean((group_types==3)*y_t3_never))/
      mean((group_types==3))

    # Get unfeasible Event-Study (ES) (assuming we see the Potential outcomes and pscores)
    # get weights
    prob_g2_p1 <- mean(pi_2A/(pi_2A + pi_3A))
    prob_g3_p1 <- mean(pi_3A/(pi_2A + pi_3A))
    # Get ES
    ES_0_unf <- att_g2_t2_unf * prob_g2_p1 + att_g3_t3_unf * prob_g3_p1
    ES_1_unf <- att_g2_t3_unf
    # ------------------------------------------------------------------------------
    # Put all data in a data frame
    dta_wide <- data.frame(cbind(id = 1:size,
                                 y_t1, y_t2 , y_t3,
                                 cohort , partition,
                                 z)) # we always observe z
    # Transform data into long format
    dta <- get_long_data(list(as.vector(y_t1), as.vector(y_t2), as.vector(y_t3)), cohort, partition, z)
    setorder(dta, "id", "time")
    # generate clusters (there's no actual within-cluster correlation) with sample with replacement by id
    # get unique ids
    unique_ids <- unique(dta$id)
    # sample cluster values for each unique id
    clusters <- data.table(id = unique_ids, cluster = sample(1:50, size = length(unique_ids), replace = TRUE))
    # merge the sampled clusters back into the original data table
    dta <- merge(dta, clusters, by = "id", all.x = TRUE)
    #-----------------------------------------------------------------------------
    # Return the data
    return(list(data = as.data.table(dta),
                data_wide = as.data.table(dta_wide),
                ES_0_unf = ES_0_unf,
                prob_g2_p1 = prob_g2_p1,
                prob_g3_p1 = prob_g3_p1))

  } else {
    stop("Invalid DGP type")
  }

}
