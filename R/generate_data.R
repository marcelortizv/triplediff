#' Function to generate a fake dataset for testing
#' @import tidyr
#' @description
#' Function to generate fake dataset to test procedure
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

#' @export

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
    dt <- get_long_data(y0, y1, state, partition, z)

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
    dt <- get_long_data(y0, y1, state, partition, z)

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
    dt <- get_long_data(y0, y1, state, partition, z)

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
    dt <- get_long_data(y0, y1, state, partition, z)

    return(list(data = dt, att = 0, att.unf = att.unf, eff = eff))

  } else {
    stop("Invalid DGP type")
  }

}

# ---------------------------------------------------------------------
# Functions to generate data for multiple time periods
# ---------------------------------------------------------------------

# Define the function for each fps based on the given W and the period g
fps_function <- function(W, g, partition) {
  if (g == 0 && partition == 0) {
    return(0.2 * W)
  } else if (g > 0 && partition == 0) {
    return(0.15 * W)
  } else if (g == 0 && partition == 1) {
    return(0.18 * W)
  } else if (g > 0 && partition == 1) {
    return(0.04 * W)
  }
}

reorder_columns <- function(matrix) {
  # Function to reorder columns based on pattern "p_a_b"
  # Extract the 'a' and 'b' parts of the column names
  extract_g_l <- function(colname) {
    parts <- strsplit(colname, "_")[[1]]
    return(c(as.numeric(parts[2]), as.numeric(parts[3])))
  }

  # Create a data frame with column names and their corresponding 'a' and 'b' values
  col_info <- do.call(rbind, lapply(colnames(matrix), function(col) {
    c(col, extract_g_l(col))
  }))
  col_info <- as.data.frame(col_info, stringsAsFactors = FALSE)
  colnames(col_info) <- c("name", "g", "l")

  # Convert 'a' and 'b' to numeric for proper sorting
  col_info$g <- as.numeric(col_info$g)
  col_info$l <- as.numeric(col_info$l)

  # Order the columns by 'a' and then by 'b'
  ordered_col_info <- col_info[order(col_info$g, col_info$l), ]

  # Reorder the columns of the probability matrix
  reordered_matrix <- matrix[, ordered_col_info$name]
  return(list(reordered_pmatrix = reordered_matrix, p_names = ordered_col_info))
}

# Function to categorize based on cumulative probabilities
categorize <- function(U_row, prob_row) {
  cum_probs <- cumsum(prob_row)
  category <- sum(U_row > cum_probs) + 1
  return(category)
}

calculate_proba_g <- function(W, time_periods) {
  # gammaG <- c(0,0,1:time_periods)/(2*time_periods)
  gammaG <- c(0,1:(time_periods+1))/(2*(time_periods+1))
  # Calculate the exp(fps) for PA = 1 and PA = 2
  exp_un_a <- exp(gammaG[1] * apply(W, 1, fps_function, g = 0, partition = 0)) # (S = \infty, L = A)
  exp_un_b <- exp(gammaG[2] * apply(W, 1, fps_function, g = 0, partition = 1)) # (S = \infty, L = B)
  exp_tr_a <- list()
  exp_tr_b <- list()
  for (g in 1:time_periods) {
    exp_tr_a[[paste0("p_", g,"_0")]] <- exp(gammaG[g+2] * apply(W, 1, fps_function, g = g, partition = 0)) # (S = g, L = A)
    exp_tr_b[[paste0("p_", g,"_1")]] <- exp(gammaG[g+2] * apply(W, 1, fps_function, g = g, partition = 1)) # (S = g, L = B)
  }

  sum_exp <- exp_un_a + exp_un_a + Reduce(`+`, exp_tr_a) + Reduce(`+`, exp_tr_b)


  p_0_0 <- exp_un_a / sum_exp # P(S = \infty, L = A | X)
  p_0_1 <- exp_un_b / sum_exp # P(S = \infty, L = B | X)

  p_g_0 <- lapply(exp_tr_a, function(x) x / sum_exp)
  p_g_1 <- lapply(exp_tr_b, function(x) x / sum_exp)

  # Convert the list of divided elements into a matrix
  p_g_0_mat <- do.call(cbind, p_g_0)
  p_g_1_mat <- do.call(cbind, p_g_1)

  # Combine the original vector with the divided matrix
  p_matrix <- cbind(p_0_0, p_0_1, p_g_0_mat, p_g_1_mat)

  # adjusting probs to sum up to 1
  row_sums <- rowSums(p_matrix[, 1:(ncol(p_matrix) - 1)])
  p_matrix[, ncol(p_matrix)] <- 1 - row_sums

  return(p_matrix)
}

set_params <- function(time_periods){
  # time fixed effect
  thet <- seq(1:time_periods)
  # coefficient on X
  bett <- seq(1:time_periods)
  theu <- thet # changing this creates violations of parallel trends
  # covariate effect
  betu <- bett # changing this creates violations of conditional parallel trends

  te.t <- thet # no calendar time effects
  te.bet.ind <- rep(1,time_periods) # no selective treatment timing
  te.bet.X <- bett #no heterogeneous effects by X
  te.e <- rep(0,time_periods) # no dynamic effects
  te <- 1 # overall basic effect
  ge1 <- 2 # group fixed effect
  ge0 <- 1 # group fixed effect

  return(list(thet = thet,
              bett = bett,
              theu = theu,
              betu = betu,
              te.t = te.t,
              te.bet.ind = te.bet.ind,
              te.bet.X = te.bet.X,
              te.e = te.e,
              te = te,
              ge1 = ge1,
              ge0 = ge0))

}

# ---------------------------------------------------------------------
# Generate data for multiple treatment periods
# ---------------------------------------------------------------------

#' Function that generates panel data with staggered treatment assignment for multiple periods
#' @description
#' Function to generate data with staggered treatment adoption.
#'
#' @param size number of units
#' @param tperiods number of periods
#' @param dgp_type type of DGP to generate.
#'           1 if both nuisance functions are correct,
#'           2 if only the outcome model is correct,
#'           3 if only the pscore is correct,
#'           4 if both nuisance functions are incorrect
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

#' @export

gen_dgp_mult_periods <- function(size, tperiods, dgp_type){

  # Simulate W (covariate)
  W <- matrix(rnorm(size), ncol = 1)

  if (dgp_type == 1){
    X_m = W
    X_ps = W
  } else if (dgp_type == 2){
    X_m = W
    X_ps = (pnorm(W) + 0.5)^2
  } else if (dgp_type == 3){
    X_m = W^2
    X_ps = W
  } else if (dgp_type == 4){
    X_m = W^2
    X_ps = (pnorm(W) + 0.5)^2
  }

  params <- set_params(tperiods)
  thet = params$thet
  bett = params$bett
  theu = params$theu
  betu = params$betu
  te.t = params$te.t
  te.bet.ind = params$te.bet.ind
  te.bet.X = params$te.bet.X
  te.e = params$te.e
  te = params$te
  ge1 = params$ge1
  ge0 = params$ge0

  # Calculate probabilities for the constant partitions
  probs <- calculate_proba_g(X_ps, time_periods = tperiods)

  # Reorder columns
  probs_ordered <- reorder_columns(probs)
  probs <- probs_ordered$reordered_pmatrix
  p_names <- probs_ordered$p_names

  # Uniform random variable U
  U <- runif(size)

  # Apply categorization for each row
  PA <- sapply(1:size, function(i) {
    categorize(U[i], probs[i, ])
  })

  # Create G and L variables based on PA
  G <- sapply(PA, function(pa) p_names$g[pa])
  L <- sapply(PA, function(pa) p_names$l[pa])

  GL <- G*L

  GLt <- GL[GL>0]
  nt <- length(GLt)
  Xt <- X_m[GL>0]

  # draw individual fixed effect
  Ct <- rnorm(nt, mean=GL)

  # ---------------------------------------------------------
  # GENERATING POTENTIAL OUTCOMES IN EACH TIME PERIODS
  # ---------------------------------------------------------

  # Untreated potential outcomes for treated group
  Y0tmat <- sapply(1:tperiods, function(t) {
    thet[t] + ge1 + Ct + Xt*bett[t] + rnorm(nt)
  })
  Y0tdf <- as.data.frame(Y0tmat)

  # Generate treated potential outcomes for treated group
  Y1tdf <- sapply(1:tperiods, function(t) {
    te.t[t] + ge1 + te.bet.ind[GLt]*Ct + Xt*te.bet.X[t] + (GLt <= t)*te.e[sapply(1:nt, function(i) max(t-GLt[i]+1,1))] + te + rnorm(nt)
  })

  # generate observed data
  Ytdf <- sapply(1:tperiods, function(t) {
    (GLt<=t)*Y1tdf[,t] + (GLt>t)*Y0tdf[,t]
  })

  Ynames <- paste0("Y",1:tperiods)
  colnames(Ytdf) <- Ynames

  # store observed data for treated group
  dft <- cbind.data.frame(G=GLt, L=1, X=W[GL>0], Ytdf)

  # ---------------------------------------------------------
  # GENERATE DATA FOR CONTROL GROUPS
  # ---------------------------------------------------------

  # Control group because not in partiton
  nt_0 <- sum((G > 0)*(L==0))
  Ct_0 <- rnorm(nt_0, mean=0)
  Xt_0 <- X_m[as.logical((G>0)*(L==0))]


  Y0t_0_mat <- sapply(1:tperiods, function(t) {
    thet[t] + Ct_0 + ge0 + Xt_0*bett[t] + rnorm(nt_0)
  })
  Y0t_0df <- as.data.frame(Y0t_0_mat)
  colnames(Y0t_0df) <- Ynames


  # Control group because never treated
  nu_0 <- sum((G == 0)*(L==0))
  nu_1 <- sum((G == 0)*(L==1))

  Xu_0 <- X_m[as.logical((G==0)*(L==0))]
  Xu_1 <- X_m[as.logical((G==0)*(L==1))]

  # draw untreated fixed effect
  Cu_0 <- rnorm(nu_0, mean=0)
  Cu_1 <- rnorm(nu_1, mean=0)

  # generate untreated potential outcomes
  Y0u_0_mat <- sapply(1:tperiods, function(t) {
    theu[t] + Cu_0 + Xu_0*betu[t] + rnorm(nu_0)
  })
  Y0u_0df <- as.data.frame(Y0u_0_mat)
  colnames(Y0u_0df) <- Ynames

  Y0u_1_mat <- sapply(1:tperiods, function(t) {
    theu[t] + Cu_1 + Xu_1*betu[t] + rnorm(nu_1)
  })
  Y0u_1df <- as.data.frame(Y0u_1_mat)
  colnames(Y0u_1df) <- Ynames

  # store observed data for control group
  dfu_1 <- cbind.data.frame(G=0, L = 1, X = W[as.logical((G==0)*(L==1))], Y0u_1df)
  dfu_0 <- cbind.data.frame(G=0, L = 0, X = W[as.logical((G==0)*(L==0))], Y0u_0df)
  dft_0 <- cbind.data.frame(G=G[as.logical((G>0)*(L==0))], L = 0, X = W[as.logical((G>0)*(L==0))], Y0t_0df)

  # store overall dataset
  df <- rbind.data.frame(dft, dft_0, dfu_1, dfu_0)

  # generate id variable
  df$id <- 1:nrow(df)
  # generate clusters (there's no actual within-cluster correlation)
  df$cluster <- sample(1:50, size=nrow(df), replace=TRUE)

  ddf <- tidyr::pivot_longer(df,
                             cols=tidyr::starts_with("Y"),
                             names_to="period",
                             names_prefix="Y",
                             values_to="Y")
  ddf$period <- as.numeric(ddf$period)
  ddf <- ddf[order(ddf$id, ddf$period),] # reorder data
  ddf <- subset(ddf, G != 1)

  return(ddf)
}


