#' Function to generate a fake dataset for testing
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
