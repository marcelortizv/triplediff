# Testing error handling in run_preprocess_2periods() function

library(data.table)

test_that("Testing error handling in run_preprocess_2periods() function", {

  generate_test_panel <- function(seed = 123,
                                  num_ids = 50,
                                  time = 2,
                                  initial.year = 2019,
                                  treatment.year = 2020) {
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
      outcome = y
    )
    return(dt)
  }

  # generating dataset without errors
  two.periods.no.errors.df = generate_test_panel()

  # ------------------------------
  # Performing tests
  # ------------------------------

  # Introducing discrepancy in the dataset
  y.not.numeric = copy(two.periods.no.errors.df)
  y.not.numeric$outcome = as.character(y.not.numeric$outcome) # Introducing error: converting numeric to character

  # Introducing missing values in the outcome
  missing.values.outcome.df = copy(two.periods.no.errors.df)
  missing.values.outcome.df[1:2, "outcome"] = NA # Introducing missing values in the outcome variable

  # Introducing missing values in the treatment
  missing.values.treat.df = copy(two.periods.no.errors.df)
  missing.values.treat.df[1:2, "treat"] = NA # Introducing missing values in the treatment variable

  # Dataset only with 1 treated unit (inference is no feasible)
  one.treated.unit.df = copy(two.periods.no.errors.df)
  one.treated.unit.df$treat = ifelse(one.treated.unit.df$treat == 1, 0, 0)
  one.treated.unit.df[1:2, "treat"] = 1

  # partition is no unique by id
  partition.not.unique.df = copy(two.periods.no.errors.df)
  partition.not.unique.df$partition[1] <- ifelse(partition.not.unique.df$partition[1] == 1, 0, 1)

  # treatment variables "dname" is not unique by id
  treat.not.unique.df = copy(two.periods.no.errors.df)
  treat.not.unique.df$treat[1] <- ifelse(treat.not.unique.df$treat[1] == 1, 0, 1)

  # covariates are varying over time
  covariates.invariant.df = copy(two.periods.no.errors.df)
  covariates.invariant.df$x1 = rnorm(nrow(covariates.invariant.df))

  # more than 2 time periods
  more.than.two.periods.df = generate_test_panel(time = 3)

  # more than 2 groups
  more.than.two.groups.df = copy(two.periods.no.errors.df)
  more.than.two.groups.df$treat[1] = 2

  # partition is not numeric
  partition.not.numeric.df = copy(two.periods.no.errors.df)
  partition.not.numeric.df$partition = as.character(partition.not.numeric.df$partition)

  # partition is not binary
  partition.not.binary.df = copy(two.periods.no.errors.df)
  partition.not.binary.df$partition = rnorm(nrow(partition.not.binary.df))

  # ------------------------------
  # Warnings
  # ------------------------------

  # Test 1: Warning for bootype
  expect_warning(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                     gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                     data = two.periods.no.errors.df, control.group = NULL, estMethod = "trad", learners = NULL,
                     weightsname = NULL, boot = TRUE, boot.type = "whatever", nboot = NULL, inffunc = FALSE))

  # Test 2: Warning for estMethod
  expect_warning(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                     gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                     data = two.periods.no.errors.df, control.group = NULL, estMethod = "whatever", learners = NULL,
                     weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE))

  # Test 3: Warning for missing values in outcome variable "yname"
  expect_warning(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                     gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                     data = missing.values.outcome.df, control.group = NULL, estMethod = "trad", learners = NULL,
                     weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE))


  # ------------------------------
  # Errors
  # ------------------------------
  # Test 4: handling of non-numeric "yname" column
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                   data = y.not.numeric, control.group = NULL, estMethod = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE))

  # Test 5: error for missing values in treatment variable "dname"
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                   data = missing.values.treat.df, control.group = NULL, estMethod = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE))

  # Test 6: error for small groups for inference (e.g. only one treated unit)
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                   data = one.treated.unit.df, control.group = NULL, estMethod = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = TRUE))

  # Test 7: error when "yname" is not in data
  expect_error(ddd(yname = "whatever", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                   data = two.periods.no.errors.df, control.group = NULL, estMethod = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE))

  # Test 8: error when "tname" is not in data
  expect_error(ddd(yname = "outcome", tname = "whatever", idname = "id", dname = "treat",
                   gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                   data = two.periods.no.errors.df, control.group = NULL, estMethod = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE))

  # Test 9: error when "dname" is not in data
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "whatever",
                   gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                   data = two.periods.no.errors.df, control.group = NULL, estMethod = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE))

  # Test 10: error when "partition" is not in data
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition.name = "whatever", xformla = ~x1 + x2,
                   data = two.periods.no.errors.df, control.group = NULL, estMethod = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE))

  # Test 11: error when "idname" is not in data
  expect_error(ddd(yname = "outcome", tname = "year", idname = "whatever", dname = "treat",
                   gname = NULL, partition.name = "whatever", xformla = ~x1 + x2,
                   data = two.periods.no.errors.df, control.group = NULL, estMethod = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE))

  # Test 12: error when "xformla" is not a formula
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition.name = "partition", xformla = "x1 - x2",
                   data = two.periods.no.errors.df, control.group = NULL, estMethod = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE))

  # Test 13: partition variable "partition" is not unique by id
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                   data = partition.not.unique.df, control.group = NULL, estMethod = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE))

  # Test 14: treatment variable "dname" is not unique by id
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                   data = treat.not.unique.df, control.group = NULL, estMethod = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE))

  # Test 15: error when covariates are not invariant
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                   data = covariates.not.invariant.df, control.group = NULL, estMethod = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE))

  # Test 16: More that 2 time periods
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                   data = more.than.two.periods.df, control.group = NULL, estMethod = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE))

  # Test 17: More than 2 groups
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = "group", partition.name = "partition", xformla = ~x1 + x2,
                   data = more.than.two.groups.df, control.group = NULL, estMethod = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE))

  # Test 18: Partition is not numeric
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                   data = partition.not.numeric.df, control.group = NULL, estMethod = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE))

  # Test 19: Partition is not numeric
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                   data = partition.not.numeric.df, control.group = NULL, estMethod = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE))

})
