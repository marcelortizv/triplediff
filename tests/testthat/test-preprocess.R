# Testing error handling in run_preprocess_2periods() function
library(data.table)

test_that("Testing error handling in run_preprocess_2periods() function", {

  # generating dataset without errors
  two_periods_no_errors_df = generate_test_panel()

  # ------------------------------
  # Performing tests
  # ------------------------------

  # Introducing discrepancy in the dataset
  y_not_numeric = copy(two_periods_no_errors_df)
  y_not_numeric$outcome = as.character(y_not_numeric$outcome) # Introducing error: converting numeric to character

  # Introducing missing values in the outcome
  missing_values_outcome_df = copy(two_periods_no_errors_df)
  missing_values_outcome_df[1:2, "outcome"] = NA # Introducing missing values in the outcome variable

  # Introducing missing values in the treatment
  missing_values_treat_df = copy(two_periods_no_errors_df)
  missing_values_treat_df[1:2, "treat"] = NA # Introducing missing values in the treatment variable

  # Dataset only with 1 treated unit (inference is no feasible)
  one_treated_unit_df = copy(two_periods_no_errors_df)
  one_treated_unit_df$treat = ifelse(one_treated_unit_df$treat == 1, 0, 0)
  one_treated_unit_df[1:2, "treat"] = 1

  # partition is no unique by id
  partition_not_unique_df = copy(two_periods_no_errors_df)
  partition_not_unique_df$partition[1] <- ifelse(partition_not_unique_df$partition[1] == 1, 0, 1)

  # treatment variables "dname" is not unique by id
  treat_not_unique_df = copy(two_periods_no_errors_df)
  treat_not_unique_df$treat[1] <- ifelse(treat_not_unique_df$treat[1] == 1, 0, 1)

  # covariates are varying over time
  covariates_invariant_df = copy(two_periods_no_errors_df)
  covariates_invariant_df$x1 = rnorm(nrow(covariates_invariant_df))

  # more than 2 time periods
  more_than_two_periods_df = generate_test_panel(time = 3)

  # more than 2 groups
  more_than_two_groups_df = copy(two_periods_no_errors_df)
  more_than_two_groups_df$treat[1] = 2

  # partition is not numeric
  partition_not_numeric_df = copy(two_periods_no_errors_df)
  partition_not_numeric_df$partition = as.character(partition_not_numeric_df$partition)

  # partition is not binary
  partition_not_binary_df = copy(two_periods_no_errors_df)
  partition_not_binary_df$partition = rnorm(nrow(partition_not_binary_df))

  # Providing weights columns and one of them is null
  weights_null_df = copy(two_periods_no_errors_df)
  weights_null_df$weights = rep(1, nrow(weights_null_df))
  weights_null_df[weights_null_df$id == 1]$weights = NA

  # ------------------------------
  # Warnings
  # ------------------------------

  # Test 1: Warning for bootype
  expect_warning(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                     gname = NULL, partition_name = "partition", xformla = ~x1 + x2,
                     data = two_periods_no_errors_df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                     weightsname = NULL, boot = TRUE, boot_type = "whatever", nboot = NULL,
                     inffunc = FALSE, skip_data_checks = FALSE))

  # Test 2: Warning for est_method
  expect_warning(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                     gname = NULL, partition_name = "partition", xformla = ~x1 + x2,
                     data = two_periods_no_errors_df, control_group = NULL, base_period = NULL, est_method = "whatever", learners = NULL,
                     weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                     inffunc = FALSE, skip_data_checks = FALSE))

  # Test 3: Warning for missing values in outcome variable "yname"
  expect_warning(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                     gname = NULL, partition_name = "partition", xformla = ~x1 + x2,
                     data = missing_values_outcome_df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                     weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                     inffunc = FALSE, skip_data_checks = FALSE))


  # ------------------------------
  # Errors
  # ------------------------------
  # Test 4: handling of non-numeric "yname" column
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition_name = "partition", xformla = ~x1 + x2,
                   data = y_not_numeric, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 5: error for missing values in treatment variable "dname"
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition_name = "partition", xformla = ~x1 + x2,
                   data = missing_values_treat_df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 6: error for small groups for inference (e.g. only one treated unit)
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition_name = "partition", xformla = ~x1 + x2,
                   data = one_treated_unit_df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                   inffunc = TRUE, skip_data_checks = FALSE))

  # Test 7: error when "yname" is not in data
  expect_error(ddd(yname = "whatever", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition_name = "partition", xformla = ~x1 + x2,
                   data = two_periods_no_errors_df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 8: error when "tname" is not in data
  expect_error(ddd(yname = "outcome", tname = "whatever", idname = "id", dname = "treat",
                   gname = NULL, partition_name = "partition", xformla = ~x1 + x2,
                   data = two_periods_no_errors_df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 9: error when "dname" is not in data
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "whatever",
                   gname = NULL, partition_name = "partition", xformla = ~x1 + x2,
                   data = two_periods_no_errors_df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 10: error when "partition" is not in data
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition_name = "whatever", xformla = ~x1 + x2,
                   data = two_periods_no_errors_df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 11: error when "idname" is not in data
  expect_error(ddd(yname = "outcome", tname = "year", idname = "whatever", dname = "treat",
                   gname = NULL, partition_name = "whatever", xformla = ~x1 + x2,
                   data = two_periods_no_errors_df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 12: error when "xformla" is not a formula
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition_name = "partition", xformla = "x1 - x2",
                   data = two_periods_no_errors_df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 13: partition variable "partition" is not unique by id
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition_name = "partition", xformla = ~x1 + x2,
                   data = partition_not_unique_df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 14: treatment variable "dname" is not unique by id
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition_name = "partition", xformla = ~x1 + x2,
                   data = treat_not_unique_df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 15: error when covariates are not invariant
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition_name = "partition", xformla = ~x1 + x2,
                   data = covariates.not.invariant.df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 16: More that 2 time periods
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition_name = "partition", xformla = ~x1 + x2,
                   data = more_than_two_periods_df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 17: More than 2 groups
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = "group", partition_name = "partition", xformla = ~x1 + x2,
                   data = more_than_two_groups_df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 18: Partition is not numeric
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition_name = "partition", xformla = ~x1 + x2,
                   data = partition_not_numeric_df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 19: Partition is not numeric
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition_name = "partition", xformla = ~x1 + x2,
                   data = partition_not_numeric_df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                   weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 20: Weights columns contains Null values
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                   gname = NULL, partition_name = "partition", xformla = ~x1 + x2,
                   data = weights_null_df, control_group = NULL, base_period = NULL, est_method = "trad", learners = NULL,
                   weightsname = "weights", boot = FALSE, boot_type = "multiplier", nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

})
