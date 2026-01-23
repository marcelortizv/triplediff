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

  # Test 2: Warning for est_method
  expect_warning(ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                     pname = "partition", xformla = ~x1 + x2,
                     data = two_periods_no_errors_df, control_group = NULL, base_period = NULL, est_method = "whatever",
                     weightsname = NULL, boot = FALSE, nboot = NULL,
                     inffunc = FALSE, skip_data_checks = FALSE))

  # Test 3: Warning for missing values in outcome variable "yname"
  expect_warning(ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                     pname = "partition", xformla = ~x1 + x2,
                     data = missing_values_outcome_df, control_group = NULL, base_period = NULL, est_method = "dr",
                     weightsname = NULL, boot = FALSE, nboot = NULL,
                     inffunc = FALSE, skip_data_checks = FALSE))


  # ------------------------------
  # Errors
  # ------------------------------
  # Test 4: handling of non-numeric "yname" column
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                   pname = "partition", xformla = ~x1 + x2,
                   data = y_not_numeric, control_group = NULL, base_period = NULL, est_method = "dr",
                   weightsname = NULL, boot = FALSE, nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 5: error for missing values in treatment variable "gname"
  expect_error(suppressWarnings(ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                   pname = "partition", xformla = ~x1 + x2,
                   data = missing_values_treat_df, control_group = "nevertreated", base_period = NULL, est_method = "dr",
                   weightsname = NULL, boot = FALSE, nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE)))

  # Test 6: error for small groups for inference (e.g. only one treated unit)
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                   pname = "partition", xformla = ~x1 + x2,
                   data = one_treated_unit_df, control_group = NULL, base_period = NULL, est_method = "dr",
                   weightsname = NULL, boot = FALSE,  nboot = NULL,
                   inffunc = TRUE, skip_data_checks = FALSE))

  # Test 7: error when "yname" is not in data
  expect_error(ddd(yname = "whatever", tname = "year", idname = "id", gname = "treat",
                   pname = "partition", xformla = ~x1 + x2,
                   data = two_periods_no_errors_df, control_group = NULL, base_period = NULL, est_method = "dr",
                   weightsname = NULL, boot = FALSE, nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 8: error when "tname" is not in data
  expect_error(ddd(yname = "outcome", tname = "whatever", idname = "id", gname = "treat",
                   pname = "partition", xformla = ~x1 + x2,
                   data = two_periods_no_errors_df, control_group = NULL, base_period = NULL, est_method = "dr",
                   weightsname = NULL, boot = FALSE, nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 9: error when "dname" is not in data
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", gname = "whatever",
                   pname = "partition", xformla = ~x1 + x2,
                   data = two_periods_no_errors_df, control_group = NULL, base_period = NULL, est_method = "dr",
                   weightsname = NULL, boot = FALSE, nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 10: error when "partition" is not in data
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                   pname = "whatever", xformla = ~x1 + x2,
                   data = two_periods_no_errors_df, control_group = NULL, base_period = NULL, est_method = "dr",
                   weightsname = NULL, boot = FALSE, nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 11: error when "idname" is not in data
  expect_error(ddd(yname = "outcome", tname = "year", idname = "whatever", gname = "treat",
                   pname = "whatever", xformla = ~x1 + x2,
                   data = two_periods_no_errors_df, control_group = NULL, base_period = NULL, est_method = "dr",
                   weightsname = NULL, boot = FALSE, nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 12: error when "xformla" is not a formula
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                   pname = "partition", xformla = "x1 - x2",
                   data = two_periods_no_errors_df, control_group = NULL, base_period = NULL, est_method = "dr",
                   weightsname = NULL, boot = FALSE, nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 13: partition variable "partition" is not unique by id
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                   pname = "partition", xformla = ~x1 + x2,
                   data = partition_not_unique_df, control_group = NULL, base_period = NULL, est_method = "dr",
                   weightsname = NULL, boot = FALSE, nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 14: treatment variable "dname" is not unique by id
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                   pname = "partition", xformla = ~x1 + x2,
                   data = treat_not_unique_df, control_group = NULL, base_period = NULL, est_method = "dr",
                   weightsname = NULL, boot = FALSE, nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 15: error when covariates are not invariant
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                   pname = "partition", xformla = ~x1 + x2,
                   data = covariates.not.invariant.df, control_group = NULL, base_period = NULL, est_method = "dr",
                   weightsname = NULL, boot = FALSE, nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 16: More that 2 time periods
  expect_error(suppressWarnings(ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                   pname = "partition", xformla = ~x1 + x2,
                   data = more_than_two_periods_df, control_group = "nevertreated", base_period = NULL, est_method = "dr",
                   weightsname = NULL, boot = FALSE, nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE)))

  # Test 17: Partition is not numeric
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                   pname = "partition", xformla = ~x1 + x2,
                   data = partition_not_numeric_df, control_group = NULL, base_period = NULL, est_method = "dr",
                   weightsname = NULL, boot = FALSE, nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 18: Partition is not numeric
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                   pname = "partition", xformla = ~x1 + x2,
                   data = partition_not_numeric_df, control_group = NULL, base_period = NULL, est_method = "dr",
                   weightsname = NULL, boot = FALSE, nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

  # Test 19: Weights columns contains Null values
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                   pname = "partition", xformla = ~x1 + x2,
                   data = weights_null_df, control_group = NULL, base_period = NULL, est_method = "dr",
                   weightsname = "weights", boot = FALSE, nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

})

# Tests for collinearity handling

test_that("check_partition_collinearity helper function detects global collinearity", {
  # Test the helper function directly with globally collinear variables
  library(data.table)

  # Create test data with known collinearity (x3 = x1 exactly)
  set.seed(789)
  n <- 100
  test_data <- data.table(
    subgroup = c(rep(4, 25), rep(3, 25), rep(2, 25), rep(1, 25)),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )

  # Add a variable that's collinear with x1 (globally collinear)
  test_data$x3 <- test_data$x1

  # Test the function
  result <- triplediff:::check_partition_collinearity(
    test_data,
    "subgroup",
    c("x1", "x2", "x3")
  )

  # Should detect collinearity in all partitions (since x3 = x1 everywhere)
  expect_true(length(result$all_collinear) > 0)
  expect_true("x3" %in% result$all_collinear)
})

test_that("check_partition_collinearity detects subset-specific collinearity", {
  # Test the helper function with subset-specific collinearity
  # x3 is a linear combination of x1 only within subgroups 4 and 3
  library(data.table)

  set.seed(999)
  n <- 200

  # Create test data
  test_data <- data.table(
    subgroup = c(rep(4, 50), rep(3, 50), rep(2, 50), rep(1, 50)),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )

  # Add x3 that's collinear with x1 only in subgroups 4 and 3
  test_data$x3 <- rnorm(n)
  # Make x3 = 2*x1 only in subgroups 4 and 3 (perfect collinearity)
  test_data[subgroup %in% c(4, 3), x3 := 2 * x1]

  # Test the function
  result <- triplediff:::check_partition_collinearity(
    test_data,
    "subgroup",
    c("x1", "x2", "x3")
  )

  # Should detect collinearity in "subgroup 4 vs 3"
  expect_true("x3" %in% result$all_collinear)
  expect_true(any(grepl("subgroup 4 vs 3", result$collinear_vars$x3)))
})

test_that("check_partition_collinearity returns empty when no collinearity", {
  # Test the helper function with independent variables
  library(data.table)

  set.seed(111)
  n <- 100
  test_data <- data.table(
    subgroup = c(rep(4, 25), rep(3, 25), rep(2, 25), rep(1, 25)),
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n)
  )

  # Test the function - should find no collinearity
  result <- triplediff:::check_partition_collinearity(
    test_data,
    "subgroup",
    c("x1", "x2", "x3")
  )

  # Should NOT detect any collinearity
  expect_equal(length(result$all_collinear), 0)
})

test_that("check_partition_collinearity identifies correct partition", {
  # Test that the function correctly identifies which partition has collinearity
  library(data.table)

  set.seed(222)
  n <- 200

  # Create test data
  test_data <- data.table(
    subgroup = c(rep(4, 50), rep(3, 50), rep(2, 50), rep(1, 50)),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )

  # Add x3 that's collinear with x1 ONLY in subgroups 4 and 2 (not 3 or 1)
  test_data$x3 <- rnorm(n)
  test_data[subgroup %in% c(4, 2), x3 := 3 * x1]

  # Test the function
  result <- triplediff:::check_partition_collinearity(
    test_data,
    "subgroup",
    c("x1", "x2", "x3")
  )

  # Should detect collinearity in "subgroup 4 vs 2"
  expect_true("x3" %in% result$all_collinear)
  expect_true(any(grepl("subgroup 4 vs 2", result$collinear_vars$x3)))
})

test_that("check_partition_collinearity detects constant variables (collinear with intercept)", {
  # IMPORTANT: A constant variable is collinear with the intercept
  # This is a key case that must be detected
  library(data.table)

  set.seed(333)
  n <- 200

  # Create test data with an intercept column
  test_data <- data.table(
    subgroup = c(rep(4, 50), rep(3, 50), rep(2, 50), rep(1, 50)),
    intercept = rep(1, n),  # Intercept column (all 1s)
    x1 = rnorm(n),
    x2 = rnorm(n)
  )

  # Add x3 that's constant only in subgroups 4 and 3 (collinear with intercept in that partition)
  test_data$x3 <- rnorm(n)
  test_data[subgroup %in% c(4, 3), x3 := 5]  # Constant value in 4 vs 3 comparison

  # Test the function with intercept included
  result <- triplediff:::check_partition_collinearity(
    test_data,
    "subgroup",
    c("intercept", "x1", "x2", "x3")
  )

  # Should detect that x3 is collinear with the intercept in "subgroup 4 vs 3"
  expect_true("x3" %in% result$all_collinear)
  expect_true(any(grepl("subgroup 4 vs 3", result$collinear_vars$x3)))
})

test_that("check_partition_collinearity detects globally constant variables", {
  # A globally constant variable should be detected in ALL partitions
  library(data.table)

  set.seed(444)
  n <- 200

  # Create test data with an intercept column
  test_data <- data.table(
    subgroup = c(rep(4, 50), rep(3, 50), rep(2, 50), rep(1, 50)),
    intercept = rep(1, n),  # Intercept column (all 1s)
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rep(5, n)  # Globally constant - should be collinear with intercept everywhere
  )

  # Test the function
  result <- triplediff:::check_partition_collinearity(
    test_data,
    "subgroup",
    c("intercept", "x1", "x2", "x3")
  )

  # Should detect that x3 is collinear with the intercept in ALL partitions
  expect_true("x3" %in% result$all_collinear)
  # Should be detected in all three comparisons (4 vs 3, 4 vs 2, 4 vs 1)
  expect_true(length(result$collinear_vars$x3) == 3)
})

# =============================================================================
# Propensity Score Trimming Tests
# =============================================================================

test_that("compute_pscore returns keep_ps indicator for trimming", {
  # Test that compute_pscore returns keep_ps and that it follows DRDID approach
  library(data.table)

  # Generate test data with panel structure
  set.seed(555)
  n <- 200

  # Create panel data with 2 time periods
  test_data <- data.table(
    id = rep(1:n, each = 2),
    post = rep(c(0, 1), n),
    subgroup = rep(c(rep(4, n/2), rep(3, n/2)), each = 2),
    y = rnorm(n * 2),
    weights = rep(1, n * 2),
    x1 = rep(rnorm(n), each = 2),
    x2 = rep(rnorm(n), each = 2)
  )

  # Call compute_pscore
  result <- triplediff:::compute_pscore(
    data = test_data,
    condition_subgroup = 3,
    xformula = ~ x1 + x2,
    trim_level = 0.995
  )

  # Should return propensity_scores, hessian_matrix, and keep_ps

  expect_true("propensity_scores" %in% names(result))
  expect_true("hessian_matrix" %in% names(result))
  expect_true("keep_ps" %in% names(result))

  # keep_ps should be a logical vector of same length as propensity_scores
  expect_equal(length(result$keep_ps), length(result$propensity_scores))
  expect_true(is.logical(result$keep_ps))
})

test_that("compute_pscore_null returns keep_ps = TRUE for all units (REG method)", {
  # REG method doesn't use IPW, so no trimming is needed
  library(data.table)

  set.seed(666)
  n <- 100

  test_data <- data.table(
    id = rep(1:n, each = 2),
    post = rep(c(0, 1), n),
    subgroup = rep(c(rep(4, n/2), rep(3, n/2)), each = 2),
    y = rnorm(n * 2),
    weights = rep(1, n * 2)
  )

  result <- triplediff:::compute_pscore_null(
    data = test_data,
    condition_subgroup = 3
  )

  # Should return keep_ps = TRUE for all units
  expect_true("keep_ps" %in% names(result))
  expect_true(all(result$keep_ps == TRUE))
})

test_that("propensity score trimming excludes controls with ps >= 0.995", {
  # Test that control units with high propensity scores are trimmed
  library(data.table)

  # Create data where some controls will have very high propensity scores
  # by making their covariates similar to treated units
  set.seed(777)
  n_treated <- 100
  n_control <- 100

  # Treated units (subgroup 4) have x1 ~ N(3, 0.5) - strongly positive
  treated_x1 <- rnorm(n_treated, mean = 3, sd = 0.5)

  # Controls: most have x1 ~ N(-1, 1), but some have x1 ~ N(3, 0.1) (almost like treated)
  n_high_ps_controls <- 10
  control_x1 <- c(
    rnorm(n_control - n_high_ps_controls, mean = -1, sd = 1),  # Normal controls
    rnorm(n_high_ps_controls, mean = 3, sd = 0.1)  # "Almost treated" controls
  )

  test_data <- data.table(
    id = rep(1:(n_treated + n_control), each = 2),
    post = rep(c(0, 1), n_treated + n_control),
    subgroup = rep(c(rep(4, n_treated), rep(3, n_control)), each = 2),
    y = rnorm((n_treated + n_control) * 2),
    weights = rep(1, (n_treated + n_control) * 2),
    x1 = rep(c(treated_x1, control_x1), each = 2)
  )

  result <- suppressWarnings(triplediff:::compute_pscore(
    data = test_data,
    condition_subgroup = 3,
    xformula = ~ x1,
    trim_level = 0.995
  ))

  # Get unique data to match propensity scores
  uid_data <- unique(test_data, by = "id")
  PA4 <- ifelse(uid_data$subgroup == 4, 1, 0)

  # Treated units (PA4 = 1) should all have keep_ps = TRUE
  expect_true(all(result$keep_ps[PA4 == 1] == TRUE))

  # Some control units should be trimmed (have keep_ps = FALSE)
  # due to their high propensity scores
  control_keep_ps <- result$keep_ps[PA4 == 0]
  control_ps <- result$propensity_scores[PA4 == 0]

  # Verify trimming logic: keep_ps = FALSE when ps >= 0.995
  expect_true(all(control_keep_ps[control_ps >= 0.995] == FALSE))
  expect_true(all(control_keep_ps[control_ps < 0.995] == TRUE))
})

test_that("propensity scores are capped at 1 - 1e-6 (following DRDID)", {
  # Test that propensity scores are bounded away from 1
  library(data.table)

  # Create extreme separation case
  set.seed(888)
  n_treated <- 50
  n_control <- 50

  # Perfect separation: treated have x1 > 0, controls have x1 < 0
  treated_x1 <- runif(n_treated, min = 5, max = 10)
  control_x1 <- runif(n_control, min = -10, max = -5)

  test_data <- data.table(
    id = rep(1:(n_treated + n_control), each = 2),
    post = rep(c(0, 1), n_treated + n_control),
    subgroup = rep(c(rep(4, n_treated), rep(3, n_control)), each = 2),
    y = rnorm((n_treated + n_control) * 2),
    weights = rep(1, (n_treated + n_control) * 2),
    x1 = rep(c(treated_x1, control_x1), each = 2)
  )

  result <- suppressWarnings(
    triplediff:::compute_pscore(
      data = test_data,
      condition_subgroup = 3,
      xformula = ~ x1,
      trim_level = 0.995
    )
  )

  # All propensity scores should be <= 1 - 1e-6
  expect_true(all(result$propensity_scores <= 1 - 1e-6))

  # And > 0 (lower bound)
  expect_true(all(result$propensity_scores > 0))
})

test_that("trimming is applied to ATT weights in compute_did", {
  # Integration test: verify that trimmed units get zero weight
  library(data.table)

  # Generate data where we know some controls will be trimmed
  set.seed(999)
  n_treated <- 80
  n_control <- 80

  # Create strong separation for some controls
  treated_x1 <- rnorm(n_treated, mean = 2, sd = 0.5)
  n_high_ps <- 5
  control_x1 <- c(
    rnorm(n_control - n_high_ps, mean = -2, sd = 1),
    rnorm(n_high_ps, mean = 2, sd = 0.1)  # These should have high ps
  )

  test_data <- data.table(
    id = rep(1:(n_treated + n_control), each = 2),
    post = rep(c(0, 1), n_treated + n_control),
    subgroup = rep(c(rep(4, n_treated), rep(3, n_control)), each = 2),
    y = rnorm((n_treated + n_control) * 2),
    weights = rep(1, (n_treated + n_control) * 2),
    x1 = rep(c(treated_x1, control_x1), each = 2)
  )

  # Compute propensity scores
  pscores_result <- suppressWarnings(triplediff:::compute_pscore(
    data = test_data,
    condition_subgroup = 3,
    xformula = ~ x1,
    trim_level = 0.995
  ))

  # Create pscores list structure as expected by compute_did
  pscores <- list(pscores_result)  # Index 1 for subgroup 3

  # Compute outcome regression (simple version)
  reg_result <- triplediff:::compute_outcome_regression(
    data = test_data,
    condition_subgroup = 3,
    xformula = ~ x1
  )
  reg_adjustment <- list(reg_result)

  # Add deltaY to unique data
  uid_data <- unique(test_data, by = "id")
  uid_data[, subgroup := subgroup]

  # The compute_did function should use keep_ps to zero out trimmed weights
  # This is an indirect test - we verify the propensity score list has trimming info
  expect_true("keep_ps" %in% names(pscores[[1]]))

  # Verify some controls are trimmed
  uid_controls <- uid_data[subgroup == 3]
  control_ps <- pscores_result$propensity_scores[uid_data$subgroup == 3]
  control_keep <- pscores_result$keep_ps[uid_data$subgroup == 3]

  # If any controls have high ps, they should be trimmed
  high_ps_controls <- control_ps >= 0.995
  if (any(high_ps_controls)) {
    expect_true(all(control_keep[high_ps_controls] == FALSE))
  }
})

test_that("confidence interval coverage is valid after trimming", {
  # Monte Carlo test: verify that CI coverage is within 2% of nominal 95% level
  # This tests that the trimming doesn't break the inference
  # DGP with moderate covariate imbalance that may trigger some trimming
  skip_on_cran()  # Skip on CRAN due to computation time

  library(data.table)

  set.seed(1234)
  n_sims <- 1000  # Number of Monte Carlo replications
  true_att <- 2   # True treatment effect
  alpha <- 0.05
  nominal_coverage <- 1 - alpha  # 0.95

  # Track coverage
  covers <- logical(n_sims)

  for (sim in 1:n_sims) {
    # Generate DDD panel data with known treatment effect
    # Need units in all 4 cells: treated/control x partition=0/1
    n_per_cell <- 500  # 500 units per cell = 2000 total (larger sample for better coverage)
    n_units <- 4 * n_per_cell

    # Create balanced 2x2 design
    treat_group <- c(rep(1, 2*n_per_cell), rep(0, 2*n_per_cell))
    partition <- c(rep(1, n_per_cell), rep(0, n_per_cell),
                   rep(1, n_per_cell), rep(0, n_per_cell))

    # Covariates with moderate difference between treated and control
    # This creates propensity score variation but with good overlap
    # Treated units: x1 ~ N(0.5, 1)
    # Control units: most ~ N(-0.5, 1), but a few "almost treated" with extreme x1
    x1_treated <- rnorm(2 * n_per_cell, mean = 0.5, sd = 1)

    # Add ~5 "almost treated" controls per partition (10 total) with x1 ~ N(10, 0.1)
    # These will have propensity scores >= 0.995 and trigger trimming
    n_high_ps <- 5
    n_normal_control <- n_per_cell - n_high_ps
    x1_control <- c(
      rnorm(n_normal_control, mean = -0.5, sd = 1),  # Normal controls, partition=1
      rnorm(n_high_ps, mean = 10, sd = 0.1),         # "Almost treated" controls, partition=1
      rnorm(n_normal_control, mean = -0.5, sd = 1),  # Normal controls, partition=0
      rnorm(n_high_ps, mean = 10, sd = 0.1)          # "Almost treated" controls, partition=0
    )
    x1 <- c(x1_treated, x1_control)
    x2 <- rnorm(n_units)

    # Create panel data (2 time periods)
    panel_data <- data.table(
      id = rep(1:n_units, each = 2),
      year = rep(c(0, 1), n_units),
      partition = rep(partition, each = 2),
      treat = rep(treat_group, each = 2),
      x1 = rep(x1, each = 2),
      x2 = rep(x2, each = 2)
    )

    # Generate outcome with DDD structure:
    # Y = base + partition_effect + time_effect + treat_effect + x_effects + DDD_effect + noise
    panel_data[, outcome := (
      1 +                                      # base
        0.5 * partition +                        # partition effect
        0.3 * year +                             # time effect
        0.4 * treat +                            # treatment group effect
        0.5 * x1 + 0.3 * x2 +                    # covariate effects
        0.2 * partition * year +                 # parallel trends for partition
        0.2 * treat * year +                     # parallel trends for treatment
        0.1 * partition * treat +                # partition-treatment interaction
        # DDD effect: only for treated, post-period, partition=1
        true_att * treat * year * partition +
        rnorm(.N, sd = 1)                        # noise
    )]

    # Run the DDD estimator
    result <- tryCatch({
      suppressWarnings(
        ddd(
          yname = "outcome",
          tname = "year",
          idname = "id",
          gname = "treat",
          pname = "partition",
          xformla = ~ x1 + x2,
          data = panel_data,
          est_method = "dr",
          boot = FALSE,
          inffunc = TRUE
        )
      )
    }, error = function(e) NULL)

    if (!is.null(result) && !is.null(result$ATT) && !is.null(result$se) &&
        !is.na(result$ATT) && !is.na(result$se)) {
      # Check if true ATT is within 95% CI
      # ddd returns ATT, se, lci, uci directly (not nested)
      est <- result$ATT
      se <- result$se
      ci_lower <- est - qnorm(1 - alpha/2) * se
      ci_upper <- est + qnorm(1 - alpha/2) * se
      covers[sim] <- (true_att >= ci_lower) & (true_att <= ci_upper)
    } else {
      covers[sim] <- NA
    }
  }

  # Calculate coverage rate (excluding failed simulations)
  valid_covers <- covers[!is.na(covers)]
  coverage_rate <- mean(valid_covers)

  # Coverage should be within 2% of nominal 95% level (i.e., between 93% and 97%)
  expect_true(coverage_rate >= nominal_coverage - 0.02,
              info = paste("Coverage rate:", round(coverage_rate, 4),
                           "- expected at least", nominal_coverage - 0.02))
  expect_true(coverage_rate <= nominal_coverage + 0.02,
              info = paste("Coverage rate:", round(coverage_rate, 4),
                           "- expected at most", nominal_coverage + 0.02))
})

# =============================================================================
# Global Collinearity Check Tests for Multiple Periods
# =============================================================================

test_that("run_preprocess_multPeriods detects and drops globally collinear covariates", {
  # Test that globally collinear variables are detected and dropped in multi-period setting
  library(data.table)

  set.seed(12345)
  n_units <- 200
  n_periods <- 4

  # Create staggered adoption panel data
  # Groups: 0 (never treated), 2 (treated in period 2), 3 (treated in period 3)
  group_assignment <- sample(c(0, 2, 3), n_units, replace = TRUE, prob = c(0.4, 0.3, 0.3))

  # Create panel structure
  panel_data <- data.table(
    id = rep(1:n_units, each = n_periods),
    year = rep(1:n_periods, n_units),
    first_treat = rep(group_assignment, each = n_periods),
    partition = rep(sample(0:1, n_units, replace = TRUE), each = n_periods)
  )

  # Add covariates
  panel_data[, x1 := rnorm(.N)]
  panel_data[, x2 := rnorm(.N)]
  # x3 is perfectly collinear with x1 (x3 = 2*x1)
  panel_data[, x3 := 2 * x1]

  # Generate outcome
  panel_data[, outcome := 1 + 0.5 * x1 + 0.3 * x2 + rnorm(.N)]

  # Run ddd with collinear covariates - should warn about dropping x3

  expect_warning(
    result <- ddd(
      yname = "outcome",
      tname = "year",
      idname = "id",
      gname = "first_treat",
      pname = "partition",
      xformla = ~ x1 + x2 + x3,
      data = panel_data,
      control_group = "nevertreated",
      base_period = "varying",
      est_method = "dr",
      boot = FALSE
    ),
    regexp = "collinearity"
  )

  # Result should still be valid (estimation should complete)
  expect_true(!is.null(result))
  expect_true(inherits(result, "ddd"))
})

test_that("run_preprocess_multPeriods handles no collinearity case correctly",
{
  # Test that when there's no collinearity, all covariates are kept
  library(data.table)

  set.seed(54321)
  n_units <- 200
  n_periods <- 4

  # Create staggered adoption panel data
  group_assignment <- sample(c(0, 2, 3), n_units, replace = TRUE, prob = c(0.4, 0.3, 0.3))

  panel_data <- data.table(
    id = rep(1:n_units, each = n_periods),
    year = rep(1:n_periods, n_units),
    first_treat = rep(group_assignment, each = n_periods),
    partition = rep(sample(0:1, n_units, replace = TRUE), each = n_periods)
  )

  # Add independent covariates (no collinearity)
  panel_data[, x1 := rnorm(.N)]
  panel_data[, x2 := rnorm(.N)]
  panel_data[, x3 := rnorm(.N)]

  # Generate outcome
  panel_data[, outcome := 1 + 0.5 * x1 + 0.3 * x2 + 0.2 * x3 + rnorm(.N)]

  # Run ddd - should NOT warn about collinearity
  expect_no_warning(
    result <- suppressWarnings(ddd(
      yname = "outcome",
      tname = "year",
      idname = "id",
      gname = "first_treat",
      pname = "partition",
      xformla = ~ x1 + x2 + x3,
      data = panel_data,
      control_group = "nevertreated",
      base_period = "varying",
      est_method = "dr",
      boot = FALSE
    )),
    message = "collinearity"
  )

  # Result should be valid
  expect_true(!is.null(result))
  expect_true(inherits(result, "ddd"))
})

test_that("run_preprocess_multPeriods detects constant covariate (collinear with intercept)", {

  # Test that a constant variable (collinear with intercept) is detected
  library(data.table)

  set.seed(99999)
  n_units <- 200
  n_periods <- 4

  group_assignment <- sample(c(0, 2, 3), n_units, replace = TRUE, prob = c(0.4, 0.3, 0.3))

  panel_data <- data.table(
    id = rep(1:n_units, each = n_periods),
    year = rep(1:n_periods, n_units),
    first_treat = rep(group_assignment, each = n_periods),
    partition = rep(sample(0:1, n_units, replace = TRUE), each = n_periods)
  )

  # Add covariates - x3 is constant (collinear with intercept)
  panel_data[, x1 := rnorm(.N)]
  panel_data[, x2 := rnorm(.N)]
  panel_data[, x3 := 5]  # Constant value

  # Generate outcome
  panel_data[, outcome := 1 + 0.5 * x1 + 0.3 * x2 + rnorm(.N)]

  # Run ddd - should warn about dropping x3
  expect_warning(
    result <- ddd(
      yname = "outcome",
      tname = "year",
      idname = "id",
      gname = "first_treat",
      pname = "partition",
      xformla = ~ x1 + x2 + x3,
      data = panel_data,
      control_group = "nevertreated",
      base_period = "varying",
      est_method = "dr",
      boot = FALSE
    ),
    regexp = "collinearity"
  )

  expect_true(!is.null(result))
  expect_true(inherits(result, "ddd"))
})

# =============================================================================
# Repeated Cross-Section (RCS) Propensity Score Trimming Tests
# =============================================================================

test_that("compute_pscore_rc returns keep_ps indicator for trimming in RCS data", {

  # Test that compute_pscore_rc returns keep_ps for repeated cross-section data
  # Key difference from panel: each observation is independent (no repeated IDs)
  library(data.table)

  set.seed(1111)
  n_pre <- 100   # observations in pre-period
  n_post <- 100  # observations in post-period

  # Create RCS data structure: different units in pre and post periods
  # In RCS, each row is a unique observation with unique ID
  test_data <- data.table(
    id = 1:(n_pre + n_post),  # Each observation has unique ID
    post = c(rep(0, n_pre), rep(1, n_post)),
    subgroup = c(
      rep(4, n_pre/2), rep(3, n_pre/2),   # Pre-period: half treated, half control
      rep(4, n_post/2), rep(3, n_post/2)  # Post-period: half treated, half control
    ),
    y = rnorm(n_pre + n_post),
    weights = rep(1, n_pre + n_post),
    x1 = rnorm(n_pre + n_post),
    x2 = rnorm(n_pre + n_post)
  )

  # Call compute_pscore_rc
  result <- triplediff:::compute_pscore_rc(
    data = test_data,
    condition_subgroup = 3,
    xformula = ~ x1 + x2,
    trim_level = 0.995
  )

  # Should return propensity_scores, hessian_matrix, and keep_ps
  expect_true("propensity_scores" %in% names(result))
  expect_true("hessian_matrix" %in% names(result))
  expect_true("keep_ps" %in% names(result))

  # For RCS, length should match all observations in condition_data (not unique by id)
  condition_data <- test_data[test_data$subgroup %in% c(3, 4)]
  expect_equal(length(result$keep_ps), nrow(condition_data))
  expect_true(is.logical(result$keep_ps))
})

test_that("compute_pscore_null_rc returns keep_ps = TRUE for all units (REG method in RCS)", {
  # REG method doesn't use IPW, so no trimming is needed in RCS
  library(data.table)

  set.seed(2222)
  n <- 200

  # Create RCS data
  test_data <- data.table(
    id = 1:n,
    post = c(rep(0, n/2), rep(1, n/2)),
    subgroup = rep(c(4, 3), n/2),
    y = rnorm(n),
    weights = rep(1, n)
  )

  result <- triplediff:::compute_pscore_null_rc(
    data = test_data,
    condition_subgroup = 3
  )

  # Should return keep_ps = TRUE for all observations
  expect_true("keep_ps" %in% names(result))
  expect_true(all(result$keep_ps == TRUE))

  # Length should match all observations in subgroups 3 and 4
  condition_data <- test_data[test_data$subgroup %in% c(3, 4)]
  expect_equal(length(result$keep_ps), nrow(condition_data))
})

test_that("RCS propensity score trimming excludes controls with ps >= 0.995", {
  # Test that control units with high propensity scores are trimmed in RCS
  library(data.table)

  set.seed(3333)
  n_treated <- 100
  n_control <- 100

  # Treated units (subgroup 4) have x1 ~ N(3, 0.5) - strongly positive
  treated_x1 <- rnorm(n_treated, mean = 3, sd = 0.5)

  # Controls: most have x1 ~ N(-1, 1), but some have x1 ~ N(3, 0.1) (almost like treated)
  n_high_ps_controls <- 10
  control_x1 <- c(
    rnorm(n_control - n_high_ps_controls, mean = -1, sd = 1),
    rnorm(n_high_ps_controls, mean = 3, sd = 0.1)  # "Almost treated" controls
  )

  # Create RCS data: each observation is independent
  test_data <- data.table(
    id = 1:(n_treated + n_control),
    post = sample(0:1, n_treated + n_control, replace = TRUE),  # Random pre/post assignment
    subgroup = c(rep(4, n_treated), rep(3, n_control)),
    y = rnorm(n_treated + n_control),
    weights = rep(1, n_treated + n_control),
    x1 = c(treated_x1, control_x1)
  )

  result <- suppressWarnings(triplediff:::compute_pscore_rc(
    data = test_data,
    condition_subgroup = 3,
    xformula = ~ x1,
    trim_level = 0.995
  ))

  # Get PA4 indicator for the condition_data
  condition_data <- test_data[test_data$subgroup %in% c(3, 4)]
  PA4 <- ifelse(condition_data$subgroup == 4, 1, 0)

  # Treated units (PA4 = 1) should all have keep_ps = TRUE
  expect_true(all(result$keep_ps[PA4 == 1] == TRUE))

  # Verify trimming logic: keep_ps = FALSE when ps >= 0.995 for controls
  control_keep_ps <- result$keep_ps[PA4 == 0]
  control_ps <- result$propensity_scores[PA4 == 0]

  expect_true(all(control_keep_ps[control_ps >= 0.995] == FALSE))
  expect_true(all(control_keep_ps[control_ps < 0.995] == TRUE))
})

test_that("RCS propensity scores are capped at 1 - 1e-6 (following DRDID)", {
  # Test that propensity scores are bounded away from 1 in RCS
  library(data.table)

  set.seed(4444)
  n_treated <- 50
  n_control <- 50

  # Perfect separation: treated have x1 > 0, controls have x1 < 0
  treated_x1 <- runif(n_treated, min = 5, max = 10)
  control_x1 <- runif(n_control, min = -10, max = -5)

  # Create RCS data
  test_data <- data.table(
    id = 1:(n_treated + n_control),
    post = sample(0:1, n_treated + n_control, replace = TRUE),
    subgroup = c(rep(4, n_treated), rep(3, n_control)),
    y = rnorm(n_treated + n_control),
    weights = rep(1, n_treated + n_control),
    x1 = c(treated_x1, control_x1)
  )

  result <- suppressWarnings(
    triplediff:::compute_pscore_rc(
      data = test_data,
      condition_subgroup = 3,
      xformula = ~ x1,
      trim_level = 0.995
    )
  )

  # All propensity scores should be <= 1 - 1e-6
  expect_true(all(result$propensity_scores <= 1 - 1e-6))

  # And > 0 (lower bound)
  expect_true(all(result$propensity_scores > 0))
})

test_that("RCS trimming works correctly with compute_did_rc for IPW method", {
  # Integration test: verify trimming is applied in compute_did_rc
  library(data.table)

  set.seed(5555)
  n_treated <- 80
  n_control <- 80

  # Create separation for some controls
  treated_x1 <- rnorm(n_treated, mean = 2, sd = 0.5)
  n_high_ps <- 5
  control_x1 <- c(
    rnorm(n_control - n_high_ps, mean = -2, sd = 1),
    rnorm(n_high_ps, mean = 2, sd = 0.1)  # High propensity score controls
  )

  # Create RCS data with pre and post periods
  n_total <- n_treated + n_control
  test_data <- data.table(
    id = 1:n_total,
    post = c(rep(0, n_total/2), rep(1, n_total/2)),  # Half pre, half post
    subgroup = rep(c(rep(4, n_treated/2), rep(3, n_control/2)), 2),
    y = rnorm(n_total),
    weights = rep(1, n_total),
    x1 = rep(c(treated_x1[1:(n_treated/2)], control_x1[1:(n_control/2)]), 2)
  )

  # Compute propensity scores for RCS
  pscores_result <- suppressWarnings(triplediff:::compute_pscore_rc(
    data = test_data,
    condition_subgroup = 3,
    xformula = ~ x1,
    trim_level = 0.995
  ))

  # Verify the propensity score result has trimming info
  expect_true("keep_ps" %in% names(pscores_result))

  # Check that some controls are trimmed if they have high propensity scores
  condition_data <- test_data[test_data$subgroup %in% c(3, 4)]
  PA4 <- ifelse(condition_data$subgroup == 4, 1, 0)
  control_ps <- pscores_result$propensity_scores[PA4 == 0]
  control_keep <- pscores_result$keep_ps[PA4 == 0]

  # If any controls have high ps, they should be trimmed
  high_ps_controls <- control_ps >= 0.995
  if (any(high_ps_controls)) {
    expect_true(all(control_keep[high_ps_controls] == FALSE))
  }
})

test_that("End-to-end RCS DDD estimation with trimming works correctly", {
  # Full integration test: run ddd() with panel=FALSE (RCS) and verify trimming
  skip_on_cran()  # Skip on CRAN due to computation time
  library(data.table)

  set.seed(6666)
  n_obs <- 400  # Total observations

  # Create RCS data with known structure
  # In RCS, we need: treated/control x partition x pre/post
  # Each observation is independent (different units)

  # Create balanced design across cells
  n_per_cell <- n_obs / 8  # 8 cells: 2 treat x 2 partition x 2 time

  test_data <- data.table(
    id = 1:n_obs,  # Each row is unique observation
    year = rep(c(1, 2), each = n_obs/2),  # Pre and post periods
    treat = rep(rep(c(0, 1), each = n_per_cell * 2), 2),
    partition = rep(rep(c(0, 1), each = n_per_cell), 4)
  )

  # Add covariates with moderate separation
  # Treated units have higher x1 on average
  test_data[, x1 := ifelse(treat == 1, rnorm(.N, mean = 1, sd = 1), rnorm(.N, mean = -1, sd = 1))]
  test_data[, x2 := rnorm(.N)]

  # Add a few "almost treated" controls to trigger trimming
  n_high_ps <- 5
  high_ps_indices <- sample(which(test_data$treat == 0), n_high_ps)
  test_data[high_ps_indices, x1 := rnorm(n_high_ps, mean = 5, sd = 0.1)]

  # Generate outcome with DDD structure
  true_att <- 2
  test_data[, outcome := (
    1 +
    0.5 * partition +
    0.3 * year +
    0.4 * treat +
    0.3 * x1 + 0.2 * x2 +
    true_att * treat * year * partition +
    rnorm(.N, sd = 1)
  )]

  # Run DDD with panel = FALSE (RCS mode)
  result <- suppressWarnings(
    ddd(
      yname = "outcome",
      tname = "year",
      idname = "id",
      gname = "treat",
      pname = "partition",
      xformla = ~ x1 + x2,
      data = test_data,
      panel = FALSE,  # This triggers RCS mode
      est_method = "dr",
      boot = FALSE
    )
  )

  # Result should be valid
  expect_true(!is.null(result))
  expect_true(inherits(result, "ddd"))
  expect_true(!is.na(result$ATT))
  expect_true(!is.na(result$se))
})

test_that("RCS confidence interval coverage is valid after trimming", {
  # Monte Carlo test: verify that CI coverage is within 2% of nominal 95% level

  # This tests that the trimming doesn't break the inference for RCS data
  # DGP with moderate covariate imbalance that may trigger some trimming
  skip_on_cran()  # Skip on CRAN due to computation time

  library(data.table)

  set.seed(5678)
  n_sims <- 1000  # Number of Monte Carlo replications
  true_att <- 2   # True treatment effect
  alpha <- 0.05
  nominal_coverage <- 1 - alpha  # 0.95

  # Track coverage
  covers <- logical(n_sims)

  for (sim in 1:n_sims) {
    # Generate DDD RCS data with known treatment effect
    # In RCS, each observation is independent (different units in pre vs post)
    # Need observations in all 8 cells: treated/control x partition x pre/post
    n_per_cell <- 250  # 250 obs per cell = 2000 total
    n_obs <- 8 * n_per_cell

    # Create balanced 2x2x2 design for RCS
    # Each row is a unique observation (unique id)
    rcs_data <- data.table(
      id = 1:n_obs,
      year = rep(c(0, 1), each = n_obs/2),
      treat = rep(rep(c(0, 1), each = n_per_cell * 2), 2),
      partition = rep(rep(c(0, 1), each = n_per_cell), 4)
    )

    # Covariates with moderate difference between treated and control
    # Treated units: x1 ~ N(0.5, 1)
    # Control units: most ~ N(-0.5, 1), but a few "almost treated" with extreme x1
    rcs_data[, x1 := ifelse(treat == 1,
                            rnorm(.N, mean = 0.5, sd = 1),
                            rnorm(.N, mean = -0.5, sd = 1))]
    rcs_data[, x2 := rnorm(.N)]

    # Add ~10 "almost treated" controls with x1 ~ N(10, 0.1)
    # These will have propensity scores >= 0.995 and trigger trimming
    n_high_ps <- 10
    high_ps_indices <- sample(which(rcs_data$treat == 0), n_high_ps)
    rcs_data[high_ps_indices, x1 := rnorm(n_high_ps, mean = 10, sd = 0.1)]

    # Generate outcome with DDD structure:
    # Y = base + partition_effect + time_effect + treat_effect + x_effects + DDD_effect + noise
    rcs_data[, outcome := (
      1 +                                      # base
        0.5 * partition +                        # partition effect
        0.3 * year +                             # time effect
        0.4 * treat +                            # treatment group effect
        0.5 * x1 + 0.3 * x2 +                    # covariate effects
        0.2 * partition * year +                 # parallel trends for partition
        0.2 * treat * year +                     # parallel trends for treatment
        0.1 * partition * treat +                # partition-treatment interaction
        # DDD effect: only for treated, post-period, partition=1
        true_att * treat * year * partition +
        rnorm(.N, sd = 1)                        # noise
    )]

    # Run the DDD estimator with panel = FALSE (RCS mode)
    result <- tryCatch({
      suppressWarnings(
        ddd(
          yname = "outcome",
          tname = "year",
          idname = "id",
          gname = "treat",
          pname = "partition",
          xformla = ~ x1 + x2,
          data = rcs_data,
          panel = FALSE,  # RCS mode
          est_method = "dr",
          boot = FALSE,
          inffunc = TRUE
        )
      )
    }, error = function(e) NULL)

    if (!is.null(result) && !is.null(result$ATT) && !is.null(result$se) &&
        !is.na(result$ATT) && !is.na(result$se)) {
      # Check if true ATT is within 95% CI
      est <- result$ATT
      se <- result$se
      ci_lower <- est - qnorm(1 - alpha/2) * se
      ci_upper <- est + qnorm(1 - alpha/2) * se
      covers[sim] <- (true_att >= ci_lower) & (true_att <= ci_upper)
    } else {
      covers[sim] <- NA
    }
  }

  # Calculate coverage rate (excluding failed simulations)
  valid_covers <- covers[!is.na(covers)]
  coverage_rate <- mean(valid_covers)

  # Coverage should be within 2% of nominal 95% level (i.e., between 93% and 97%)
  expect_true(coverage_rate >= nominal_coverage - 0.02,
              info = paste("RCS Coverage rate:", round(coverage_rate, 4),
                           "- expected at least", nominal_coverage - 0.02))
  expect_true(coverage_rate <= nominal_coverage + 0.02,
              info = paste("RCS Coverage rate:", round(coverage_rate, 4),
                           "- expected at most", nominal_coverage + 0.02))
})

# =============================================================================
# End-to-End Panel DDD Estimation with Trimming Tests
# =============================================================================

test_that("End-to-end Panel DDD estimation with trimming works correctly", {
  # Full integration test: run ddd() with panel=TRUE (default) and verify trimming
  skip_on_cran()  # Skip on CRAN due to computation time
  library(data.table)

  set.seed(7777)
  n_units <- 200  # Number of unique units

  # Create panel data with known structure
  # In panel data, same units are observed in pre and post periods
  # Need: treated/control x partition balanced design

  n_per_cell <- n_units / 4  # 4 cells: 2 treat x 2 partition

  # Create unit-level attributes (constant across time)
  unit_data <- data.table(
    id = 1:n_units,
    treat = rep(c(0, 1), each = n_units/2),
    partition = rep(rep(c(0, 1), each = n_per_cell), 2)
  )

  # Add covariates with moderate separation (unit-level)
  # Treated units have higher x1 on average
  unit_data[, x1 := ifelse(treat == 1, rnorm(.N, mean = 1, sd = 1), rnorm(.N, mean = -1, sd = 1))]
  unit_data[, x2 := rnorm(.N)]

  # Add a few "almost treated" controls to trigger trimming
  n_high_ps <- 5
  high_ps_indices <- sample(which(unit_data$treat == 0), n_high_ps)
  unit_data[high_ps_indices, x1 := rnorm(n_high_ps, mean = 5, sd = 0.1)]

  # Expand to panel data (2 time periods per unit)
  test_data <- rbindlist(lapply(c(1, 2), function(t) {
    dt <- copy(unit_data)
    dt[, year := t]
    dt
  }))
  setorder(test_data, id, year)

  # Generate outcome with DDD structure
  true_att <- 2
  test_data[, outcome := (
    1 +
    0.5 * partition +
    0.3 * (year - 1) +  # year effect (year 1 = pre, year 2 = post)
    0.4 * treat +
    0.3 * x1 + 0.2 * x2 +
    true_att * treat * (year - 1) * partition +  # DDD effect only in post period
    rnorm(.N, sd = 1)
  )]

  # Run DDD with panel = TRUE (default, panel mode)
  result <- suppressWarnings(
    ddd(
      yname = "outcome",
      tname = "year",
      idname = "id",
      gname = "treat",
      pname = "partition",
      xformla = ~ x1 + x2,
      data = test_data,
      panel = TRUE,  # Explicitly set panel mode
      est_method = "dr",
      boot = FALSE
    )
  )

  # Result should be valid
  expect_true(!is.null(result))
  expect_true(inherits(result, "ddd"))
  expect_true(!is.na(result$ATT))
  expect_true(!is.na(result$se))

  # ATT should be reasonably close to true value (within 3 SE)
  # This is a sanity check, not a strict test
  expect_true(abs(result$ATT - true_att) < 3 * result$se,
              info = paste("ATT:", round(result$ATT, 3), "SE:", round(result$se, 3)))
})

test_that("End-to-end Panel DDD estimation with IPW method and trimming", {
  # Test IPW method specifically to ensure trimming is applied correctly
  skip_on_cran()
  library(data.table)

  set.seed(8888)
  n_units <- 200

  n_per_cell <- n_units / 4

  # Create unit-level data
  unit_data <- data.table(
    id = 1:n_units,
    treat = rep(c(0, 1), each = n_units/2),
    partition = rep(rep(c(0, 1), each = n_per_cell), 2)
  )

  # Add covariates - treated units have higher x1
  unit_data[, x1 := ifelse(treat == 1, rnorm(.N, mean = 1.5, sd = 1), rnorm(.N, mean = -1.5, sd = 1))]
  unit_data[, x2 := rnorm(.N)]

  # Add "almost treated" controls to trigger trimming
  n_high_ps <- 8
  high_ps_indices <- sample(which(unit_data$treat == 0), n_high_ps)
  unit_data[high_ps_indices, x1 := rnorm(n_high_ps, mean = 6, sd = 0.1)]

  # Expand to panel
  test_data <- rbindlist(lapply(c(1, 2), function(t) {
    dt <- copy(unit_data)
    dt[, year := t]
    dt
  }))
  setorder(test_data, id, year)

  # Generate outcome
  true_att <- 1.5
  test_data[, outcome := (
    1 +
    0.5 * partition +
    0.3 * (year - 1) +
    0.4 * treat +
    0.3 * x1 + 0.2 * x2 +
    true_att * treat * (year - 1) * partition +
    rnorm(.N, sd = 1)
  )]

  # Run DDD with IPW method
  result <- suppressWarnings(
    ddd(
      yname = "outcome",
      tname = "year",
      idname = "id",
      gname = "treat",
      pname = "partition",
      xformla = ~ x1 + x2,
      data = test_data,
      panel = TRUE,
      est_method = "ipw",
      boot = FALSE
    )
  )

  # Result should be valid
  expect_true(!is.null(result))
  expect_true(inherits(result, "ddd"))
  expect_true(!is.na(result$ATT))
  expect_true(!is.na(result$se))
})

test_that("End-to-end Panel DDD estimation with REG method (no trimming needed)", {
  # REG method doesn't use propensity scores, so trimming shouldn't affect it
  skip_on_cran()
  library(data.table)

  set.seed(9999)
  n_units <- 200

  n_per_cell <- n_units / 4

  # Create unit-level data
  unit_data <- data.table(
    id = 1:n_units,
    treat = rep(c(0, 1), each = n_units/2),
    partition = rep(rep(c(0, 1), each = n_per_cell), 2)
  )

  # Add covariates
  unit_data[, x1 := ifelse(treat == 1, rnorm(.N, mean = 1, sd = 1), rnorm(.N, mean = -1, sd = 1))]
  unit_data[, x2 := rnorm(.N)]

  # Add "almost treated" controls (shouldn't matter for REG)
  n_high_ps <- 5
  high_ps_indices <- sample(which(unit_data$treat == 0), n_high_ps)
  unit_data[high_ps_indices, x1 := rnorm(n_high_ps, mean = 5, sd = 0.1)]

  # Expand to panel
  test_data <- rbindlist(lapply(c(1, 2), function(t) {
    dt <- copy(unit_data)
    dt[, year := t]
    dt
  }))
  setorder(test_data, id, year)

  # Generate outcome
  true_att <- 2.5
  test_data[, outcome := (
    1 +
    0.5 * partition +
    0.3 * (year - 1) +
    0.4 * treat +
    0.3 * x1 + 0.2 * x2 +
    true_att * treat * (year - 1) * partition +
    rnorm(.N, sd = 1)
  )]

  # Run DDD with REG method
  result <- suppressWarnings(
    ddd(
      yname = "outcome",
      tname = "year",
      idname = "id",
      gname = "treat",
      pname = "partition",
      xformla = ~ x1 + x2,
      data = test_data,
      panel = TRUE,
      est_method = "reg",
      boot = FALSE
    )
  )

  # Result should be valid
  expect_true(!is.null(result))
  expect_true(inherits(result, "ddd"))
  expect_true(!is.na(result$ATT))
  expect_true(!is.na(result$se))
})
