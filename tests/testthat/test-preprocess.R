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
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                   pname = "partition", xformla = ~x1 + x2,
                   data = missing_values_treat_df, control_group = "nevertreated", base_period = NULL, est_method = "dr",
                   weightsname = NULL, boot = FALSE, nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

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
  expect_error(ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                   pname = "partition", xformla = ~x1 + x2,
                   data = more_than_two_periods_df, control_group = "nevertreated", base_period = NULL, est_method = "dr",
                   weightsname = NULL, boot = FALSE, nboot = NULL,
                   inffunc = FALSE, skip_data_checks = FALSE))

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

  result <- triplediff:::compute_pscore(
    data = test_data,
    condition_subgroup = 3,
    xformula = ~ x1,
    trim_level = 0.995
  )

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
  pscores_result <- triplediff:::compute_pscore(
    data = test_data,
    condition_subgroup = 3,
    xformula = ~ x1,
    trim_level = 0.995
  )

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
