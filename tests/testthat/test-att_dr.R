# Testing that results with analytical standard errors and bootstrap are comparable
test_that("multiplication works", {
  # generating dataset without errors
  test_panel = generate_test_panel()

  # ------------------------------
  # Performing tests
  # ------------------------------

  ddd_analytical <- ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                 pname = "partition", xformla = ~x1 + x2,
                  data = test_panel, control_group = NULL, base_period = NULL, est_method = "dr", learners = NULL, n_folds = NULL,
                  weightsname = NULL, boot = TRUE, nboot = 1000,
                  inffunc = FALSE, skip_data_checks = FALSE)

  ddd_boostrap <- ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                  pname = "partition", xformla = ~x1 + x2,
                  data = test_panel, control_group = NULL, base_period = NULL, est_method = "dr", learners = NULL, n_folds = NULL,
                  weightsname = NULL, boot = FALSE, nboot = NULL,
                  inffunc = FALSE, skip_data_checks = FALSE)

  # Check that point estimates are the same
  expect_equal(ddd_analytical$ATT, ddd_boostrap$ATT)

  # Check that standard errors are comparable
  expect_equal(ddd_analytical$se, ddd_boostrap$se, tolerance = 0.5)
})
