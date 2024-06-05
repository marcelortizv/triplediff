# Testing if main function in generating output
test_that("Testing generation of output in main function", {

  # generating dataset without errors
  test_panel = generate_test_panel()

  # ------------------------------
  # Performing tests
  # ------------------------------

  ddd_test <- ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                  gname = NULL, partition_name = "partition", xformla = ~x1 + x2,
                  data = test_panel, control_group = NULL, est_method = "trad", learners = NULL,
                  weightsname = NULL, boot = FALSE, boot_type = "multiplier", nboot = NULL,
                  inffunc = FALSE, skip_data_checks = FALSE)

  expect_output(summary(ddd_test))
  expect_output(print(ddd_test))

})
