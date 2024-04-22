# Testing that results with analytical standard errors and bootstrap are comparable
test_that("multiplication works", {
  # generating dataset without errors
  test_panel = generate_test_panel()

  # ------------------------------
  # Performing tests
  # ------------------------------

  ddd_analytical <- ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                  gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                  data = test_panel, control.group = NULL, estMethod = "trad", learners = NULL,
                  weightsname = NULL, boot = TRUE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE)

  ddd_boostrap <- ddd(yname = "outcome", tname = "year", idname = "id", dname = "treat",
                  gname = NULL, partition.name = "partition", xformla = ~x1 + x2,
                  data = test_panel, control.group = NULL, estMethod = "trad", learners = NULL,
                  weightsname = NULL, boot = FALSE, boot.type = "multiplier", nboot = NULL, inffunc = FALSE)

  # Check that point estimates are the same
  expect_equal(ddd_analytical$ATT, ddd_boostrap$ATT)

  # Check that standard errors are comparable
  expect_equal(ddd_analytical$se, ddd_boostrap$se, tolerance = 0.5)
})
