# Testing that results with analytical standard errors and bootstrap are comparable
test_that("multiplication works", {
  # generating dataset without errors
  test_panel = generate_test_panel()

  # ------------------------------
  # Performing tests
  # ------------------------------

  ddd_boostrap <- ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                 pname = "partition", xformla = ~x1 + x2,
                  data = test_panel, control_group = NULL, base_period = NULL, est_method = "dr",
                  weightsname = NULL, boot = TRUE, nboot = 1000, cband = TRUE,
                  inffunc = FALSE, skip_data_checks = FALSE)

  ddd_analytical <- ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                  pname = "partition", xformla = ~x1 + x2,
                  data = test_panel, control_group = NULL, base_period = NULL, est_method = "dr",
                  weightsname = NULL, boot = FALSE, nboot = NULL,
                  inffunc = FALSE, skip_data_checks = FALSE)

  # Check that point estimates are the same
  expect_equal(ddd_analytical$ATT, ddd_boostrap$ATT)

  # Check that standard errors are comparable
  expect_equal(ddd_analytical$se, ddd_boostrap$se, tolerance = 0.5)
})

# Testing clustered standard error is working correctly
test_that("clustered standard errors are working correctly", {
  # generating dataset without errors
  test_panel = gen_dgp_2periods(size = 5000, dgp_type = 1)$data

  # ------------------------------
  # Performing tests
  # ------------------------------

  att_nocluster <- ddd(yname = "y", tname = "time", idname = "id", gname = "state",
                pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4, base_period = "universal",
                data = test_panel, control_group = "nevertreated", est_method = "dr")

  att_cluster <-  ddd(yname = "y", tname = "time", idname = "id", gname = "state",
                      pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
                      data = test_panel, control_group = "nevertreated", boot = TRUE, nboot = 1000, cband = TRUE,
                      base_period = "universal", est_method = "dr", cluster = "cluster")

  # Check that point estimates are the same
  expect_equal(att_nocluster$ATT, att_cluster$ATT)

  # Check that standard errors are different
  expect_false(isTRUE(all.equal(att_nocluster$se, att_cluster$se)))
})
