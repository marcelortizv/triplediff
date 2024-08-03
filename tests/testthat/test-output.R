# Testing if main function in generating output
test_that("Testing generation of output in main function", {

  # generating dataset without errors
  test_panel = generate_test_panel()
  # data from dgp 2 periods
  data_2periods <- gen_dgp_2periods(size = 1000, dgp_type = 1)[["data"]]
  # data from dgp multiple periods
  data_multperiods <- gen_dgp_mult_periods(size = 1000, dgp_type = 1)[["data"]]

  # ------------------------------
  # Performing tests
  # ------------------------------

  ddd_test <- ddd(yname = "outcome", tname = "year", idname = "id", gname = "treat",
                  pname = "partition", xformla = ~x1 + x2,
                  data = test_panel, control_group = NULL, base_period = NULL, est_method = "dr", learners = NULL,
                  weightsname = NULL, boot = FALSE, nboot = NULL,
                  inffunc = FALSE, skip_data_checks = FALSE)

  # performing DR estimatio with 2 periods
  ddd_dr = ddd(yname = "y", tname = "time", idname = "id", gname = "state",
               pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
               data = data_2periods, alpha = 0.05, est_method = "dr")

  # performing DR estimation with 2 periods and bootstrapping + clustered std. errors.
  ddd_boot_cluster = ddd(yname = "y", tname = "time", idname = "id", gname = "state",
                         pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
                         data = data_2periods, alpha = 0.05, boot = TRUE, nboot = 999,
                         est_method = "dr", cluster = "cluster", cband = TRUE)


  # performing DR estimation with multiple periods
  ddd_gt <- ddd(yname = "y", tname = "time", idname = "id",
                gname = "state", pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
                data = data_multperiods, control_group = "nevertreated", base_period = "universal",
                est_method = "dr")

  # performing DR estimation with multiple periods and bootstrapping
  ddd_gt_boot_cluster <- ddd(yname = "y", tname = "time", idname = "id",
                             gname = "state", pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
                             data = data_multperiods, control_group = "nevertreated", base_period = "universal",
                             est_method = "dr", boot = TRUE, nboot = 999, cluster = "cluster", cband = TRUE)

  # expecting results
  expect_output(summary(ddd_test))
  expect_output(print(ddd_test))
  expect_output(summary(ddd_dr))
  expect_output(summary(ddd_boot_cluster))
  expect_output(summary(ddd_gt))
  expect_output(summary(ddd_gt_boot_cluster))


})
