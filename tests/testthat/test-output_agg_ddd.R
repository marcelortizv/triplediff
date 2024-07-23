# Testing if agg_ddd in generating output
test_that("Testing generation of output in aggregation function", {

  data <- gen_dgp_mult_periods(size = 1000, tperiods = 4, dgp_type = 1)

  # Performing simple tests
  out <- ddd(yname = "Y", tname = "period", idname = "id", dname = NULL,
              gname = "G", partition_name = "L", xformla = ~X,
              data = data, control_group = "nevertreated", base_period = "varying",
              est_method = "dr")
  # Simple aggregation
  simple.agg <- agg_ddd(out, type = "simple", alpha = 0.10)

  # Event study aggregation
  es.agg <- agg_ddd(out, type = "eventstudy", alpha = 0.10)

  # Group aggregation
  group.agg <- agg_ddd(out, type = "group", alpha = 0.10)

  # Calendar aggregation
  calendar.agg <- agg_ddd(out, type = "calendar", alpha = 0.10)



  # Performing tests with bootstrap + clustered standard errors
  ddd_gt_boot_cluster <- ddd(yname = "Y", tname = "period", idname = "id", dname = NULL,
                             gname = "G", partition_name = "L", xformla = ~X,
                             data = data, control_group = "nevertreated", base_period = "varying",
                             est_method = "dr", boot = TRUE, nboot = 999, cluster = "cluster", cband = TRUE)

  # Simple aggregation boot
  simple.agg.bott <- agg_ddd(ddd_gt_boot_cluster, type = "simple", alpha = 0.10)

  # Event study aggregation boot
  es.agg.boot <- agg_ddd(ddd_gt_boot_cluster, type = "eventstudy", alpha = 0.10)

  # Group aggregation boot
  group.agg.boot <- agg_ddd(ddd_gt_boot_cluster, type = "group", alpha = 0.10)

  # Calendar aggregation boot
  calendar.agg.boot <- agg_ddd(ddd_gt_boot_cluster, type = "calendar", alpha = 0.10)


  # expecting results
  expect_output(summary(simple.agg))
  expect_output(summary(es.agg))
  expect_output(summary(group.agg))
  expect_output(summary(calendar.agg))
  expect_output(summary(simple.agg.bott))
  expect_output(summary(es.agg.boot))
  expect_output(summary(group.agg.boot))
  expect_output(summary(calendar.agg.boot))

})
