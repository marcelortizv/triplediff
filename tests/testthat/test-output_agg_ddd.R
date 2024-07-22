# Testing if agg_ddd in generating output
test_that("Testing generation of output in aggregation function", {

  data <- gen_dgp_mult_periods(size = 1000, tperiods = 4, dgp_type = 1)

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

  expect_output(summary(simple.agg))
  expect_output(summary(es.agg))
  expect_output(summary(group.agg))
  expect_output(summary(calendar.agg))

})
