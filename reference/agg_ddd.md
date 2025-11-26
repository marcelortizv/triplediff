# Aggregate Group-Time Average Treatment Effects in Staggered Triple-Differences Designs.

`agg_ddd` is a function that take group-time average treatment effects
and aggregate them into a smaller number of summary parameters in
staggered triple differences designs. There are several possible
aggregations including `"simple"`, `"eventstudy"`, `"group"`, and
`"calendar"`. Default is `"eventstudy"`.

## Usage

``` r
agg_ddd(
  ddd_obj,
  type = "eventstudy",
  balance_e = NULL,
  min_e = -Inf,
  max_e = Inf,
  na.rm = FALSE,
  boot = NULL,
  nboot = NULL,
  cband = NULL,
  alpha = 0.05
)
```

## Arguments

- ddd_obj:

  a `ddd` object (i.e., the results of the
  [`ddd()`](http://marcelortiz.com/triplediff/reference/ddd.md)
  function)

- type:

  Which type of aggregated treatment effect parameter to compute.
  `"simple"` just computes a weighted average of all group-time average
  treatment effects with weights proportional to group size.
  `"eventstudy"` computes average effects across different lengths of
  exposure to the treatment (event times). Here the overall effect
  averages the effect of the treatment across the positive lengths of
  exposure. This is the default option; `"group"` computes average
  treatment effects across different groups/cohorts; here the overall
  effect averages the effect across different groups using group size as
  weights; `"calendar"` computes average treatment effects across
  different time periods, with weights proportional to the group size;
  here the overall effect averages the effect across each time period.

- balance_e:

  If set (and if one computes event study), it balances the sample with
  respect to event time. For example, if `balance_e=2`, `agg_ddd` will
  drop groups that are not exposed to treatment for at least three
  periods, the initial period `e=0` as well as the next two periods,
  `e=1` and `e=2`. This ensures that the composition of groups does not
  change when event time changes.

- min_e:

  For event studies, this is the smallest event time to compute dynamic
  effects for. By default, `min_e = -Inf` so that effects at all lengths
  of exposure are computed.

- max_e:

  For event studies, this is the largest event time to compute dynamic
  effects for. By default, `max_e = Inf` so that effects at all lengths
  of exposure are computed.

- na.rm:

  Logical value if we are to remove missing Values from analyses.
  Defaults is FALSE.

- boot:

  Boolean for whether or not to compute standard errors using the
  multiplier bootstrap. If standard errors are clustered, then one must
  set `boot=TRUE`. Default is value set in the ddd object. If
  `boot = FALSE`, then analytical standard errors are reported.

- nboot:

  The number of bootstrap iterations to use. The default is the value
  set in the ddd object, and this is only applicable if `boot=TRUE`.

- cband:

  Boolean for whether or not to compute a uniform confidence band that
  covers all of the group-time average treatment effects with fixed
  probability `0.95`. In order to compute uniform confidence bands,
  `boot` must also be set to `TRUE`. The default is the value set in the
  ddd object

- alpha:

  The level of confidence for the confidence intervals. The default is
  0.05. Otherwise, it will use the value set in the ddd object.

## Value

A object (list) of class `agg_ddd` that holds the results from the
aggregation step.

## Examples

``` r
#----------------------------------------------------------
# Triple Diff with multiple time periods
#----------------------------------------------------------

data <- gen_dgp_mult_periods(size = 500, dgp_type = 1)[["data"]]

out <- ddd(yname = "y", tname = "time", idname = "id",
            gname = "state", pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
            data = data, control_group = "nevertreated", base_period = "varying",
            est_method = "dr")
# Simple aggregation
agg_ddd(out, type = "simple", alpha = 0.10)
#>  Call:
#> agg_ddd(ddd_obj = out, type = "simple", alpha = 0.1)
#> ========================= DDD Aggregation ============================
#>  Overall ATT:
#>     ATT    Std. Error     [ 90%  Conf. Int.]  
#>  19.005        0.3718    18.3935     19.6166 *
#> 
#>  Note: * indicates that the confidence interval does not contain zero.
#>  --------------------------- Data Info   -----------------------------
#>  Outcome variable: y
#>  Qualification variable: partition
#>  Control group:  Never Treated
#>  --------------------------- Std. Errors  ----------------------------
#>  Level of significance:  0.1
#>  Analytical standard errors.
#>  =====================================================================
#>  See Ortiz-Villavicencio and Sant'Anna (2025) for details.

# Event study aggregation
agg_ddd(out, type = "eventstudy", alpha = 0.10)
#>  Call:
#> agg_ddd(ddd_obj = out, type = "eventstudy", alpha = 0.1)
#> ========================= DDD Aggregation ============================
#>  Overall summary of ATT's based on event-study/dynamic aggregation: 
#>      ATT    Std. Error     [ 90%  Conf. Int.]  
#>  19.2985        0.3653    18.6976     19.8993 *
#> 
#>  Event Study:
#>  Event time Estimate Std. Error [90% Ptwise.  Conf. Band]  
#>          -1  -0.2292     0.4289       -0.9346      0.4762  
#>           0  18.5301     0.4589       17.7752     19.2850 *
#>           1  20.0668     0.4283       19.3623     20.7714 *
#> 
#>  Note: * indicates that the confidence interval does not contain zero.
#>  --------------------------- Data Info   -----------------------------
#>  Outcome variable: y
#>  Qualification variable: partition
#>  Control group:  Never Treated
#>  --------------------------- Std. Errors  ----------------------------
#>  Level of significance:  0.1
#>  Analytical standard errors.
#>  =====================================================================
#>  See Ortiz-Villavicencio and Sant'Anna (2025) for details.

# Group aggregation
agg_ddd(out, type = "group", alpha = 0.10)
#>  Call:
#> agg_ddd(ddd_obj = out, type = "group", alpha = 0.1)
#> ========================= DDD Aggregation ============================
#>  Overall summary of ATT's based on group/cohort aggregation: 
#>      ATT    Std. Error     [ 90%  Conf. Int.]  
#>  20.7971        0.2922    20.3166     21.2777 *
#> 
#>  Group Effects:
#>  Group Estimate Std. Error [90% Ptwise.  Conf. Band]  
#>      2  14.9986     0.3623       14.4026     15.5946 *
#>      3  25.4900     0.4594       24.7344     26.2456 *
#> 
#>  Note: * indicates that the confidence interval does not contain zero.
#>  --------------------------- Data Info   -----------------------------
#>  Outcome variable: y
#>  Qualification variable: partition
#>  Control group:  Never Treated
#>  --------------------------- Std. Errors  ----------------------------
#>  Level of significance:  0.1
#>  Analytical standard errors.
#>  =====================================================================
#>  See Ortiz-Villavicencio and Sant'Anna (2025) for details.

# Calendar aggregation
agg_ddd(out, type = "calendar", alpha = 0.10)
#>  Call:
#> agg_ddd(ddd_obj = out, type = "calendar", alpha = 0.1)
#> ========================= DDD Aggregation ============================
#>  Overall summary of ATT's based on calendar time aggregation: 
#>      ATT    Std. Error     [ 90%  Conf. Int.]  
#>  16.4972        0.2783    16.0395      16.955 *
#> 
#>  Calendar Effects:
#>  Time Estimate Std. Error [90% Ptwise.  Conf. Band]  
#>     2   9.9303     0.4422        9.2029     10.6577 *
#>     3  23.0642     0.3908       22.4214     23.7069 *
#> 
#>  Note: * indicates that the confidence interval does not contain zero.
#>  --------------------------- Data Info   -----------------------------
#>  Outcome variable: y
#>  Qualification variable: partition
#>  Control group:  Never Treated
#>  --------------------------- Std. Errors  ----------------------------
#>  Level of significance:  0.1
#>  Analytical standard errors.
#>  =====================================================================
#>  See Ortiz-Villavicencio and Sant'Anna (2025) for details.

```
