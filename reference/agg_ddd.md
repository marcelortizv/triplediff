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
  cluster = NULL,
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

- cluster:

  The name of the variable to be used for clustering. The maximum number
  of cluster variables is 1. Default is `NULL`.

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
#>      ATT    Std. Error     [ 90%  Conf. Int.]  
#>  18.6383        0.3347    18.0878     19.1888 *
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
#>  18.9605        0.3254    18.4253     19.4956 *
#> 
#>  Event Study:
#>  Event time Estimate Std. Error [90% Ptwise.  Conf. Band]  
#>          -1  -0.1004     0.4030       -0.7633      0.5625  
#>           0  18.0946     0.4280       17.3907     18.7986 *
#>           1  19.8263     0.3742       19.2108     20.4417 *
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
#>  20.2332        0.2551    19.8135     20.6528 *
#> 
#>  Group Effects:
#>  Group Estimate Std. Error [90% Ptwise.  Conf. Band]  
#>      2  15.1527     0.3421        14.590     15.7155 *
#>      3  24.5191     0.3921        23.874     25.1641 *
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
#>  16.4255        0.2721    15.9779      16.873 *
#> 
#>  Calendar Effects:
#>  Time Estimate Std. Error [90% Ptwise.  Conf. Band]  
#>     2  10.4792     0.4294        9.7729     11.1855 *
#>     3  22.3717     0.3237       21.8393     22.9041 *
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
