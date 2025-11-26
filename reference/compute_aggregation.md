# Compute Aggregated Treatment Effect Parameters

Does the heavy lifting on computing aggregated group-time average
treatment effects

## Usage

``` r
compute_aggregation(
  ddd_obj,
  type = "simple",
  cluster = NULL,
  balance_e = NULL,
  min_e = -Inf,
  max_e = Inf,
  na.rm = FALSE,
  boot = FALSE,
  nboot = NULL,
  cband = NULL,
  alpha = 0.05
)
```

## Arguments

- ddd_obj:

  a ddd object (i.e., the results of the
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
  probability `1 - alpha`. In order to compute uniform confidence bands,
  `boot` must also be set to `TRUE`. The default is the value set in the
  ddd object

- alpha:

  The level of confidence for the confidence intervals. The default is
  0.05. Otherwise, it will use the value set in the ddd object.

## Value

Aggregation object (list) of class
[`agg_ddd`](http://marcelortiz.com/triplediff/reference/agg_ddd.md)
