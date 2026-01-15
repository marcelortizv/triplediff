# Doubly Robust DDD estimators for the group-time average treatment effects.

`ddd` is the main function for computing the Doubly Robust DDD
estimators for the ATT, with balanced panel data. It can be used with
covariates and/or under multiple time periods. At its core, `triplediff`
employs the doubly robust estimator for the ATT, which is a combination
of the propensity score weighting and the outcome regression.
Furthermore, this package supports the application of machine learning
methods for the estimation of the nuisance parameters.

## Usage

``` r
ddd(
  yname,
  tname,
  idname = NULL,
  gname,
  pname,
  xformla,
  data,
  control_group = NULL,
  base_period = NULL,
  est_method = "dr",
  panel = TRUE,
  allow_unbalanced_panel = FALSE,
  weightsname = NULL,
  boot = FALSE,
  nboot = NULL,
  cluster = NULL,
  cband = FALSE,
  alpha = 0.05,
  use_parallel = FALSE,
  cores = 1,
  inffunc = FALSE,
  skip_data_checks = FALSE
)
```

## Arguments

- yname:

  The name of the outcome variable.

- tname:

  The name of the column containing the time periods.

- idname:

  The name of the column containing the unit id.

- gname:

  The name of the column containing the first period when a particular
  observation is treated. It is a positive number for treated units and
  defines which group the unit belongs to. It takes value 0 or Inf for
  untreated units.

- pname:

  The name of the column containing the partition variable (e.g., the
  subgroup identifier). This is an indicator variable that is 1 for the
  units eligible for treatment and 0 otherwise.

- xformla:

  The formula for the covariates to be included in the model. It should
  be of the form `~ x1 + x2`. Default is `xformla = ~1` (no covariates).

- data:

  A data frame or data table containing the data.

- control_group:

  Valid for multiple periods only. The control group to be used in the
  estimation. Default is `control_group = "notyettreated"` which sets as
  control group the units that have not yet participated in the
  treatment. The alternative is `control_group = "nevertreated"` which
  sets as control group the units that never participate in the
  treatment and does not change across groups or time periods.

- base_period:

  Valid for multiple periods. Choose between a "varying" or "universal"
  base period. Both yield the same post-treatment ATT(g,t) estimates.
  Varying base period: Computes pseudo-ATT in pre-treatment periods by
  comparing outcome changes for a group to its comparison group from t-1
  to t, repeatedly changing t. Universal base period: Fixes the base
  period to (g-1), reporting average changes from t to (g-1) for a group
  relative to its comparison group, similar to event study regressions.
  Varying base period reports ATT(g,t) right before treatment. Universal
  base period normalizes the estimate before treatment to be 0, adding
  one extra estimate in an earlier period.

- est_method:

  The estimation method to be used. Default is `"dr"` (doubly robust).
  It computes propensity score using logistic regression and outcome
  regression using OLS. The alternative are `c("reg", "ipw")`.

- panel:

  Logical. If `TRUE` (default), the data is treated as panel data where
  each unit is observed in all time periods. If `FALSE`, the data is
  treated as repeated cross-sections (RCS) where each observation may
  represent a different unit. For RCS data, `idname` can be omitted or
  set to `NULL`, and the function will automatically create unique IDs
  for each observation.

- allow_unbalanced_panel:

  Logical. If `TRUE`, allows for unbalanced panel data where units may
  not be observed in all time periods. Default is `FALSE`. Note: This
  parameter requires `panel = TRUE` and a valid `idname`.

- weightsname:

  The name of the column containing the weights. Default is `NULL`. As
  part of data processing, weights are enforced to be normalized and
  have mean 1 across all observations.

- boot:

  Logical. If `TRUE`, the function computes standard errors using the
  multiplier bootstrap. Default is `FALSE`.

- nboot:

  The number of bootstrap samples to be used. Default is `NULL`. If
  `boot = TRUE`, the default is `nboot = 999`.

- cluster:

  The name of the variable to be used for clustering. The maximum number
  of cluster variables is 1. Default is `NULL`. If `boot = TRUE`, the
  function computes the bootstrap standard errors clustering at the unit
  level setting as cluster variable the one in `idname`.

- cband:

  Logical. If `TRUE`, the function computes a uniform confidence band
  that covers all of the average treatment effects with fixed
  probability `1-alpha`. In order to compute uniform confidence bands,
  `boot` must also be set to `TRUE`. The default is `FALSE`.

- alpha:

  The level of significance for the confidence intervals. Default is
  `0.05`.

- use_parallel:

  Logical. If `TRUE`, the function runs in parallel processing. Valid
  only when `boot = TRUE`. Default is `FALSE`.

- cores:

  The number of cores to be used in the parallel processing. Default is
  `cores = 1`.

- inffunc:

  Logical. If `TRUE`, the function returns the influence function.
  Default is `FALSE`.

- skip_data_checks:

  Logical. If `TRUE`, the function skips data validation checks and
  proceeds directly to estimation. This can improve performance when you
  are confident the data is correctly formatted. Default is `FALSE`. Use
  with caution as skipping checks may lead to unexpected errors if data
  is malformed.

## Value

A `ddd` object with the following basic elements:

- ATT:

  The average treatment effect on the treated.

- se:

  The standard error of the ATT.

- uci:

  The upper confidence interval of the ATT.

- lci:

  The lower confidence interval of the ATT.

- inf_func:

  The estimate of the influence function.

## Examples

``` r
#----------------------------------------------------------
# Triple Diff with covariates and 2 time periods
#----------------------------------------------------------
set.seed(1234) # Set seed for reproducibility
# Simulate data for a two-periods DDD setup
df <- gen_dgp_2periods(size = 5000, dgp_type = 1)$data

head(df)
#> Key: <id>
#>       id state partition  time        y        cov1       cov2      cov3
#>    <int> <num>     <num> <int>    <num>       <num>      <num>     <num>
#> 1:     1     0         0     1 209.9152 -0.97080934 -1.1726958 2.3893945
#> 2:     1     0         0     2 417.5260 -0.97080934 -1.1726958 2.3893945
#> 3:     2     0         0     1 211.4919  0.02591115  0.2763066 0.1063123
#> 4:     2     0         0     2 420.3656  0.02591115  0.2763066 0.1063123
#> 5:     3     0         0     1 221.9431  0.97147321 -0.4292088 0.5012794
#> 6:     3     0         0     2 440.9623  0.97147321 -0.4292088 0.5012794
#>          cov4 cluster
#>         <num>   <int>
#> 1:  0.2174955      39
#> 2:  0.2174955      39
#> 3: -0.1922253      29
#> 4: -0.1922253      29
#> 5:  1.1027248      44
#> 6:  1.1027248      44

att_22 <- ddd(yname = "y", tname = "time", idname = "id", gname = "state",
              pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
             data = df, control_group = "nevertreated", est_method = "dr")

summary(att_22)
#>  Call:
#> ddd(yname = "y", tname = "time", idname = "id", gname = "state", 
#>     pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4, 
#>     data = df, control_group = "nevertreated", est_method = "dr")
#> =========================== DDD Summary ==============================
#>  DR-DDD estimation for the ATT: 
#>      ATT       Std. Error    Pr(>|t|)  [95% Ptwise. Conf. Band]              
#>     -0.0780       0.0828       0.3463      -0.2404       0.0843              
#> 
#>  Note: * indicates that the confidence interval does not contain zero.
#>  --------------------------- Data Info   -----------------------------
#>  Panel Data
#>  Outcome variable: y
#>  Qualification variable: partition
#>  No. of units at each subgroup:
#>    treated-and-eligible: 1232
#>    treated-but-ineligible: 1285
#>    eligible-but-untreated: 1256
#>    untreated-and-ineligible: 1227
#>  --------------------------- Algorithms ------------------------------
#>  Outcome Regression estimated using: OLS
#>  Propensity score estimated using: Maximum Likelihood
#>  --------------------------- Std. Errors  ----------------------------
#>  Level of significance:  0.05
#>  Analytical standard errors.
#>  Type of confidence band:  Pointwise Confidence Interval
#>  =====================================================================
#>  See Ortiz-Villavicencio and Sant'Anna (2025) for details.

# Performing clustered standard errors with mutiplier bootstrap

att_cluster <-  ddd(yname = "y", tname = "time", idname = "id", gname = "state",
                    pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
                    data = df, control_group = "nevertreated",
                    base_period = "universal", est_method = "dr", 
                    boot = TRUE, nboot = 500, cband = TRUE, cluster = "cluster")

summary(att_cluster)
#>  Call:
#> ddd(yname = "y", tname = "time", idname = "id", gname = "state", 
#>     pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4, 
#>     data = df, control_group = "nevertreated", base_period = "universal", 
#>     est_method = "dr", boot = TRUE, nboot = 500, cluster = "cluster", 
#>     cband = TRUE)
#> =========================== DDD Summary ==============================
#>  DR-DDD estimation for the ATT: 
#>      ATT       Std. Error    Pr(>|t|)  [95% Simult. Conf. Band]              
#>     -0.0780       0.0944       0.4088      -0.2647       0.1087              
#> 
#>  Note: * indicates that the confidence interval does not contain zero.
#>  --------------------------- Data Info   -----------------------------
#>  Panel Data
#>  Outcome variable: y
#>  Qualification variable: partition
#>  No. of units at each subgroup:
#>    treated-and-eligible: 1232
#>    treated-but-ineligible: 1285
#>    eligible-but-untreated: 1256
#>    untreated-and-ineligible: 1227
#>  --------------------------- Algorithms ------------------------------
#>  Outcome Regression estimated using: OLS
#>  Propensity score estimated using: Maximum Likelihood
#>  --------------------------- Std. Errors  ----------------------------
#>  Level of significance:  0.05
#>  Boostrapped standard error based on 500 reps. 
#>  Method: Multiplier Bootstrap.
#>  Type of confidence band:  Uniform Confidence Band 
#>  Clustering Std. Errors by: cluster
#>  =====================================================================
#>  See Ortiz-Villavicencio and Sant'Anna (2025) for details.

#----------------------------------------------------------
# Triple Diff with multiple time periods
#----------------------------------------------------------
data <- gen_dgp_mult_periods(size = 1000, dgp_type = 1)[["data"]]

ddd(yname = "y", tname = "time", idname = "id",
     gname = "state", pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4,
     data = data, control_group = "nevertreated", base_period = "varying",
     est_method = "dr")
#>  Call:
#> ddd(yname = "y", tname = "time", idname = "id", gname = "state", 
#>     pname = "partition", xformla = ~cov1 + cov2 + cov3 + cov4, 
#>     data = data, control_group = "nevertreated", base_period = "varying", 
#>     est_method = "dr")
#> =========================== DDD Summary ==============================
#>  DR-DDD estimation for the ATT(g,t): 
#> Group Time  ATT(g,t)  Std. Error [95% Pointwise  Conf. Band]  
#>   2    2      10.1689     0.2879       9.6047       10.7331  *
#>   2    3      20.0451     0.2776      19.5010       20.5893  *
#>   3    2      -0.1446     0.2736      -0.6808        0.3916   
#>   3    3      25.2267     0.3106      24.6178       25.8355  *
#> 
#>  Note: * indicates that the confidence interval does not contain zero.
#>  --------------------------- Data Info   -----------------------------
#>  Panel Data
#>  Outcome variable: y
#>  Qualification variable: partition
#>  Control group: Never Treated
#>  No. of units per treatment group:
#>   Units enabling treatment at period 3: 443
#>   Units enabling treatment at period 2: 356
#>   Units never enabling treatment: 201
#>  --------------------------- Algorithms ------------------------------
#>  Outcome Regression estimated using: OLS
#>  Propensity score estimated using: Maximum Likelihood
#>  --------------------------- Std. Errors  ----------------------------
#>  Level of significance:  0.05
#>  Analytical standard errors.
#>  Type of confidence band:  Pointwise Confidence Interval
#>  =====================================================================
#>  See Ortiz-Villavicencio and Sant'Anna (2025) for details.
```
