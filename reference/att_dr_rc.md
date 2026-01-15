# Doubly robust DDD estimator for ATT, with repeated cross-section data and 2 periods

This function implements a doubly robust estimator for assessing the
average treatment effect on the treated (ATT) using a triple differences
(DDD) approach in repeated cross-section data settings across two time
periods. The function takes preprocessed data structured specifically
for this analysis.

## Usage

``` r
att_dr_rc(did_preprocessed)
```

## Arguments

- did_preprocessed:

  A list containing preprocessed data and specifications for the DDD
  estimation. Expected elements include: - `preprocessed_data`: A data
  table containing the data with variables needed for the analysis. -
  `est_method`: The estimation method to be used. Default is
  `est_method = "dr"`. - `xformula`: The formula for the covariates to
  be included in the model. It should be of the form `~ x1 + x2`.
  Default is `xformla = ~1` (no covariates). - `boot`: Logical. If
  `TRUE`, the function use the multiplier bootstrap to compute standard
  errors. Default is `FALSE`. - `nboot`: The number of bootstrap samples
  to be used. Default is `NULL`. If `boot = TRUE`, the default is
  `nboot = 999`. - `subgroup_counts`: A matrix containing the number of
  observations in each subgroup. - `alpha` The level of significance for
  the confidence intervals. Default is `0.05`. - `inffunc`: Logical. If
  `TRUE`, the function returns the influence function. Default is
  `FALSE`. - `use_parallel`: Boolean of whether or not to use parallel
  processing in the multiplier bootstrap, default is
  `use_parallel=FALSE` - `cores`: the number of cores to use with
  parallel processing, default is `cores=1` - `cband`: Boolean of
  whether or not to compute simultaneous confidence bands, default is
  `cband=FALSE`

## Value

A list with the estimated ATT, standard error, upper and lower
confidence intervals, and influence function.
