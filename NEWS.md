# triplediff 0.1.0

  * Initial release of triplediff in alpha stage, functions for computing group-time average treatment effects in DDD and combining them into a smaller number of parameters are available.

# triplediff 0.1.1

  * Bug fix in `cluster` parameter. When user specifies a cluster variable, the function now correctly uses it for clustering standard errors performing Multiplier Bootstrap.

# triplediff 0.1.2

  * Bug fix in `preprocess` when checking for small groups.

# triplediff 0.2.0

  * Replaced `parglm` with `fastglm` to avoid issues related to parglm's scheduled archival on 2026-01-29.
  * Added support for unbalanced panel data and repeated cross-sectional data by properly implementing the `allow_unbalanced_panel` parameter across all functions.

# triplediff 0.2.1
  * Add asymmetric propensity score trimming for control units with pscore >= 0.995. 
  * Add partition-specific collinearity detection with two-stage checking.
  * Add comprehensive test suite including Monte Carlo coverage test when trimming. 


# triplediff 0.2.2

  * Track BMisc (>= 1.4.9) API rename: replaced `makeBalancedPanel` with `make_balanced_panel` and `rhs.vars` with `rhs_vars` in internal preprocessing. No user-visible behavior change. (#34)
  * Replaced the remaining deprecated `BMisc::getListElement` call with `BMisc::get_list_element` to silence deprecation warnings (follow-up to #34).

# triplediff 0.2.3

  * Added analytical cluster-robust standard errors without the bootstrap in the multiple-period path. Calling `ddd()` with `cluster = <var>` and `boot = FALSE` now returns analytical cluster-robust standard errors (cluster-sum CRVE on the influence function) instead of requiring the bootstrap. The `ddd` object now carries `cluster_vector` and `cluster_var`. Two-period designs still require `boot = TRUE` for clustered inference.
  * Added a `cluster` argument to `agg_ddd()`. Aggregated parameters (simple, event-study, group, and calendar) now report analytical cluster-robust standard errors. If clustering is requested on a different variable than `ddd()` used (or on an object built without clustering), `agg_ddd()` warns and falls back to i.i.d. standard errors instead of silently mis-reporting.
  * Behavior change: the clustered multiplier bootstrap now follows Callaway & Sant'Anna (2021, Remark 10), applying one multiplier per cluster to the influence function aggregated to cluster *sums* rather than cluster *means*. Clustered bootstrap standard errors change for unbalanced clusters and repeated cross-sections; equal-sized clusters are unaffected.

# triplediff 0.2.4

  * Resubmission. Simplified the `ddd()` examples by removing the bootstrap-based clustered standard errors example, so the package examples stay within CRAN's limit on the number of cores used during checks. No changes to package functionality.

# triplediff 0.3.0 (development version)

  * Added the `gruber1994` dataset: the Gruber (1994) maternity-mandate data. It holds repeated cross-sections of workers from the May Current Population Survey for 1974, 1975, 1977 and 1978 in 8 states (27,575 rows, 29 columns). In `ddd()`, `treatment_starts` plays the role of `gname` and `treat_group` the role of `pname`; use `panel = FALSE`. See `?gruber1994`.
  * Added `plot()` methods for `ddd` and `agg_ddd` objects. They return a `ggplot` object and do not print it. On a `ddd()` fit, `plot()` draws a coefficient plot with one row per ATT(g,t). Each estimate is a point with its confidence interval, and a dashed line marks zero. The zero line is always in view. On an `agg_ddd()` result, `plot()` draws the event-study, group or calendar estimates with their intervals. In an event study, pre- and post-treatment estimates have different colors. The plots draw the intervals that are stored in the object. If `cband = TRUE` was used, this is the uniform confidence band. Otherwise, it is pointwise. ggplot2 is a suggested package, not a required one. Install it to use these plots.

## Breaking changes

  * `agg_ddd()` now weights each cohort by its share of treated units, P(S = g, Q = 1), in the simple, event-study, group and calendar aggregations. A treated unit is a unit of cohort g that is also eligible (`pname` equal to 1). These are the weights of the event-study parameter in Ortiz-Villavicencio and Sant'Anna (2025). Before, cohorts were weighted by P(S = g), which also counted the ineligible units of each cohort. Aggregates that average more than one cohort, and their standard errors, change. Estimates change when eligibility shares differ across cohorts. The ATT(g,t) estimates and the aggregates of a single cohort do not change. See "Cohort weights" in `?agg_ddd`.
  * With multiple periods and `base_period = "universal"`, pre-treatment ATT(g,t) with t < g - 1 now use the covariates of period t in both nuisance models. These models are the outcome regression and the propensity score. Period t is the earlier of the two periods, and the did package uses the same rule. Before, the outcome regression used the covariates of period g - 1, and the propensity score used those of period t. The `"dr"` and `"reg"` estimates of these cells and their standard errors change when covariates vary over time. Event-study estimates for event times -2 and earlier change with them. This also makes the influence function of these cells exact, so their standard errors now match the estimator. All other cells, `est_method = "ipw"`, `base_period = "varying"`, time-invariant covariates, repeated cross-sections and two-period designs do not change.

## Bug fixes

  * With `control_group = "notyettreated"` and a balanced panel, some standard errors were too small. This affected cells that combine several comparison cohorts (GMM) and leave out at least one treated cohort. The reported standard error was about sqrt(n_cell / n) times the correct one, where n_cell is the number of units in the cell and n is the number of units. The influence function of each comparison cohort went to the wrong units. It now goes to the units of that cohort. The GMM weights change too, so the point estimates of these cells also move slightly, and so do the aggregates that use them. Other cells, never-treated fits and repeated cross-sections do not change.
  * The overall standard error of `agg_ddd(type = "group")` was wrong when cohort effects differ. The weight term of every cohort used the indicator of the first cohort. Each cohort now uses its own indicator. The estimates and the other aggregation types do not change.
  * `agg_ddd()` now uses the `alpha` that is stored in the `ddd` object, as its help page said. The default is now `alpha = NULL`. Before, the default was 0.05 for every `ddd` object. Results change only when `ddd()` was called with an `alpha` other than 0.05. If you give `alpha` to `agg_ddd()`, it uses that value.
  * Two-period panels: `ddd()` now sorts the rows by unit and period and balances the panel before it checks that the covariates are time invariant. Before, a balanced panel with rows out of unit order stopped with a false "Covariates are not time invariant" error. An unbalanced panel with covariates stopped with "non-conformable arrays". With `allow_unbalanced_panel = TRUE`, the panel estimator sometimes paired the wrong rows and returned a wrong estimate without an error. By default, `ddd()` now drops the units that are not observed in both periods, with a warning. With `allow_unbalanced_panel = TRUE`, it now uses the repeated cross-section estimator. Its standard errors sum the influence function by unit, so a unit observed in both periods counts once. Balanced panels in sorted order give the same results as before.
  * `ddd()` checks more of its input. It stops when a unit has more than one row in the same period. Before, this check did not look at the data. `pname` must now be coded 0 and 1. Before, other codings were misread, and `est_method = "ipw"` returned NaN without an error. With two periods, `gname = Inf` now marks never-treated units, as with multiple periods, and any other never-treated code stops with an error. Before, a coding such as {Inf, g} or two nonzero values silently flipped the sign of the estimate.
  * `cluster = idname` now works and is the same as `cluster = NULL`, because standard errors are always clustered at the unit level. Before, it stopped with "argument is of length zero". `cluster = c(idname, "x")` is now the same as `cluster = "x"`, and `agg_ddd()` then clusters on `x`. With multiple periods, sampling weights must now be the same in every period for each unit, as with two periods. Before, each ATT(g,t) silently used the weights of its base period.
  * With multiple periods, repeated cross-sections without an id (`panel = FALSE`, `idname = NULL`) now work with `cluster`. Before, every such call stopped with "object '.row_id' not found". The results are the same as for a fit with a unique row id.
  * With no never-treated group, the last cohort to be treated serves only as a comparison group. If some units were treated in the first period, `ddd()` dropped those units and also that comparison cohort. The result was a malformed object (for example, one ATT next to four standard errors), or an error with `base_period = "universal"`. Now `ddd()` drops only the units treated in the first period. Also, `ddd(cband = TRUE, boot = FALSE)` now warns, sets `cband = FALSE` and reports pointwise intervals. Before, it gave no warning, and every later `agg_ddd()` call with default arguments stopped with "cband is only available when boot = TRUE".
  * With two periods, `panel = FALSE` and `idname = NULL`, `ddd()` no longer adds a `.row_id` column to a data.table passed as `data`. Before, a second call on the same table stopped with an error.
  * With two periods, `panel = FALSE`, `idname = NULL` and `skip_data_checks = TRUE`, `ddd()` now runs and gives the same result as with the data checks. Before, it stopped with an error about a missing `.row_id` column.
