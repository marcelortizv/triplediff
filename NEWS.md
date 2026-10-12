# triplediff 0.1.0

- Initial release of triplediff in alpha stage, functions for computing group-time average treatment effects in DDD and combining them into a smaller number of parameters are available.

# triplediff 0.1.1

- Bug fix in `cluster` parameter. When user specifies a cluster variable, the function now correctly uses it for clustering standard errors performing Multiplier Bootstrap.

# triplediff 0.1.2

- Bug fix in `preprocess` when checking for small groups.

# triplediff 0.2.0

- Replaced `parglm` with `fastglm` to avoid issues related to parglm's scheduled archival on 2026-01-29.
- Added support for unbalanced panel data and repeated cross-sectional data by properly implementing the `allow_unbalanced_panel` parameter across all functions.

# triplediff 0.2.1

- Add asymmetric propensity score trimming for control units with pscore &gt;= 0.995. 
- Add partition-specific collinearity detection with two-stage checking.
- Add comprehensive test suite including Monte Carlo coverage test when trimming.

# triplediff 0.2.2

- Track BMisc (&gt;= 1.4.9) API rename: replaced `makeBalancedPanel` with `make_balanced_panel` and `rhs.vars` with `rhs_vars` in internal preprocessing. No user-visible behavior change. (#34)
- Replaced the remaining deprecated `BMisc::getListElement` call with `BMisc::get_list_element` to silence deprecation warnings (follow-up to #34).

# triplediff 0.2.3

- Added analytical cluster-robust standard errors without the bootstrap in the multiple-period path. Calling `ddd()` with `cluster = <var>` and `boot = FALSE` now returns analytical cluster-robust standard errors (cluster-sum CRVE on the influence function) instead of requiring the bootstrap. The `ddd` object now carries `cluster_vector` and `cluster_var`. Two-period designs still require `boot = TRUE` for clustered inference.
- Added a `cluster` argument to `agg_ddd()`. Aggregated parameters (simple, event-study, group, and calendar) now report analytical cluster-robust standard errors. If clustering is requested on a different variable than `ddd()` used (or on an object built without clustering), `agg_ddd()` warns and falls back to i.i.d. standard errors instead of silently mis-reporting.
- Behavior change: the clustered multiplier bootstrap now follows Callaway &amp; Sant'Anna (2021, Remark 10), applying one multiplier per cluster to the influence function aggregated to cluster *sums* rather than cluster *means*. Clustered bootstrap standard errors change for unbalanced clusters and repeated cross-sections; equal-sized clusters are unaffected.

# triplediff 0.2.4

- Resubmission. Simplified the `ddd()` examples by removing the bootstrap-based clustered standard errors example, so the package examples stay within CRAN's limit on the number of cores used during checks. No changes to package functionality.

# triplediff 0.3.0 

- Added the `gruber1994` dataset with the Gruber (1994) maternity-mandate data. See `?gruber1994`.
- Added `plot()` methods for `ddd` and `agg_ddd` objects that return a `ggplot` object (requires ggplot2). See `?plot.ddd` and `?plot.agg_ddd`.

## Bug fixes

- `agg_ddd()`  weights cohorts by P(S = g, Q = 1) in all aggregation types. Aggregates over more than one cohort change. See "Cohort weights" in `?agg_ddd`.
- With `base_period = "universal"`, pre-treatment ATT(g,t) with t &lt; g - 1 now use the covariates of period t in both nuisance models. 
- Fixed standard errors  for GMM cells with `control_group = "notyettreated"` in balanced panels. The point estimates of these cells and the aggregates that use them also move slightly.
- Fixed the overall standard error of `agg_ddd(type = "group")` when cohort effects differ.
- `agg_ddd()` now uses the `alpha` stored in the `ddd` object, and the default is `alpha = NULL`.
- Fixed `ddd()` with unsorted or unbalanced two-period panels, which gave false errors or wrong estimates. It now drops units not observed in both periods with a warning.
- `ddd()` now stops when a unit has more than one row in the same period.
- `ddd()` now requires `pname` to be coded 0 and 1, which fixes `est_method = "ipw"` returning NaN.
- With two periods, `ddd()` now requires `gname = Inf` for never-treated units and stops on any other never-treated code, which fixes a flipped sign of the estimate.
- `cluster = idname` now works and equals `cluster = NULL`. `cluster = c(idname, "x")` now equals `cluster = "x"`, and `agg_ddd()` then clusters on `x`.
- With multiple periods, `ddd()` now requires sampling weights to be the same in every period for each unit.
- Fixed `ddd()` with `cluster` for multiple-period repeated cross-sections without an id (`panel = FALSE`, `idname = NULL`).
- Fixed `ddd()` without a never-treated group when some units are treated in the first period. It now drops only those units and keeps the last-treated cohort as the comparison group.
- `ddd(cband = TRUE, boot = FALSE)` now warns and sets `cband = FALSE`, so later `agg_ddd()` calls no longer stop.
- With two periods, `panel = FALSE` and `idname = NULL`, `ddd()` no longer adds a `.row_id` column to a `data.table` passed as `data` and now runs with `skip_data_checks = TRUE`.