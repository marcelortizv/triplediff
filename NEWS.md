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
  * Added analytical cluster-robust standard errors without the bootstrap in the multiple-period path. Calling `ddd()` with `cluster = <var>` and `boot = FALSE` now returns analytical cluster-robust standard errors (cluster-sum CRVE on the influence function) instead of requiring the bootstrap. The `ddd` object now carries `cluster_vector` and `cluster_var`. Two-period designs still require `boot = TRUE` for clustered inference.
  * Added a `cluster` argument to `agg_ddd()`. Aggregated parameters (simple, event-study, group, and calendar) now report analytical cluster-robust standard errors. If clustering is requested on a different variable than `ddd()` used (or on an object built without clustering), `agg_ddd()` warns and falls back to i.i.d. standard errors instead of silently mis-reporting.
  * Behavior change: the clustered multiplier bootstrap now follows Callaway & Sant'Anna (2021, Remark 10), applying one multiplier per cluster to the influence function aggregated to cluster *sums* rather than cluster *means*. Clustered bootstrap standard errors change for unbalanced clusters and repeated cross-sections; equal-sized clusters are unaffected.

