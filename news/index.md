# Changelog

## triplediff 0.1.0

CRAN release: 2025-07-23

- Initial release of triplediff in alpha stage, functions for computing
  group-time average treatment effects in DDD and combining them into a
  smaller number of parameters are available.

## triplediff 0.1.1

- Bug fix in `cluster` parameter. When user specifies a cluster
  variable, the function now correctly uses it for clustering standard
  errors performing Multiplier Bootstrap.

## triplediff 0.1.2

- Bug fix in `preprocess` when checking for small groups.

## triplediff 0.2.0

CRAN release: 2026-01-16

- Replaced `parglm` with `fastglm` to avoid issues related to parglm’s
  scheduled archival on 2026-01-29.
- Added support for unbalanced panel data and repeated cross-sectional
  data by properly implementing the `allow_unbalanced_panel` parameter
  across all functions.

## triplediff 0.2.1

- Add asymmetric propensity score trimming for control units with pscore
  \>= 0.995.
- Add partition-specific collinearity detection with two-stage checking.
- Add comprehensive test suite including Monte Carlo coverage test when
  trimming.

## triplediff 0.2.2

- Track BMisc (\>= 1.4.9) API rename: replaced `makeBalancedPanel` with
  `make_balanced_panel` and `rhs.vars` with `rhs_vars` in internal
  preprocessing. No user-visible behavior change. (#34)
