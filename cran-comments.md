## Submission

This release fixes the deprecation WARNING reported on the CRAN checks for
version 0.2.0 (the `examples` check emitted
`'BMisc::rhs.vars' is deprecated`). All internal calls now use the current
`BMisc` API (`BMisc (>= 1.4.9)`, available on CRAN), so no deprecated functions
are called. The previously used `Remotes:` field has been removed now that the
required `BMisc` version is on CRAN.

### Changes in this version (0.2.3):

1. Fixed BMisc deprecation warnings: replaced the deprecated `rhs.vars`,
   `makeBalancedPanel`, and `getListElement` calls with their current
   equivalents (`rhs_vars`, `make_balanced_panel`, `get_list_element`). This
   resolves the `examples` WARNING reported for 0.2.0.

2. New functionality: added analytical cluster-robust standard errors (without
   the bootstrap) in the multiple-period path; added a `cluster` argument to
   `agg_ddd()` for clustered aggregated standard errors; and updated the
   clustered multiplier bootstrap to follow Callaway & Sant'Anna (2021,
   Remark 10).

## Test environments

* local macOS, R 4.x: R CMD check --as-cran
* win-builder: R-devel and R-release
* GitHub Actions: Ubuntu (R-release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies for this package.
