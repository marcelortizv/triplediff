## Resubmission

This is a resubmission of 0.2.x. The previous upload (0.2.3) was flagged on
`r-devel-linux-x86_64-debian-gcc` with the NOTE "Examples with CPU time > 2.5
times elapsed time" for the `ddd` example, caused by a multi-threaded BLAS
parallelising the multiplier-bootstrap matrix algebra during the example. In
this version (0.2.4) the bootstrap-based example has been removed so the
examples no longer trigger multi-core execution. No package functionality
changed.

This release also resolves the deprecation WARNING reported for 0.2.0 (the
`examples` check emitted `'BMisc::rhs.vars' is deprecated`): all internal calls
now use the current `BMisc` API (`BMisc (>= 1.4.9)`, available on CRAN).

## Test environments

* local macOS, R 4.x: R CMD check --as-cran
* win-builder: R-devel and R-release
* GitHub Actions: Ubuntu (R-release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies for this package.
