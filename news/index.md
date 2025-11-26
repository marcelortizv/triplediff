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
