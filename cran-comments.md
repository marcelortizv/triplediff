## Resubmission
This is a resubmission to address the scheduled archival of the `parglm` package
and to fix a bug discovered in the previous version.

### Changes in this version:

1. **Removed parglm dependency**: Replaced `parglm` with `fastglm` to avoid issues
   related to parglm's scheduled archival on 2026-01-29. This change ensures
   long-term stability and removes the strong reverse dependency on parglm.

2. **Enhanced functionality**: Added support for unbalanced panel data and
   repeated cross-sectional data by properly implementing the `allow_unbalanced_panel`
   parameter across all functions.

## Test environments
* local macOS 14, R 4.4.0
* Ubuntu 24.04 (GitHub Actions), R 4.3.3 and 4.4.0
* Windows Server 2025 R-devel (win-builder)
* Windows Server 2025 R-release (win-builder)
* Windows Server 2025 R-oldrelease (win-builder)

## R CMD check results
0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Reverse dependencies
There are no reverse dependencies for this package.
