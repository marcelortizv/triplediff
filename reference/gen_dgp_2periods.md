# Function that generates panel data with single treatment date assignment and two time periods.

Generate panel data with a single treatment date and two periods

## Usage

``` r
gen_dgp_2periods(size, dgp_type)
```

## Arguments

- size:

  Integer. Number of units.

- dgp_type:

  Integer in {1,2,3,4}. 1 = both nuisance functions correct; 2 = only
  the outcome model correct; 3 = only the propensity score correct; 4 =
  both nuisance functions incorrect.

## Value

A list with the following elements:

- data:

  A `data.table` in long format with columns:

  - `id`: unit identifier

  - `state`: state variable

  - `time`: time variable

  - `partition`: partition assignment

  - `x1`, `x2`, `x3`, `x4`: covariates

  - `y`: outcome variable

  - `cluster`: cluster ID (no within-cluster correlation)

- att:

  True average treatment effect on the treated (ATT), set to 0.

- att.unf:

  Oracle ATT computed under the unfeasible specification.

- eff:

  Theoretical efficiency bound for the estimator.
