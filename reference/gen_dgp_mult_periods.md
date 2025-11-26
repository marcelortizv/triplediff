# Generate panel data with staggered treatment adoption (three periods)

Generate panel data where units adopt treatment at different times
across three periods.

## Usage

``` r
gen_dgp_mult_periods(size, dgp_type = 1)
```

## Arguments

- size:

  Integer. Number of units to simulate.

- dgp_type:

  Integer in {1,2,3,4}. 1 = both nuisance functions correct; 2 = only
  the outcome model correct; 3 = only the propensity-score model
  correct; 4 = both nuisance functions misspecified.

## Value

A named list with components:

- data:

  A `data.table` in long format with columns:

  - `id`: unit identifier

  - `cohort`: first period when treatment is assigned

  - `partition`: partition indicator

  - `x1`, `x2`, `x3`, `x4`: covariates

  - `cluster`: cluster identifier (no within-cluster correlation)

  - `time`: time period index

  - `y`: observed outcome

- data_wide:

  A `data.table` in wide format (one row per `id`) with columns:

  - `id`, `cohort`, `partition`, `x1`, `x2`, `x3`, `x4`, `cluster`

  - `y_t0`, `y_t1`, `y_t2`: outcomes in periods 0, 1, and 2

- ES_0_unf:

  Unfeasible (oracle) event-study parameter at time 0.

- prob_g2_p1:

  Proportion of units with `cohort == 2` and eligibility in period 1.

- prob_g3_p1:

  Proportion of units with `cohort == 3` and eligibility in period 1.
