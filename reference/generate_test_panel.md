# Function to generate a fake dataset for testing purposes only.

Function to generate fake dataset to test internal procedures.

## Usage

``` r
generate_test_panel(
  seed = 123,
  num_ids = 100,
  time = 2,
  initial.year = 2019,
  treatment.year = 2020
)
```

## Arguments

- seed:

  Seed for reproducibility

- num_ids:

  Number of IDs

- time:

  Number of time periods

- initial.year:

  Initial year

- treatment.year:

  Treatment year

## Value

A data.table with the following columns:

- id: ID

- state: State variable

- year: Time variable

- partition: Partition variable

- x1: Covariate 1

- x2: Covariate 2

- treat: Treatment variable

- outcome: Outcome variable
