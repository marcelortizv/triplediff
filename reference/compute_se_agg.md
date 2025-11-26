# Take influence function and compute standard errors

Function to take an nx1 influence function and return a standard error

## Usage

``` r
compute_se_agg(influence_function, boot = FALSE, boot_std_errors = NA)
```

## Arguments

- influence_function:

  An influence function

- boot:

  a boolean indicating whether bootstrapping was performed

- boot_std_errors:

  a vector of bootstrapped standard errors

## Value

scalar standard error
