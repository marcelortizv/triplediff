# Multiplier Bootstrap

This function take an influence function and use the multiplier
bootstrap to compute standard errors and critical values for uniform
confidence bands.

## Usage

``` r
mboot(inf_func, did_preprocessed, use_parallel = FALSE, cores = 1)
```

## Arguments

- inf_func:

  an influence function

- did_preprocessed:

  A `dp` object obtained after preprocess

- use_parallel:

  Boolean of whether or not to use parallel processing in the multiplier
  bootstrap, default is `use_parallel=FALSE`

- cores:

  the number of cores to use with parallel processing, default is
  `cores=1`

## Value

A list with following elements

- bres:

  results from each bootstrap iteration.

- V:

  variance matrix.

- se:

  standard errors.

- crit_val:

  a critical value for computing uniform confidence bands.
