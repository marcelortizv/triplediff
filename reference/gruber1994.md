# Gruber (1994) maternity mandates data

Repeated cross-sections of workers from the May Current Population
Survey (CPS) for 1974, 1975, 1977 and 1978, used by Gruber (1994) to
study who pays for state laws that required health insurance to cover
maternity. It is a classic triple-differences (DDD) application: mandate
states versus comparison states, before versus after the laws, and
married women aged 20 to 40 versus a comparison group of workers.

## Usage

``` r
gruber1994
```

## Format

A data frame with 27,653 rows and 29 columns:

- id:

  Integer. Row identifier, unique for each row. Each row is a different
  worker.

- year:

  Integer. Year of the May CPS: 1974 or 1975 (before the mandates) and
  1977 or 1978 (after).

- state:

  Character. State: IL, NJ or NY (mandate states) or OH, IN, CT, MA or
  NC (comparison states).

- treatment_starts:

  Integer. First treated period: 1976 for IL, NJ and NY, whose mandates
  took effect between July 1976 and January 1977; 0 for the comparison
  states (never treated).

- log_wage:

  Numeric. Log of `hourly_wage`.

- log_hours:

  Numeric. Log of `hours_per_week`.

- hourly_wage:

  Numeric. Hourly wage in 1978 dollars (see Details).

- hours_per_week:

  Numeric. Hours worked per week.

- education:

  Numeric. Years of education completed: the highest grade attended
  minus 1, minus 1 more if that grade was not completed, and 0 if
  negative.

- experience:

  Numeric. Potential experience in years: age minus education minus 6,
  set to 0 if negative.

- nonwhite:

  Integer. 1 if the worker is not white, 0 otherwise.

- union:

  Integer. 1 if the worker is a union member, 0 otherwise.

- treat:

  Integer. Same as `treat_group`.

- exp:

  Integer. Same as `exp_state`.

- post:

  Integer. Same as `post_period`.

- treat_group:

  Integer. 1 if the worker is a married woman aged 20 to 40 (spouse
  present or armed-forces spouse present), the pregnancy-relevant group;
  0 for the comparison group (workers older than 40 and never-married
  men aged 20 to 40).

- exp_state:

  Integer. 1 if the state passed a maternity mandate (IL, NJ, NY), 0
  otherwise.

- post_period:

  Integer. 1 for 1977 and 1978 (after the mandates), 0 for 1974 and
  1975.

- industry_3:

  Integer. 1 if the worker is in construction, 0 otherwise.

- industry_4:

  Integer. 1 if the worker is in manufacturing, 0 otherwise.

- industry_5:

  Integer. 1 if the worker is in transportation, 0 otherwise.

- industry_6:

  Integer. 1 if the worker is in wholesale trade, 0 otherwise.

- industry_7:

  Integer. 1 if the worker is in retail trade, 0 otherwise.

- industry_8:

  Integer. 1 if the worker is in finance, 0 otherwise.

- industry_9:

  Integer. 1 if the worker is in business services, 0 otherwise.

- industry_10:

  Integer. 1 if the worker is in personal services, 0 otherwise.

- industry_11:

  Integer. 1 if the worker is in entertainment, 0 otherwise.

- industry_12:

  Integer. 1 if the worker is in professional services, 0 otherwise.

- industry_13:

  Integer. 1 if the worker is in public administration, 0 otherwise.

## Source

NBER extracts of the May Current Population Survey for 1974, 1975, 1977
and 1978. Prepared by the package authors with the sample restrictions
of Gruber (1994).

## Details

Illinois, New Jersey and New York passed laws that required health
insurance to cover maternity like other conditions. The laws took effect
between July 1976 and January 1977. Gruber (1994) compares the change in
wages of married women aged 20 to 40 with the change for a comparison
group of workers, in the mandate states and in comparison states (OH and
IN for IL; CT, MA and NC for NJ and NY).

The sample follows Gruber (1994, footnote 15): workers aged 20 to 65,
not self-employed, with an hourly wage of at least \$1 and at most \$100
in 1978 dollars, and with non-missing education, experience, race, union
status and industry; only the treatment and comparison groups; only the
eight states above.

For workers paid by the hour, `hourly_wage` is the reported hourly rate;
for other workers it is weekly earnings divided by weekly hours. Nominal
wages are converted to 1978 dollars with the consumer price index (49.3
in 1974, 53.8 in 1975, 60.6 in 1977 and 65.2 in 1978).

The industry columns group 1970 Census industry codes into major
industries. Agriculture and mining are the omitted groups: a worker with
all industry columns equal to 0 is in one of these two industries.

## References

Gruber, J. (1994). The Incidence of Mandated Maternity Benefits.
*American Economic Review*, 84(3), 622-641.

## Examples

``` r
# DR-DDD with repeated cross-sections (one row per worker)
out <- ddd(yname = "log_wage", tname = "year", idname = "id",
           gname = "treatment_starts", pname = "treat_group",
           xformla = ~ nonwhite + union + education + experience,
           data = gruber1994, control_group = "nevertreated",
           base_period = "universal", est_method = "dr", panel = FALSE)
agg_ddd(out, type = "simple")
#>  Call:
#> agg_ddd(ddd_obj = out, type = "simple")
#> ========================= DDD Aggregation ============================
#>  Overall ATT:
#>    ATT    Std. Error     [ 95%  Conf. Int.] 
#>  0.004        0.0387    -0.0717      0.0798 
#> 
#>  Note: * indicates that the confidence interval does not contain zero.
#>  --------------------------- Data Info   -----------------------------
#>  Outcome variable: log_wage
#>  Qualification variable: treat_group
#>  Control group:  Never Treated
#>  --------------------------- Std. Errors  ----------------------------
#>  Level of significance:  0.05
#>  Analytical standard errors.
#>  =====================================================================
#>  See Ortiz-Villavicencio and Sant'Anna (2025) for details.
```
