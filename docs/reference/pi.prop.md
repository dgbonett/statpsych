# Prediction interval for a sample proportion in a future study

Computes an approximate one-sided or two-sided prediction interval for
the estimated proportion in a future study with a planned sample size of
n. The prediction interval uses a proportion estimate from a prior study
that had a sample size of n0.

Several confidence interval sample size functions in this package
require a planning value of the expected sample value of a proportion in
the planned study. A one-sided proportion prediction limit is useful as
a proportion planning value for a conservatively large sample size
required to obtain a confidence interval with desired width. This
strategy for specifying a proportion planning value is useful in
applications where the population proportion in the prior study is
assumed to be very similar to the population proportion in the planned
study.

For sample size planning, use an upper prediction limit if the
population proportion is assumed to be less than .5; and if the upper
prediction limit is greater than .5, set the proportion planning value
to .5. Use a lower prediction limit if the population proportion is
asumed to be greater than .5; and if the lower prediction limit is less
than .5, set the proportion planning value to .5.

For more details, see Section 1.16 of Bonett (2021, Volume 3)

## Usage

``` r
pi.prop(alpha, prop, n0, n, type)
```

## Arguments

- alpha:

  alpha value for 1-alpha confidence

- prop:

  estimated proportion from prior study

- n0:

  sample size used to estimate proportion in prior study

- n:

  planned sample size of future study

- type:

  - set to 1 for two-sided prediction interval

  - set to 2 for one-sided upper prediction limit

  - set to 3 for one-sided lower prediction limit

## Value

Returns one-sided or two-sided prediction limit(s) for an estimated
proportion in a future study

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
pi.prop(.1, .225, 80, 120, 1)
#>         LL       UL
#>  0.1390955 0.337095

# Should return:
#         LL       UL
#  0.1390955 0.337095
 
```
