# Prediction limits for a sample correlation in a future study

Computes approximate one-sided or two-sided prediction limits for the
estimated Pearson correlation in a future study with a planned sample
size of n. The prediction interval uses a correlation estimate from a
prior study that had a sample size of n0.

Several confidence interval sample size functions in this package
require a planning value of the expected sample value of a Pearson
correlation in the planned study. A one-sided lower correlation
prediction limit is useful as a correlation planning value for a
conservatively large sample size required to obtain a confidence
interval with desired width. This strategy for specifying a correlation
planning value is useful in applications where the population
correlation in the prior study is assumed to be very similar to the
population correlation in the planned study.

For more details, see Section 1.26 of Bonett (2021, Volume 2)

## Usage

``` r
pi.cor(alpha, cor, n0, n, type)
```

## Arguments

- alpha:

  alpha value for 1-alpha confidence

- cor:

  estimated Pearson correlation from prior study

- n0:

  sample size used to estimate the correlation in prior study

- n:

  planned sample size of future study

- type:

  - set to 1 for two-sided prediction interval

  - set to 2 for one-sided upper prediction limit

  - set to 3 for one-sided lower prediction limit

## Value

Returns one-sided or two-sided prediction limit(s) of an estimated
Pearson correlation in a future study

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
pi.cor(.1, .761, 50, 100, 1)
#>      LL     UL
#>  0.6034 0.8573

# Should return:
#      LL     UL
#  0.6034 0.8573
 
pi.cor(.1, .761, 50, 100, 3)
#>      LL
#>  0.6429

# Should return:
#      LL
#  0.6429
 
```
