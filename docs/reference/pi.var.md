# Prediction limits for a sample variance in a future study

Computes a two-sided or one-sided prediction limit for the estimated
variance in a future study for a planned sample size of n. The
prediction limit uses a variance estimate from a prior study of size n0.

Several confidence interval sample size functions in this package
require a planning value of the expected sample variance in the planned
study. A one-sided upper variance prediction limit is useful as a
variance planning value for a conservatively large sample size required
to obtain a confidence interval with desired width. This strategy for
specifying a variance planning value is useful in applications where the
population variance in the prior study is assumed to be very similar to
the population variance in the planned study.

For more details, see Section 1.31 of Bonett (2021, Volume 1)

## Usage

``` r
pi.var(alpha, var, n0, n, type)
```

## Arguments

- alpha:

  alpha value for upper 1-alpha confidence

- var:

  estimated variance from prior study

- n0:

  sample size used to estimate the variance

- n:

  planned sample size of future study

- type:

  - set to 1 for two-sided prediction interval

  - set to 2 for one-sided upper prediction limit

  - set to 3 for one-sided lower prediction limit

## Value

Returns two-sided or one-sided prediction limit(s) of an estimate
variance in a future study

## References

Hahn GJ (1972). “Simultaneous prediction intervals to contain the
standard deviations or ranges of future samples from a normal
population.” *Journal of the American Statistical Association*,
**67**(340), 938–942.
[doi:10.1080/01621459.1972.10481322](https://doi.org/10.1080/01621459.1972.10481322)
.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
pi.var(.05, 15, 40, 100, 2)
#>       UL
#>  23.9724

# Should return:
#      UL
# 23.9724
 
```
