# Confidence interval for an unknown population size

Computes a Wald confidence interval for an unknown population size using
mark-recapture sampling. This method assumes independence of the two
samples. This function requires the frequency counts from an incomplete
2 x 2 contingency table for the two samples (f11 is the unknown number
of people who were not observed in either sample). This method sets the
estimated odds ratio (with .5 added to each cell) to 1 and solves for
unobserved cell frequency. An approximate standard error is recovered
from the confidence interval.

For more details, see Section 3.7 of Bonett (2021, Volume 3)

## Usage

``` r
ci.popsize(alpha, f00, f01, f10)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f00:

  number of people observed in both samples

- f01:

  number of people observed in first sample but not second sample

- f10:

  number of people observed in second sample but not first sample

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimate of the unknown population size

- SE - recovered standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
ci.popsize(.05, 85, 196, 184)
#>  Estimate    SE  LL   UL
#>       889 67.35 777 1041

# Should return:
# Estimate    SE  LL   UL
#      889 67.35 777 1041

```
