# Confidence interval for an odds ratio

Computes a confidence interval for an odds ratio with .5 added to each
cell frequency. This function requires the frequency counts from a 2 x 2
contingency table for two dichotomous variables.

For more details, see Section 3.4 of Bonett (2021, Volume 3)

## Usage

``` r
ci.oddsratio(alpha, f00, f01, f10, f11)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f00:

  number of participants with y = 0 and x = 0

- f01:

  number of participants with y = 0 and x = 1

- f10:

  number of participants with y = 1 and x = 0

- f11:

  number of participants with y = 1 and x = 1

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimate of odds ratio

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Fleiss JL, Paik MC (2003). *Statistical Methods for Rates and
Proportions*, 3rd edition. Wiley. Bonett DG (2021url). *Statistical
Methods for Psychologists*.

## Examples

``` r
ci.oddsratio(.05, 229, 28, 96, 24)
#>  Estimate        SE       LL       UL
#>  2.044451 0.6154578 1.133267 3.688254

# Should return:
#  Estimate        SE       LL       UL
#  2.044451 0.6154578 1.133267 3.688254

```
