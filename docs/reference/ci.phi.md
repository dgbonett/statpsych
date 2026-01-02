# Confidence interval for a phi correlation

Computes a Fisher confidence interval for a population phi correlation.
This function requires the frequency counts from a 2 x 2 contingency
table for two dichotomous variables. This measure of association is
usually most appropriate when both dichotomous variables are naturally
dichotomous.

For more details, see Section 3.4 of Bonett (2021, Volume 3)

## Usage

``` r
ci.phi(alpha, f00, f01, f10, f11)
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

- Estimate - estimate of phi correlation

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bishop YMM, Fienberg SE, Holland PW (1975). *Discrete Multivariate
Analysis*. MIT Press. Bonett DG (2021url). *Statistical Methods for
Psychologists*.

## Examples

``` r
ci.phi(.05, 229, 28, 96, 24)
#>  Estimate     SE    LL    UL
#>     0.123 0.0548 0.015 0.229

# Should return:
#  Estimate     SE    LL    UL
#     0.123 0.0548 0.015 0.229

```
