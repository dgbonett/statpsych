# Confidence interval for a median

Computes a distribution-free confidence interval for a population
median. Tied scores are assumed to be rare.

For more details, see Section 1.25 of Bonett (2021, Volume 1)

## Usage

``` r
ci.median(alpha, y)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- y:

  vector of scores

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated median

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Snedecor GW, Cochran WG (1989). *Statistical Methods*, 8th edition. ISU
University Pres, Ames, Iowa. Bonett DG (2021url). *Statistical Methods
for Psychologists*.

## Examples

``` r
y <- c(25, 29, 35, 36, 36, 40, 41, 43, 44, 54)
ci.median(.05, y)
#>  Estimate       SE LL UL
#>        38 3.261774 29 44

# Should return:
#  Estimate       SE  LL  UL
#        38 3.261774  29  44

```
