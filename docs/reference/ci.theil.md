# Theil-Sen estimate and confidence interval for slope

Computes a Theil-Sen estimate and distribution-free confidence interval
for the population slope in a simple linear regression model. An
approximate standard error is recovered from the confidence interval.

For more details, see Section 1.32 of Bonett (2021, Volume 2)

## Usage

``` r
ci.theil(alpha, y, x)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- y:

  vector of response variable scores

- x:

  vector of predictor variable scores (paired with y)

## Value

Returns a 1-row matrix. The columns are:

- Estimate - Theil-Sen estimate of population slope

- SE - recovered standard error

- LL - lower limit of confidence interval

- UL - upper limit of confidence interval

## References

Hollander M, Wolf DA (1999). *Nonparametric Statistical Methods*. Wiley.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
y <- c(21, 4, 9, 12, 35, 18, 10, 22, 24, 1, 6, 8, 13, 16, 19)
x <- c(67, 28, 30, 28, 52, 40, 25, 37, 44, 10, 14, 20, 28, 40, 51)
ci.theil(.05, y, x)
#>  Estimate        SE        LL   UL
#>       0.5 0.1085927 0.3243243 0.75

# Should return:
# Estimate        SE        LL   UL
#      0.5 0.1085927 0.3243243 0.75

```
