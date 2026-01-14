# Confidence interval for an indirect effect

Computes a Monte Carlo confidence interval (500,000 trials) for a
population unstandardized or standardized indirect effect in a path
model and a Sobel standard error. This function is not recommended for a
standardized indirect if the standardized slopes are greater than .4 The
Monte Carlo method is general in that the slope estimates and standard
errors do not need to be OLS estimates with homoscedastic standard
errors. For example, LAD slope estimates and their standard errors, OLS
slope estimates and heteroscedastic-consistent standard errors also
could be used. In models with no direct effects, distribution-free
Theil-Sen slope estimates with recovered standard errors (see
[ci.theil](https://dgbonett.github.io/statpsych/reference/ci.theil.md))
also could be used.

For more details, see Section 1.18 of Bonett (2021, Volume 4)

## Usage

``` r
ci.indirect(alpha, b1, b2, se1, se2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- b1:

  slope estimate for first path

- b2:

  slope estimate for second path

- se1:

  standard error for b1

- se2:

  standard error for b2

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated indirect effect

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.indirect (.05, 2.48, 1.92, .586, .379)
#>  Estimate       SE       LL       UL
#>    4.7616 1.466064 2.165953 7.955538

# Should return (within sampling error):
# Estimate       SE       LL       UL
#   4.7616 1.625282 2.178812 7.972262
 
```
