# Confidence interval for squared multiple correlation

Computes an approximate confidence interval for a population squared
multiple correlation in a linear model with random predictor variables.
This function uses the scaled central F approximation method. An
approximate standard error is recovered from the confidence interval.

For more details, see Section 2.4 of Bonett (2021, Volume 2)

## Usage

``` r
ci.rsqr(alpha, r2, s, n)
```

## Arguments

- alpha:

  alpha value for 1-alpha confidence

- r2:

  estimated unadjusted squared multiple correlation

- s:

  number of predictor variables

- n:

  sample size

## Value

Returns a 1-row matrix. The columns are:

- R-squared - estimate of unadjusted R-squared (from input)

- adj R-squared - bias adjusted R-squared estimate

- SE - recovered standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Helland IS (1987). “On the interpretation and use of R2 in regression
analysis.” *Biometrics*, **43**(1), 61–69.
[doi:10.2307/2531949](https://doi.org/10.2307/2531949) .

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.rsqr(.05, .247, 4, 150)
#>  R-squared adj R-squared      SE     LL     UL
#>      0.247        0.2262 0.06024 0.1152 0.3514

# Should return:
# R-squared adj R-squared      SE      LL     UL  
#     0.247        0.2262 0.06024  0.1152 0.3514
 
```
