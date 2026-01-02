# Sample size for a squared multiple correlation confidence interval

Computes the sample size required to estimate a population squared
multiple correlation in a random-x regression model with desired
confidence interval precision. Set the planning value of the squared
multiple correlation to 1/3 for a conservatively large sample size.

For more details, see Section 2.28 of Bonett (2021, Volume 2)

## Usage

``` r
size.ci.rsqr(alpha, r2, s, w)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- r2:

  planning value of squared multiple correlation

- s:

  number of predictor variables in model

- w:

  desired confidence interval width

## Value

Returns the required sample size

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.ci.rsqr(.05, .3, 5, .2)
#>  Sample size
#>          227

# Should return:
# Sample size
#         227
 
```
