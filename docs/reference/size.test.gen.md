# Sample size for a test of any type of parameter

Computes the sample size required to test a single population parameter
with desired power using a standard error for the parameter estimate
from a prior or pilot study. This function can be used with any type of
parameter where the standard error of the parameter estimate is a
function of the square root of the sample size (most parameter estimates
have this property). This function also assumes that the sampling
distribution of the parameter estimate is approximately normal in large
samples.

For more details, see Section 2.29 of Bonett (2021, Volume 2)

## Usage

``` r
size.test.gen(alpha, pow, se, n0, es)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- se:

  standard error of parameter estimate from prior/pilot study

- n0:

  sample size of prior/pilot study

- es:

  planning value of parameter minus null hypothesis value

## Value

Returns the required sample size

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.test.gen(.05, .8, 2.89, 30, 5)
#>  Sample size
#>           79

# Should return:
# Sample size
#          79
 
```
