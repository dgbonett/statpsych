# Sample size for a test of a within-subjects mean linear contrast

Computes the sample size required to test a linear contrast of
population means with desired power in a within-subjects design. Set the
average variance planning value to the largest value within a plausible
range for a conservatively large sample size. Set the average
correlation planning value to the smallest value within a plausible
range for a conservatively large sample size.

For more details, see Section 4.27 of Bonett (2021, Volume 1)

## Usage

``` r
size.test.lc.mean.ws(alpha, pow, var, es, cor, q)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- var:

  planning value of average variance of measurements

- es:

  planning value of linear contrast of means

- cor:

  planning value of average correlation between measurements

- q:

  vector of with-subjects contrast coefficients

## Value

Returns the required sample size

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
q <- c(1, -1, -1, 1)
size.test.lc.mean.ws(.05, .95, 15.0, 3, .8, q)
#>  Sample size
#>           20

# Should return:
# Sample size
#          20
 
```
