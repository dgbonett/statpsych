# Sample size for a test of a between-subjects mean linear contrast

Computes the sample size in each group (assuming equal sample sizes)
required to test a linear contrast of population means with desired
power in a between-subjects design. Set the variance planning value to
the largest value within a plausible range for a conservatively large
sample size.

For more details, see Section 3.25 of Bonett (2021, Volume 1)

## Usage

``` r
size.test.lc.mean.bs(alpha, pow, var, es, v)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- var:

  planning value of average within-group variance

- es:

  planning value of linear contrast of means

- v:

  vector of between-subjects contrast coefficients

## Value

Returns the required sample size for each group

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
v <- c(1/4, 1/4, 1/4, 1/4, -1)
size.test.lc.mean.bs(.05, .90, 1, .5, v)
#>  Sample size per group
#>                     53

# Should return:
# Sample size per group
#                    53

v <- c(1, -1, -1, 1)
size.test.lc.mean.bs(.05, .90, 27.5, 5, v)
#>  Sample size per group
#>                     47

# Should return:
# Sample size per group
#                    47
 
```
