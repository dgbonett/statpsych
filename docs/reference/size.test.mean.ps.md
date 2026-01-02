# Sample size for a test of a paired-samples mean difference

Computes the sample size required to test a difference in population
means with desired power in a paired-samples design. Set the Pearson
correlation planning value to the smallest value within a plausible
range, and set the variance planning value to the largest value within a
plausible range for a conservatively large sample size.

For more details, see Section 4.27 of Bonett (2021, Volume 1)

## Usage

``` r
size.test.mean.ps(alpha, pow, var, es, cor)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- var:

  planning value of average variance of the two measurements

- es:

  planning value of mean difference

- cor:

  planning value of correlation

## Value

Returns the required sample size

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.test.mean.ps(.05, .80, 1.25, .5, .75) 
#>  Sample size
#>           22

# Should return:
# Sample size
#          22
 
```
