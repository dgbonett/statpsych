# Sample size for a test of a paired-samples proportion difference

Computes the sample size required to test a difference in population
proportions with desired power in a paired-samples design. This function
requires planning values for both proportions and a phi coefficient that
describes the correlation between the two dichotomous measurements. The
proportion planning values can be set to .5 for a conservatively large
sample size. The planning value for the effect size (proportion
difference) could be set equal to the difference of the two proportion
planning values or it could be set equal to a minimally interesting
effect size. Set the phi correlation planning value to the smallest
absolute value within a plausible range for a conservatively large
sample size.

For more details, see Section 3.13 of Bonett (2021, Volume 3)

## Usage

``` r
size.test.prop.ps(alpha, pow, p1, p2, phi, es)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- p1:

  planning value of proportion for measurement 1

- p2:

  planning value of proportion for measurement 2

- phi:

  planning value of phi correlation

- es:

  planning value of proportion difference

## Value

Returns the required sample size

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.test.prop.ps(.05, .90, .6, .7, .5, .1)
#>  Sample size
#>          237

# Should return:
# Sample size
#         237

```
