# Sample size for a paired-samples mean equivalence test

Computes the sample size required to perform an equivalence test for the
difference in population means with desired power in a paired-samples
design. The value of h specifies a range of practical equivalence, -h to
h, for the difference in population means. The planning value for the
absolute mean difference must be less than h. Equivalence tests often
require a very large sample size. Equivalence tests usually use 2 x
alpha rather than alpha (e.g., use alpha = .10 rather alpha = .05). Set
the Pearson correlation value to the smallest value within a plausible
range, and set the variance planning value to the largest value within a
plausible range for a conservatively large sample size.

For more details, see Section 4.27 of Bonett (2021, Volume 1)

## Usage

``` r
size.equiv.mean.ps(alpha, pow, var, es, cor, h)
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

  planning value of the correlation between measurements

- h:

  upper limit for range of practical equivalence

## Value

Returns the required sample size

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.equiv.mean.ps(.10, .90, 25, .5, .75, 2)
#>  Sample size
#>           49

# Should return:
# Sample size
#          49
 
```
