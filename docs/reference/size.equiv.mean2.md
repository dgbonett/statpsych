# Sample size for a 2-group mean equivalence test

Computes the sample size in each group (assuming equal sample sizes)
required to perform an equivalence test for the difference in population
means with desired power in a 2-group design. The value of h specifies a
range of practical equivalence, -h to h, for the difference in
population means. The planning value for the absolute mean difference
must be less than h. Equivalence tests often require a very large sample
size. Equivalence tests usually use 2 x alpha rather than alpha (e.g.,
use alpha = .10 rather alpha = .05). Set the variance planning value to
the largest value within a plausible range for a conservatively large
sample size.

For more details, see Section 2.14 of Bonett (2021, Volume 1)

## Usage

``` r
size.equiv.mean2(alpha, pow, var, es, h)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- var:

  planning value of average within-group variance

- es:

  planning value of mean difference

- h:

  upper limit for range of practical equivalence

## Value

Returns the required sample size for each group

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.equiv.mean2(.10, .80, 15, 2, 4)
#>  Sample size per group
#>                     50

# Should return:
# Sample size per group 
#                    50
 
```
