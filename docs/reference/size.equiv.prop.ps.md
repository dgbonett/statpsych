# Sample size for a paired-samples proportion equivalence test

Computes the sample size required to perform an equivalence test for the
difference in population proportions with desired power in a
paired-samples design. The value of h specifies a range of practical
equivalence, -h to h, for the difference in population proportions. The
absolute difference in the proportion planning values must be less than
h. Equivalence tests often require a very large sample size. Equivalence
tests usually use 2 x alpha rather than alpha (e.g., use alpha = .10
rather alpha = .05). This function sets the effect size equal to the
difference in proportion planning values. Set the phi correlation
planning value to the smallest absolute value within a plausible range
for a conservatively large sample size.

For more details, see Section 3.13 of Bonett (2021, Volume 3)

## Usage

``` r
size.equiv.prop.ps(alpha, pow, p1, p2, phi, h)
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

- h:

  upper limit for range of practical equivalence

## Value

Returns the required sample size

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.equiv.prop.ps(.1, .8, .30, .35, .40, .15)
#>  Sample size
#>          173

# Should return:
# Sample size
#         173

```
