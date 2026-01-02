# Sample size for a 2-group proportion equivalence test

Computes the sample size in each group (assuming equal sample sizes)
required to perform an equivalence test for the difference in population
proportions with desired power in a 2-group design. The value of h
specifies a range of practical equivalence, -h to h, for the difference
in population proportions. The absolute difference in the proportion
planning values must be less than h. Equivalence tests often require a
very large sample size. Equivalence tests usually use 2 x alpha rather
than alpha (e.g., use alpha = .10 rather than alpha = .05). This
function sets the effect size equal to the difference in proportion
planning values.

For more details, see Section 2.23 of Bonett (2021, Volume 3)

## Usage

``` r
size.equiv.prop2(alpha, pow, p1, p2, h)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- p1:

  planning value of proportion for group 1

- p2:

  planning value of proportion for group 2

- h:

  upper limit for range of practical equivalence

## Value

Returns the required sample size for each group

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.equiv.prop2(.1, .9, .8, .8, .075)
#>  Sample size per group
#>                    488

# Should return:
# Sample size per group
#                   488

```
