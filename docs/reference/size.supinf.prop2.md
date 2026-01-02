# Sample size for a 2-group superiority or inferiority test of proportions

Computes the sample size in each group (assuming equal sample sizes)
required to perform a superiority or inferiority test for the difference
in population proportions with desired power in a 2-group design. For a
superiority test, specify the upper limit (h) for the range of practical
equivalence and specify values of p1 and p2 such that p1 - p2 \> h. For
an inferiority test, specify the lower limit (-h) for the range of
practical equivalence and specify values of p1 and p2 such that p1 - p2
\> -h. This function sets the effect size equal to p1 - p2.

For more details, see Section 2.23 of Bonett (2021, Volume 3)

## Usage

``` r
size.supinf.prop2(alpha, pow, p1, p2, h)
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

  lower or upper limit for range of practical equivalence

## Value

Returns the required sample size for each group

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.supinf.prop2(.05, .9, .35, .20, .05)
#>  Sample size per group
#>                    408

# Should return:
# Sample size per group
#                   408

```
