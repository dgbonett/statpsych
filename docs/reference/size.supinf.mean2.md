# Sample size for a 2-group mean superiority or noninferiority test

Computes the sample size in each group (assuming equal sample sizes)
required to perform a superiority or noninferiority test for the
difference in population means with desired power in a 2-group design.
For a superiority test, specify the upper limit (h) for the range of
practical equivalence and specify an effect size (es) such that es \> h.
For a noninferiority test, specify the lower limit (-h) for the range of
practical equivalence and specify an effect size such that es \> -h. Set
the variance planning value to the largest value within a plausible
range for a conservatively large sample size.

For more details, see Section 2.14 of Bonett (2021, Volume 1)

## Usage

``` r
size.supinf.mean2(alpha, pow, var, es, h)
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

  upper or lower limit for range of practical equivalence

## Value

Returns the required sample size for each group

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.supinf.mean2(.05, .80, 225, 9, 4)
#>  Sample size per group
#>                    143

# Should return:
# Sample size per group 
#                   143
 
```
