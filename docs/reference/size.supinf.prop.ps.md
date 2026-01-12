# Sample size for a paired-samples superiority or inferiority test of proportions

Computes the sample size required to perform a superiority or
inferiority test for the difference in population proportions with
desired power in a paired-samples design. For a superiority test,
specify the upper limit (h) for the range of practical equivalence and
specify values of p1 and p2 such that p1 - p2 \> h. For an inferiority
test, specify the lower limit (-h) for the range of practical
equivalence and specify values of p1 and p2 such that p1 - p2 \> -h.
This function sets the effect size equal to p1 - p2. Set the phi
correlation planning value to the smallest absolute value within a
plausible range for a conservatively large sample size.

For more details, see Section 3.13 of Bonett (2021, Volume 3)

## Usage

``` r
size.supinf.prop.ps(alpha, pow, p1, p2, phi, h)
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

  lower or upper limit for range of practical equivalence

## Value

Returns the required sample size

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.supinf.prop.ps(.05, .9, .35, .20, .45, .05)
#>  Sample size
#>          227

# Should return:
# Sample size
#         227

```
