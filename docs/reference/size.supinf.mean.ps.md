# Sample size for a paired-samples mean superiority or noninferiority test

Computes the sample size required to perform a superiority or
noninferiority test for the difference in population means with desired
power in a paired-samples design. For a superiority test, specify the
upper limit (h) for the range of practical equivalence and specify an
effect size (es) such that es \> h. For a noninferiority test, specify
the lower limit (-h) for the range of practical equivalence and specify
an effect size such that es \> -h. Set the Pearson correlation planning
value to the smallest value within a plausible range, and set the
variance planning value to the largest value within a plausible range
for a conservatively large sample size.

For more details, see Section 4.27 of Bonett (2021, Volume 1)

## Usage

``` r
size.supinf.mean.ps(alpha, pow, var, es, cor, h)
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

  upper or lower limit for range of practical equivalence

## Value

Returns the required sample size

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.supinf.mean.ps(.05, .80, 225, 9, .75, 4)
#>  Sample size
#>           38

# Should return:
# Sample size
#          38
 
```
