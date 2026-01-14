# Sample size for a test of a 2-group proportion difference

Computes the sample size in each group required to test a difference in
population proportions with desired power (using a continuity
correction) in a 2-group design. This function requires planning values
for both proportions. Set each proportion planning value to .5 for a
conservatively large sample size requirement. This function does not
require the planning value for the proportion difference (effect size)
to equal the difference of the two proportion planning values; for
example, the planning value of the proportion difference could be set
equal to a minimally interesting effect size.

For more details, see Section 2.23 of Bonett (2021, Volume 3)

## Usage

``` r
size.test.prop2(alpha, pow, p1, p2, es)
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

- es:

  planning value of proportion difference (effect size)

## Value

Returns the required sample size for each group

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.test.prop2(.05, .8, .6, .4, .2)
#>  Sample size per group
#>                    106
# Should return:
# Sample size per group
#                   106

size.test.prop2(.05, .8, .5, .5, .2)
#>  Sample size per group
#>                    109
# Should return:
# Sample size per group
#                   109

```
