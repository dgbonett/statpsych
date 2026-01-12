# Sample size for a 2-group proportion difference confidence interval

Computes the sample size in each group required to estimate a difference
of proportions with desired confidence interval precision in a 2-group
design. Set the proportion planning values to .5 for a conservatively
large sample size. Set R = 1 for equal sample sizes.

For more details, see Section 2.22 of Bonett (2021, Volume 3)

## Usage

``` r
size.ci.prop2(alpha, p1, p2, w, R)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- p1:

  planning value of proportion for group 1

- p2:

  planning value of proportion for group 2

- w:

  desired confidence interval width

- R:

  n2/n1 ratio

## Value

Returns the required sample size for each group

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.ci.prop2(.05, .05, .15, .1, 1)
#>   n1  n2
#>  269 269

# Should return:
#   n1  n2
#  269 269

size.ci.prop2(.05, .4, .2, .15, .5)
#>   n1  n2
#>  383 192

# Should return:
#   n1  n2
#  383 192

```
