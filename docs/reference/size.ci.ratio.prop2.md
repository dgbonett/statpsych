# Sample size for a 2-group proportion ratio confidence interval

Computes the sample size in each group required to estimate a ratio of
proportions with desired confidence interval precision in a 2-group
design. Set R = 1 for equal sample sizes.

For more details, see Section 2.22 of Bonett (2021, Volume 3)

## Usage

``` r
size.ci.ratio.prop2(alpha, p1, p2, r, R)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- p1:

  planning value of proportion for group 1

- p2:

  planning value of proportion for group 2

- r:

  desired upper to lower confidence interval endpoint ratio

- R:

  n2/n1 ratio

## Value

Returns the required sample size for each group

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.ci.ratio.prop2(.05, .9, .7, 1.5, 1)
#>  n1 n2
#>  51 51

# Should return:
#   n1  n2
#   51  51

size.ci.ratio.prop2(.05, .9, .7, 1.5, .5)
#>  n1 n2
#>  91 46

# Should return:
#   n1  n2
#   91  46

```
