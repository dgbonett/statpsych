# Sample size for biserial-phi correlation confidence interval

Computes the sample size required to estimate a biserial-phi correlation
with desired confidence interval precision. Set the biserial-phi
planning value to the smallest absolute value within a plausible range
for a conservatively large sample size. The column variable is assumed
to be naturally dichotomous and the row variable is assumed to be
artificially dichotomous.

## Usage

``` r
size.ci.biphi(alpha, p1, p2, cor, w)
```

## Arguments

- alpha:

  alpha level for 1 - alpha confidence

- p1:

  planning value for row 1 marginal proportion

- p2:

  planning value for column 1 marginal proportion

- cor:

  planning value for biserial-phi correlation

- w:

  desired confidence interval width

## Value

Returns the required sample size

## Examples

``` r
size.ci.biphi(.05, .2, .5, .3, .4)
#>  Sample size
#>          195

# Should return:
#  Sample size
#          195

```
