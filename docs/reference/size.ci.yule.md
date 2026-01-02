# Sample size for a Yule's Q confidence interval

Computes the sample size required to estimate Yule's Q with desired
confidence interval precision. Set the Yule's Q planning value to the
smallest absolute value within a plausible range for a conservatively
large sample size.

## Usage

``` r
size.ci.yule(alpha, p1, p2, Q, w)
```

## Arguments

- alpha:

  alpha level for 1 - alpha confidence

- p1:

  planning value for row 1 marginal proportion

- p2:

  planning value for column 1 marginal proportion

- Q:

  planning value of Yule's Q

- w:

  desired confidence interval width

## Value

Returns the required sample size

## Examples

``` r
size.ci.yule(.05, .3, .2, .5, .4)
#>  Sample size
#>          354

# Should return:
#  Sample size
#          354

```
