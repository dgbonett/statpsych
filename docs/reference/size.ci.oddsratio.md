# Sample size for an odds ratio confidence interval

Computes the sample size required to estimate an odds ratio with desired
confidence interval precision.

## Usage

``` r
size.ci.oddsratio(alpha, p1, p2, or, r)
```

## Arguments

- alpha:

  alpha level for 1 - alpha confidence

- p1:

  planning value for row 1 marginal proportion

- p2:

  planning value for column 1 marginal proportion

- or:

  planning value of odds ratio

- r:

  desired upper to lower confidence interval endpoint ratio

## Value

Returns the required sample size

## Examples

``` r
size.ci.oddsratio(.05, .3, .2, 5.5, 3.0)
#>  Sample size
#>          356

# Should return:
#  Sample size
#          356

```
