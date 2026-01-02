# Sample size for phi correlation confidence interval

Computes the sample size required to estimate a phi correlation with
desired confidence interval precision. Set the phi correlation planning
value to the smallest absolute value within a plausible range for a
conservatively large sample size.

## Usage

``` r
size.ci.phi(alpha, p1, p2, phi, w)
```

## Arguments

- alpha:

  alpha level for 1 - alpha confidence

- p1:

  planning value for row 1 marginal proportion

- p2:

  planning value for column 1 marginal proportion

- phi:

  planning value for phi correlation

- w:

  desired confidence interval width

## Value

Returns the required sample size

## Examples

``` r
size.ci.phi(.05, .7, .8, .35, .2)
#>  Sample size
#>          418

# Should return:
#  Sample size
#          416

```
