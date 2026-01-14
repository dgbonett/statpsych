# Sample size for a coefficient of variation

Computes an approximate sample size required to estimate a population
coefficient of variation (CV) with desired confidence interval
precision. Set the CV planning value to the largest value within a
plausible range for a conservatively large sample size.

## Usage

``` r
size.ci.cv(alpha, CV, w)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- CV:

  planning value of coefficient of variation

- w:

  desired confidence interval width

## Value

Returns the required sample size

## Examples

``` r
size.ci.cv(.05, .25, .10)
#>  Sample size
#>           60

# Should return:
# Sample size
#          60
 
```
