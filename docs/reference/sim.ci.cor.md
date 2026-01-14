# Simulates confidence interval coverage probability for a Pearson correlation

Performs a computer simulation of confidence interval performance for a
Pearson correlation. A bias adjustment is used to reduce the bias of the
Fisher transformed Pearson correlation. Sample data can be generated
from bivariate population distributions with five different marginal
distributions. All distributions are scaled to have standard deviations
of 1.0. Bivariate random data with specified marginal skewness and
kurtosis are generated using the unonr function in the mnonr package.

## Usage

``` r
sim.ci.cor(alpha, n, cor, dist1, dist2, rep)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  sample size

- cor:

  population Pearson correlation

- dist1:

  type of distribution for variable 1 (1, 2, 3, 4, or 5)

- dist2:

  type of distribution for variable 2 (1, 2, 3, 4, or 5)

  - 1 = Gaussian (skewness = 0 and excess kurtosis = 0)

  - 2 = platykurtic (skewness = 0 and excess kurtosis = -1.2)

  - 3 = leptokurtic (skewness = 0 and excess kurtosis = 6)

  - 4 = moderate skew (skewness = 1 and excess kurtosis = 1.5)

  - 5 = large skew (skewness = 2 and excess kurtosis = 6)

- rep:

  number of Monte Carlo samples

## Value

Returns a 1-row matrix. The columns are:

- Coverage - probability of confidence interval including population
  correlation

- Lower Error - probability of lower limit greater than population
  correlation

- Upper Error - probability of upper limit less than population
  correlation

- Ave CI Width - average confidence interval width

## Examples

``` r
sim.ci.cor(.05, 30, .7, 4, 5, 1000)
#>      Coverage Lower Error Upper Error Ave CI Width
#> [1,]    0.934       0.046        0.02    0.3881714

# Should return (within sampling error):
#      Coverage Lower Error Upper Error Ave CI Width
# [1,]  0.93815     0.05125      0.0106    0.7778518

```
