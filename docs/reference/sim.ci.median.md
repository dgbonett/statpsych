# Simulates confidence interval coverage probability for a median

Performs a computer simulation of the confidence interval performance
for a population median. Sample data can be generated from five
different population distributions. All distributions are scaled to have
a standard deviation of 1.0.

## Usage

``` r
sim.ci.median(alpha, n, dist, rep)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  sample size

- dist:

  type of distribution (1, 2, 3, 4, or 5)

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
  median

- Lower Error - probability of lower limit greater than population
  median

- Upper Error - probability of upper limit less than population median

- Ave CI Width - average confidence interval width

## Examples

``` r
sim.ci.median(.05, 20, 5, 1000)
#>  Coverage Lower Error Upper Error Ave CI Width
#>     0.947       0.032       0.021    0.9687742

# Should return (within sampling error):
# Coverage Lower Error Upper Error Ave CI Width
#   0.9589      0.0216      0.0195    0.9735528

```
