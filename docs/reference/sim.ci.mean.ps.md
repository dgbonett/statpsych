# Simulates confidence interval coverage probability for a paired-samples mean difference

Performs a computer simulation of confidence interval performance for a
population mean difference in a paired-samples design. Sample data for
the two levels of the within-subjects factor can be generated from
bivariate population distributions with five different marginal
distributions. All distributions are scaled to have a standard deviation
of 1.0 at level 1. Bivariate random data with specified marginal
skewness and kurtosis are generated using the unonr function in the
mnonr package.

## Usage

``` r
sim.ci.mean.ps(alpha, n, sd2, cor, dist1, dist2, rep)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  sample size

- sd2:

  population standard deviation at level 2

- cor:

  population correlation of paired observations

- dist1:

  type of distribution at level 1 (1, 2, 3, 4, or 5)

- dist2:

  type of distribution at level 2 (1, 2, 3, 4, or 5)

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
  mean difference

- Lower Error - probability of lower limit greater than population mean
  difference

- Upper Error - probability of upper limit less than population mean
  difference

- Ave CI Width - average confidence interval width

## Examples

``` r
sim.ci.mean.ps(.05, 30, 1.5, .7, 4, 5, 1000)
#>  Coverage Lower Error Upper Error Ave CI Width
#>      0.95       0.041       0.009    0.7842072

# Should return (within sampling error):
# Coverage Lower Error Upper Error Ave CI Width
#  0.94415     0.04525      0.0106    0.7818518

```
