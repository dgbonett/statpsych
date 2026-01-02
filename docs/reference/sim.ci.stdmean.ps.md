# Simulates confidence interval coverage probability for a standardized mean difference in a paired-samples design

Performs a computer simulation of confidence interval performance for
two types of standardized mean differences in a paired-samples design
(see ci.stdmean.ps). Sample data for the two levels of the
within-subjects factor can be generated from five different population
distributions. All distributions are scaled to have a standard deviation
of 1.0 at level 1. Bivariate random data with specified marginal
skewness and kurtosis are generated using the unonr function in the
mnonr package.

## Usage

``` r
sim.ci.stdmean.ps(alpha, n, sd2, cor, dist1, dist2, d, rep)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  sample size

- sd2:

  population standard deviation at level 2

- cor:

  correlation between paired measurements

- dist1:

  type of distribution at level 1 (1, 2, 3, 4, or 5)

- dist2:

  type of distribution at level 2 (1, 2, 3, 4, or 5)

  - 1 = Gaussian (skewness = 0 and excess kurtosis = 0)

  - 2 = platykurtic (skewness = 0 and excess kurtosis = -1.2)

  - 3 = leptokurtic (skewness = 0 and excess kurtosis = 6)

  - 4 = moderate skew (skewness = 1 and excess kurtosis = 1.5)

  - 5 = large skew (skewness = 2 and excess kurtosis = 6)

- d:

  population standardized mean difference

- rep:

  number of Monte Carlo samples

## Value

Returns a 1-row matrix. The columns are:

- Coverage - Probability of confidence interval including population std
  mean difference

- Lower Error - Probability of lower limit greater than population std
  mean difference

- Upper Error - Probability of upper limit less than population std mean
  difference

- Ave CI Width - Average confidence interval width

## Examples

``` r
sim.ci.stdmean.ps(.05, 20, 1.5, .8, 4, 4, .5, 2000)
#>                         Coverage Lower Error Upper Error Ave CI Width   Ave Est
#> Unweighted Standardizer   0.9120      0.0545      0.0335    0.7351741 0.5186650
#> Level 1 Standardizer      0.9405      0.0355      0.0240    0.9301617 0.5056967

# Should return (within sampling error):
#                         Coverage Lower Error Upper Error Ave CI Width   Ave Est
# Unweighted Standardizer   0.9095      0.0555       0.035    0.7354865 0.5186796
# Level 1 Standardizer      0.9525      0.0255       0.022    0.9330036 0.5058198

```
