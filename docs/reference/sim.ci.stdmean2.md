# Simulates confidence interval coverage probability for a standardized mean difference in a 2-group design

Performs a computer simulation of confidence interval performance for
two types of standardized mean differences in a 2-group design (see
ci.stdmean2). Sample data for each group can be generated from five
different population distributions. All distributions are scaled to have
a standard deviation of 1.0 for group 1.

## Usage

``` r
sim.ci.stdmean2(alpha, n1, n2, sd2, dist1, dist2, d, rep)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n1:

  sample size for group 1

- n2:

  sample size for group 2

- sd2:

  population standard deviation for group 2

- dist1:

  type of distribution for group 1 (1, 2, 3, 4, or 5)

- dist2:

  type of distribution for group 2 (1, 2, 3, 4, or 5)

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
sim.ci.stdmean2(.05, 20, 20, 1.5, 3, 4, .75, 5000)
#>                         Coverage Lower Error Upper Error Ave CI Width   Ave Est
#> Unweighted Standardizer   0.9178       0.052      0.0302     1.341802 0.7837455
#> Group 1 Standardizer      0.9476       0.028      0.0244     1.814108 0.7831089

# Should return (within sampling error):
#                         Coverage Lower Error Upper Error Ave CI Width   Ave Est
# Unweighted Standardizer   0.9058      0.0610      0.0332     1.342560 0.7838679
# Group 1 Standardizer      0.9450      0.0322      0.0228     1.827583 0.7862640

```
