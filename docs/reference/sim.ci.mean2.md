# Simulates confidence interval coverage probability for a 2-group mean difference

Performs a computer simulation of separate variance and pooled variance
confidence interval performance for a population mean difference in a
2-group design. Sample data within each group can be generated from five
different population distributions. All distributions are scaled to have
a standard deviation of 1.0 for group 1.

## Usage

``` r
sim.ci.mean2(alpha, n1, n2, sd2, dist1, dist2, rep)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n1:

  sample size in group 1

- n2:

  sample size in group 2

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
sim.ci.mean2(.05, 30, 25, 1.5, 1, 1, 1000)
#>                              Coverage Lower Error Upper Error Ave CI Width
#> Equal Variances Assumed:        0.934       0.026       0.040      1.34483
#> Equal Variances Not Assumed:    0.946       0.019       0.035      1.40174

# Should return (within sampling error):
#                              Coverage Lower Error Upper Error Ave CI Width
# Equal Variances Assumed:      0.93988      0.0322     0.02792     1.354437
# Equal Variances Not Assumed:  0.94904      0.0262     0.02476     1.411305

sim.ci.mean2(.05, 30, 25, 1.5, 4, 5, 1000)
#>                              Coverage Lower Error Upper Error Ave CI Width
#> Equal Variances Assumed:        0.937       0.049       0.014     1.351246
#> Equal Variances Not Assumed:    0.942       0.048       0.010     1.408630

# Should return (within sampling error):
#                              Coverage Lower Error Upper Error Ave CI Width
# Equal Variances Assumed:      0.93986     0.04022     0.01992     1.344437
# Equal Variances Not Assumed:  0.94762     0.03862     0.01376     1.401305

```
