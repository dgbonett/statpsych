# Simulates confidence interval coverage probability for a median difference in a 2-group design

Performs a computer simulation of the confidence interval performance
for a difference of population medians in a 2-group design. Sample data
for each group can be generated from five different population
distributions. All distributions are scaled to have a standard deviation
of 1.0 for group 1.

## Usage

``` r
sim.ci.median2(alpha, n1, n2, sd2, dist1, dist2, rep)
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

- rep:

  number of Monte Carlo samples

## Value

Returns a 1-row matrix. The columns are:

- Coverage - Probability of confidence interval including population
  median difference

- Lower Error - Probability of lower limit greater than population
  median difference

- Upper Error - Probability of upper limit less than population median
  difference

- Ave CI Width - Average confidence interval width

## Examples

``` r
sim.ci.median2(.05, 20, 20, 2, 5, 4, 5000)
#>  Coverage Lower Error Upper Error Ave CI Width
#>    0.9478      0.0306      0.0216     2.386471

# Should return (within sampling error):
# Coverage Lower Error Upper Error Ave CI Width
#    0.952       0.027       0.021     2.368914

```
