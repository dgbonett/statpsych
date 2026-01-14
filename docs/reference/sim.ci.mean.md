# Simulates confidence interval coverage probability for a mean

Performs a computer simulation of the confidence interval performance
for a population mean. Sample data can be generated from five different
population distributions. All distributions are scaled to have a
standard deviation of 1.0.

## Usage

``` r
sim.ci.mean(alpha, n, dist, rep)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- n:

  sample size

- dist:

  type of distribution (1, 2, 3, 4,or 5)

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
  mean

- Lower Error - probability of lower limit greater than population mean

- Upper Error - probability of upper limit less than population mean

- Ave CI Width - average confidence interval width

## Examples

``` r
sim.ci.mean(.05, 10, 1, 5000)
#>  Coverage Lower Error Upper Error Ave CI Width
#>    0.9476      0.0262      0.0262     1.381611

# Should return (within sampling error):
# Coverage Lower Error Upper Error Ave CI Width
#   0.9484      0.0264      0.0252     1.392041

sim.ci.mean(.05, 40, 4, 1000)
#>  Coverage Lower Error Upper Error Ave CI Width
#>     0.929       0.017       0.054    0.6315832

# Should return (within sampling error):
# Coverage Lower Error Upper Error Ave CI Width
#  0.94722     0.01738      0.0354    0.6333067

```
