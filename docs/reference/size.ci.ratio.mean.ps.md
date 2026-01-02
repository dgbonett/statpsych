# Sample size for a paired-samples mean ratio confidence interval

Computes the sample size required to estimate a ratio of population
means with desired confidence interval precision in a paired-samples
design. Set the correlation planning value to the smallest value within
a plausible range for a conservatively large sample size. This function
requires planning values for each mean and the sample size requirement
is very sensitive to these planning values. Set the variance planning
value to the largest value within a plausible range for a conservatively
large sample size.

For more details, see Section 4.26 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.ratio.mean.ps(alpha, var, m1, m2, cor, r)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- var:

  planning value of average variance of the two measurements

- m1:

  planning value of mean for measurement 1

- m2:

  planning value of mean for measurement 2

- cor:

  planning value for correlation between measurements

- r:

  desired upper to lower confidence interval endpoint ratio

## Value

Returns the required sample size

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.ci.ratio.mean.ps(.05, 400, 150, 100, .7, 1.2)
#>  Sample size
#>           21

# Should return:
# Sample size
#          21
 
```
