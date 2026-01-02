# Sample size for a 2-group mean ratio confidence interval

Computes the sample size in each group required to estimate a ratio of
population means with desired confidence interval precision in a 2-group
design. This function requires planning values for each mean and the
sample size requirement is very sensitive to these planning values. Set
the variance planning value to the largest value within a plausible
range for a conservatively large sample size. Set R = 1 for equal sample
sizes. For unequal sample sizes, this function assumes approximately
equal population variances.

For more details, see Section 2.13 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.ratio.mean2(alpha, var, m1, m2, r, R)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- var:

  planning value of average within-group variance

- m1:

  planning value of mean for group 1

- m2:

  planning value of mean for group 2

- r:

  desired upper to lower confidence interval endpoint ratio

- R:

  n2/n1 ratio

## Value

Returns the required sample size for each group

## References

Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.ci.ratio.mean2(.05, .4, 3.5, 3.1, 1.2, 1)
#>  n1 n2
#>  70 70

# Should return:
# n1   n2
# 70   70

size.ci.ratio.mean2(.05, .4, 3.5, 3.1, 1.2, 2)
#>  n1  n2
#>  53 106

# Should return:
# n1   n2
# 53  106
 
```
