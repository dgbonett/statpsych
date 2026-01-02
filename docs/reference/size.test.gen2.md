# Sample size for a test of 2-group difference for any type of parameter

Computes the sample size per group required to test a difference in two
population parameters with desired power using a standard error for a
single parameter estimate from a prior or pilot study. This function can
be used with any type of parameter where the standard error of the
parameter estimate is a function of the square root of the sample size
(most parameter estimates have this property). This function also
assumes that the sampling distribution of the parameter estimate is
approximately normal in large samples. Set R = 1 for equal sample sizes.

## Usage

``` r
size.test.gen2(alpha, pow, se, n0, es, R)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- se:

  standard error of parameter estimate from prior/pilot study

- n0:

  sample size of prior/pilot study

- es:

  planning value of parameter difference

- R:

  n2/n1 ratio

## Value

Returns the required sample size for each group

## Examples

``` r
size.test.gen2(.05, .85, .175, 30, .5, 1)
#>  n1 n2
#>  66 66

# Should return:
# n1  n2
# 66  66
 
```
