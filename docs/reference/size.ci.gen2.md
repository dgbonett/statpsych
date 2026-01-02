# Sample size for a confidence interval for the difference of any type of parameter

Computes the sample size required to estimate a difference in population
parameters with desired precision in a 2-group design using a standard
error for a parameter estimate from a prior or pilot study. This
function can be used with any type of parameter where the standard error
of the parameter estimate is a function of the square root of the sample
size (most parameter estimates have this property). This function also
assumes that the sampling distribution of the parameter estimate is
approximately normal in large samples. Set R = 1 for equal sample sizes.

## Usage

``` r
size.ci.gen2(alpha, se, n0, w, R)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- se:

  standard error of parameter estimate from prior/pilot study

- n0:

  sample size of prior/pilot study

- w:

  desired confidence interval width

- R:

  n2/n1 ratio

## Value

Returns the required sample size for each group

## Examples

``` r
size.ci.gen2(.05, .175, 30, .8, 1)
#>  n1 n2
#>  45 45

# Should return:
# n1  n2
# 45  45
 
```
