# Sample size for a confidence interval for any type of parameter

Computes the sample size required to estimate a single population
parameter with desired precision using a standard error for the
parameter estimate from a prior or pilot study. This function can be
used with any type of parameter where the standard error of the
parameter estimate is a function of the square root of the sample size
(most parameter estimates have this property). This function also
assumes that the sampling distribution of the parameter estimate is
approximately normal in large samples.

## Usage

``` r
size.ci.gen(alpha, se, n0, w)
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

## Value

Returns the required sample size

## Examples

``` r
size.ci.gen(.05, 2.89, 30, 8)
#>  Sample size
#>           61

# Should return:
# Sample size
#          61
 
```
