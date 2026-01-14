# Sample size for a standard deviation confidence interval

Computes the sample size required to estimate a population standard
deviation with desired confidence interval precision. This function
assumes that the traditional confidence interval for a population
standard deviation will be used. The traditional confidence interval
assumes that the response variable has an approximate normal
distribution and can be highly innacurate when this assumption is not
satisfied.

## Usage

``` r
size.ci.sd(alpha, r)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- r:

  desired upper to lower confidence interval endpoint ratio

## Value

Returns the required sample size

## Examples

``` r
size.ci.sd(.05, 1.5)
#>  Sample size
#>           49

# Should return:
# Sample size
#          49
 
```
