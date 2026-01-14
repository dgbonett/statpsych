# Confidence interval for a standard deviation

Computes the traditional confidence interval for a population standard
deviation using a sample estimate of the standard deviation. The
traditional confidence interval assumes normality and is hypersensitive
to minor violations of this assumption. This function should be used
only if the data appear to come from an approximate normal or mildly
platykurtic distribution.

## Usage

``` r
ci.sd(alpha, sd, n)
```

## Arguments

- alpha:

  alpha value for 1-alpha confidence

- sd:

  estimated standard deviation

- n:

  sample size

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated standard deviation (from input)

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## Examples

``` r
ci.sd(.05, 4.65, 50)
#>  Estimate       LL      UL
#>      4.65 3.884303 5.79452

# Should return:
#  Estimate       LL       UL
#      4.65 3.884303  5.79452
 
```
