# Sample size for a mean absolute prediction error confidence interval

Computes the sample size required to estimate a population mean absolute
prediction error for a general linear model with desired confidence
interval precision. Setting s = 0 gives the sample size requirement for
a mean absolute deviation in a one-group design. This function assumes
that the prediction errors have an approximate normal distribution.

## Usage

``` r
size.ci.mape(alpha, mape, s, w)
```

## Arguments

- alpha:

  alpha value for 1-alpha confidence

- mape:

  mean absolute prediction error planning value

- s:

  number of predictor variables

- w:

  desired confidence interval width

## Value

Returns the required sample size

## Examples

``` r
size.ci.mape(.05, 4.5, 5, 2)
#>  Sample size
#>           57

# Should return:
# Sample size
#          57
 
```
