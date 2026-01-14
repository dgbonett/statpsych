# Confidence interval for a linear contrast of parameters in a between-subjects design

Computes the estimate, standard error, and approximate confidence
interval for a linear contrast of any type of parameter where each
parameter value has been estimated from a different sample. The
parameter values are assumed to be of the same type and their sampling
distributions are assumed to be approximately normal.

## Usage

``` r
ci.lc.gen.bs(alpha, est, se, v)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- est:

  vector of parameter estimates

- se:

  vector of standard errors

- v:

  vector of contrast coefficients

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimate of linear contrast

- SE - standard error of linear contrast

- LL - lower limit of confidence interval

- UL - upper limit of confidence interval

## Examples

``` r
est <- c(3.86, 4.57, 2.29, 2.88)
se <- c(0.185, 0.365, 0.275, 0.148)
v <- c(.5, .5, -.5, -.5)
ci.lc.gen.bs(.05, est, se, v)
#>  Estimate        SE       LL       UL
#>      1.63 0.2573806 1.125543 2.134457

# Should return:
# Estimate        SE       LL       UL
#     1.63 0.2573806 1.125543 2.134457

```
