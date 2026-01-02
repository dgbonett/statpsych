# Confidence interval for a ratio of coefficients of variation in a 2-group design

Computes a confidence interval for a ratio of population coefficients of
variation (CV) in a 2-group design. This confidence interval uses the
confidence interval for each CV and then uses the MOVER-DL method (see
Newcombe, page 138) to obtain a confidence interval for CV1/CV2. The CV
assumes ratio-scale scores.

## Usage

``` r
ci.ratio.cv2(alpha, m1, m2, sd1, sd2, n1, n2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m1:

  estimated mean for group 1

- m2:

  estimated mean for group 2

- sd1:

  estimated standard deviation for group 1

- sd2:

  estimated standard deviation for group 2

- n1:

  sample size for group 1

- n2:

  sample size for group 2

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated ratio of coefficients of variation

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Newcombe RG (2013). *Confidence Interval for Proportions and Related
Measures of Effect Size*. CRC Press.

## Examples

``` r
ci.ratio.cv2(.05, 34.5, 26.1, 4.15, 2.26, 50, 50)
#>  Estimate       LL       UL
#>  1.389188 1.041478 1.854101

# Should return:
# Estimate       LL       UL
# 1.389188 1.041478 1.854101
 
```
