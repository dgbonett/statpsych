# Confidence interval for a standardized mean

Computes a confidence interval for a population standardized mean
difference from a hypothesized value. If the hypothesized value is set
to 0, the reciprocals of the confidence interval endpoints gives a
confidence interval for the coefficient of variation (see
[ci.cv](https://dgbonett.github.io/statpsych/reference/ci.cv.md)).

## Usage

``` r
ci.stdmean(alpha, m, sd, n, h)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m:

  estimated mean

- sd:

  estimated standard deviation

- n:

  sample size

- h:

  hypothesized value of mean

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated standardized mean difference

- adj Estimate - bias adjusted standardized mean difference estimate

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2008). “Confidence intervals for standardized linear
contrasts of means.” *Psychological Methods*, **13**(2), 99–109. ISSN
1939-1463,
[doi:10.1037/1082-989X.13.2.99](https://doi.org/10.1037/1082-989X.13.2.99)
.

## Examples

``` r
ci.stdmean(.05, 24.5, 3.65, 40, 20)
#>  Estimate adj Estimate      SE     LL     UL
#>    1.2329        1.209 0.21244 0.8165 1.6493

# Should return:
# Estimate adj Estimate      SE     LL    UL
#   1.2329        1.209 0.21244 0.8165 1.6493
 
```
