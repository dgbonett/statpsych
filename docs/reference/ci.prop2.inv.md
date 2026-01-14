# Confidence interval for a 2-group proportion difference using inverse sampling

Computes an approximate confidence interval for a population proportion
difference when inverse sampling has been used. An approximate standard
error is recovered from the confidence interval. With inverse sampling,
the number of participants who have the attribute within group 1 (f1)
and group 2 (f2) are predetermined, and sampling continues within each
group until f1 and f2 attain their prespecified values. With inverse
sampling, the sample sizes (n1 and n2) will not be known in advance.

## Usage

``` r
ci.prop2.inv(alpha, f1, f2, n1, n2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f1:

  number of participants in group 1 who have the attribute (fixed)

- f2:

  number of participants in group 2 who have the attribute (fixed)

- n1:

  sample size for group 1 (random)

- n2:

  sample size for group 2 (random)

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimate of proportion difference

- SE - recovered standard error

- LL - lower limit of confidence interval

- UL - upper limit of confidence interval

## References

Zou GY (2010). “Confidence interval estimation under inverse sampling.”
*Computational Statistics & Data Analysis*, **54**(1), 55–64. ISSN
0167-9473,
[doi:10.1016/j.csda.2005.05.007](https://doi.org/10.1016/j.csda.2005.05.007)
.

## Examples

``` r
ci.prop2.inv(.05, 10, 10, 48, 213)
#>  Estimate         SE         LL        UL
#>  0.161385 0.05997618 0.05288277 0.2879851

# Should return:
#  Estimate         SE         LL        UL
#  0.161385 0.05997618 0.05288277 0.2879851

```
