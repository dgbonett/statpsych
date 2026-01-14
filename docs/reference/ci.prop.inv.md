# Confidence interval for a proportion using inverse sampling

Computes an exact confidence interval for a population proportion when
inverse sampling has been used. An approximate standard error is
recovered from the confidence interval. With inverse sampling, the
number of participants who have the attribute (f) is predetermined and
sampling continues until f attains its prespecified value. With inverse
sampling, the sample size (n) will not be known in advance.

For more details, see Section 1.19 of Bonett (2021, Volume 3)

## Usage

``` r
ci.prop.inv(alpha, f, n)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- f:

  number of participants who have the attribute (fixed)

- n:

  sample size (random)

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimate of proportion

- SE - recovered standard error

- LL - lower limit of confidence interval

- UL - upper limit of confidence interval

## References

Zou GY (2010). “Confidence interval estimation under inverse sampling.”
*Computational Statistics & Data Analysis*, **54**(1), 55–64. ISSN
0167-9473,
[doi:10.1016/j.csda.2005.05.007](https://doi.org/10.1016/j.csda.2005.05.007)
.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.prop.inv(.05, 10, 87)
#>   Estimate         SE         LL        UL
#>  0.1149425 0.03389574 0.05651668 0.1893855

# Should return:
#  Estimate         SE         LL        UL
# 0.1149425 0.03389574 0.05651668 0.1893855

```
